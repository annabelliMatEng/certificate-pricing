function [Upfront_X, MC_StdErr] = UpfrontPricingGamma(numSim, Fwd, DF_coupon, DF_libor, ...
                                      delta_NIG, delta_libor,delta_coupon, params, Strike, ...
                                      Spread, BoolControl )
% UPFRONTPRICING  Computes the upfront X paid by Party B at t0 such that NPV = 0.
%
% INPUTS:
%   numSim      - Number of Monte Carlo simulation paths
%   Fwd         - [F(t0,t1), F(t0,t2)]  Forward prices at year 1 and year 2
%   DF_coupon   - [B(t0,t1), B(t0,t2)]  Discount factors at coupon dates
%   DF_libor    - [B(t0,t0), ..., B(t0,t8)]  Discount factors at quarterly dates
%   delta_NIG   - [dt1, dt2]  ACT/365 time steps for NIG simulation
%   delta_libor - [dQ1,...,dQ8]  ACT/360 fractions for floating leg
%   params      - Struct with NIG params: sigma, kappa, eta
%   Strike      - Barrier level K (3200 from Annex 1)
%   Spread      - Floating spread over Euribor (0.013 = 1.30%)
%        - Contract   in EUR
%
% OUTPUTS:
%   Upfront_X - Value of X in EUR (divide by   for percentage)
%   MC_StdErr - Monte Carlo standard error on PV_CouponLeg

    %% 0. Unpack NIG parameters
    sigma = params.sigma;
    kappa = params.kappa;
    eta   = params.eta;
    t = params.time_to_maturity;

    % ------------------------------------------------------------------ %
    %% 1. NIG Martingale Correction (omega)
    % ------------------------------------------------------------------ %
    % We compute omega once (it is time-homogeneous for NIG).
    L_sub = @(u) (-t/kappa) * log(1+kappa*u);  % QUESTO VA CAMBIATO
    omega = -L_sub(eta * sigma^2);   

    % ------------------------------------------------------------------ %
    %% 2. Monte Carlo: Simulate S at t1 (Year 1 reset date)
    % ------------------------------------------------------------------ %
    dt1       = delta_NIG(1);
    a = dt1/kappa;
    b = kappa;

    G1  = random('Gamma', a, b, [numSim, 1]);
    Z1  = randn(numSim, 1);

    dX1 = omega * dt1 + (-(0.5 + eta) * sigma^2) .* G1 + sigma * sqrt(G1) .* Z1;
    S1  = Fwd(1) * exp(dX1);

    % ------------------------------------------------------------------ %
    %% 3. Estimate Digital Probabilities at t1 via Monte Carlo
    % ------------------------------------------------------------------ %

    digital_mask_t1 = (S1 <= Strike);           % logical vector: 1 if triggered

    % ------------------------------------------------------------------ %
    %% 4. Early Redemption Logic (Autocall / Trigger)
    % ------------------------------------------------------------------ %

    is_auto  = digital_mask_t1;    % boolean mask for autocalled paths
    survived = ~is_auto;           % boolean mask for surviving paths

    % ------------------------------------------------------------------ %
    %% 5. Monte Carlo: Simulate S at t2 (Year 2 reset date) - Survived paths only
    % ------------------------------------------------------------------ %
    %
    % The starting point uses the FORWARD MEASURE chaining:
    %   S2 = S1 * (F(t0,t2) / F(t0,t1)) * exp(dX2)

    dt2     = delta_NIG(2);
    a = dt2/kappa;
    b = kappa;

    % Simulate for ALL paths (vectorized); we will only use survived ones
    G2  = random('Gamma', a, b, [numSim, 1]);
    Z2  = randn(numSim, 1);

    dX2 = omega * dt2 + (-(0.5 + eta) * sigma^2) .* G2 + sigma * sqrt(G2) .* Z2;
    S2  = S1 .* (Fwd(2) / Fwd(1)) .* exp(dX2);  % chained from S1
    % ------------------------------------------------------------------ %
    %% 5b. Define Logical Masks for Path Outcomes
    % ------------------------------------------------------------------ %
    % We already have 'is_auto' (S1 <= K) and 'survived' (S1 > K) from Year 1.
    % Now we define the Year 2 outcomes.
    
    idx_digital_t2 = survived & (S2 <= Strike); % Survives t1, hits barrier at t2
    idx_vanilla_t2 = survived & (S2 >  Strike); % Survives t1, misses barrier at t2
    
    % Probabilities for diagnostics
    prob_auto_t1   = mean(is_auto);
    prob_digital_t2_conditional = sum(idx_digital_t2) / sum(survived);

    % ------------------------------------------------------------------ %
    %% 6. Path-by-Path Present Value of the Floating Leg (Party A)
    % ------------------------------------------------------------------ %
    P_start = DF_libor(1:end-1);   % B(t0, t_{i-1})
    P_end   = DF_libor(2:end);     % B(t0, t_i)
    
    % All 8 quarterly discounted cash flows (assuming Notional is defined, e.g., 100e6)
    CF_float_all =  1 * ((P_start - P_end) + Spread .* delta_libor .* P_end);
    
    PV_float_1Y = sum(CF_float_all(1:4)); % 4 quarters if autocalled at t1
    PV_float_2Y = sum(CF_float_all(1:8)); % 8 quarters if contract runs to t2
    PV_float_3Y = sum(CF_float_all);
    
    % Allocate floating payoffs to each path
    PV_paths_floating = zeros(numSim, 1);
    PV_paths_floating(is_auto)  = PV_float_1Y;
    PV_paths_floating(survived) = PV_float_2Y;

    if BoolControl
        PV_paths_floating(is_auto)  = PV_float_1Y;
        PV_paths_floating(idx_digital_t2) = PV_float_2Y;
        PV_paths_floating(idx_vanilla_t2) = PV_float_3Y;
    end
    
    % ------------------------------------------------------------------ %
    %% 7. Path-by-Path Present Value of the Coupon Leg (Party B)
    % ------------------------------------------------------------------ %
    PV_paths_coupon = zeros(numSim, 1);
    
    % SCENARIO 1: Autocall at t1
    % Pays 6% at t1, contract ends.
    PV_paths_coupon(is_auto) =  1 * 0.06 * delta_coupon(1) * DF_coupon(1);
    
    % SCENARIOS 2 & 3: Survival to t2
    if BoolControl 
        % Point (d): Digital active at t2
        
        % Scenario 2: Digital triggered at t2
        % Pays 2% (accrued at t1) + 2% fixed (at t2) + 6% digital (at t2)
        PV_paths_coupon(idx_digital_t2) =  1 * (...
            0.00 * delta_coupon(1) * DF_coupon(1) + ...
            0.06 * delta_coupon(2) * DF_coupon(2) ); 
        
        % Scenario 3: Vanilla survival at t2
        % Pays 2% (accrued at t1) + 2% fixed (at t2)
        PV_paths_coupon(idx_vanilla_t2) =  1 * (...
            0.00 * delta_coupon(1) * DF_coupon(1) + ...
            0.00 * delta_coupon(2) * DF_coupon(2) + ...
            0.02 * delta_coupon(3) * DF_coupon(3));
    else
        % Point (a): No digital at t2, just standard fixed coupon
        % Pays 2% (accrued at t1) + 2% fixed (at t2) for ALL survived paths
        PV_paths_coupon(survived) =  1 * (...
            0.00 * delta_coupon(1) * DF_coupon(1) + ...
            0.02 * delta_coupon(2) * DF_coupon(2) );
    end

    % ------------------------------------------------------------------ %
    %% 8. Upfront Premium X & Monte Carlo Standard Error
    % ------------------------------------------------------------------ %
    % The Upfront is the net difference between the legs for each path
    % NPV = Upfront_X + PV_Floating - PV_Coupon = 0 
    % Upfront_X = PV_Coupon - PV_Floating
    
    X_paths = PV_paths_floating - PV_paths_coupon;
    
    % Final expected values
    Upfront_X      = mean(X_paths);
    PV_CouponLeg   = mean(PV_paths_coupon);
    PV_FloatingLeg = mean(PV_paths_floating);
    
    % Standard Error calculated on the final net Upfront variable
    MC_StdErr = std(X_paths) / sqrt(numSim);

    % ------------------------------------------------------------------ %
    %% 9. Display intermediate results for transparency
    % ------------------------------------------------------------------ %
    fprintf('--- Monte Carlo Diagnostics ---\n');
    fprintf('  P^Q(S1 <= K=3200)              : %.4f%%\n', prob_auto_t1 * 100);
    fprintf('  P^Q(S2 <= K | survived to t1)  : %.4f%%\n', prob_digital_t2_conditional * 100);
    fprintf('  PV Coupon Leg (Avg)            : %.2f EUR\n', PV_CouponLeg);
    fprintf('  PV Floating Leg (Avg)          : %.2f EUR\n', PV_FloatingLeg);
    fprintf('  UPFRONT PREMIUM (X)            : %.2f EUR\n', Upfront_X);
    fprintf('  Upfront as %% of  1      : %.4f%%\n', (Upfront_X /  1) * 100);
    fprintf('  Monte Carlo Standard Error     : %.2f EUR\n', MC_StdErr);
    fprintf('-------------------------------\n');


end