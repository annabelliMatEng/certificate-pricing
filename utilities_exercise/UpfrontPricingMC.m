function [Upfront_X, MC_StdErr] = UpfrontPricingMC(numSim, Fwd, DF_coupon, DF_libor, ...
                                      delta_NIG, delta_libor,delta_coupon, params, Strike, ...
                                      Spread, BoolControl, model )
% UPFRONTPRICING  Computes the upfront X paid by Party B at t0 such that NPV = 0.
%
% INPUTS:
%   numSim      - Number of Monte Carlo simulation paths
%   Fwd         - Forward prices for years of interest
%   DF_coupon   - Discount factors at coupon dates
%   DF_libor    - Discount factors at quarterly dates
%   delta_NIG   - [ACT/365 time steps for underlying simulation
%   delta_libor - ACT/360 fractions for floating leg
%   params      - Struct with model params: sigma, kappa, eta
%   Strike      - Barrier level K 
%   Spread      - Floating spread over Euribor 
%   BoolControl - 0 = 2y expiry ; 1 = 3y expiry
%   model       - 'NIG' or 'VG'
%
% OUTPUTS:
%   Upfront_X - Value of X in pct 
%   MC_StdErr - Monte Carlo standard error

    %% 0. Unpack NIG parameters
    rng(8); % setting the seed
    sigma = params.sigma;
    kappa = params.kappa;
    eta   = params.eta;

    % ------------------------------------------------------------------ %
    %% 1. omega computation
    % ------------------------------------------------------------------ %
    if strcmpi(model, 'NIG')
        L_sub = @(u) (1 - sqrt(1 + 2 * kappa * u)) / kappa;  % Laplace exp of IG
        omega = -L_sub(eta * sigma^2);
    end
    if strcmpi(model, 'VG')
        L_sub = @(u) (-1/kappa) * log(1+kappa*u); % Laplace exp of VG
        omega = -L_sub(eta * sigma^2); 
    end


    % ------------------------------------------------------------------ %
    %% 2. Monte Carlo: Simulate S at t1 (Year 1 reset date)
    % ------------------------------------------------------------------ %
    if strcmpi(model, 'NIG')
        dt1       = delta_NIG(1);
        mu_IG1    = dt1;
        lam_IG1   = dt1^2 / kappa;
    
        G1  = random('InverseGaussian', mu_IG1, lam_IG1, [numSim, 1]);
    end
    if strcmpi(model, 'VG')
        dt1       = delta_NIG(1);
        a = dt1/kappa;
        b = kappa;

        G1  = random('Gamma', a, b, [numSim, 1]);
    end

    Z1  = randn(numSim, 1);

    dX1 = omega + (-(0.5 + eta) * sigma^2) .* G1*dt1 + sigma * sqrt(G1*dt1) .* Z1;
    S1  = Fwd(1) * exp(dX1);

    % ------------------------------------------------------------------ %
    %% 3. Estimate Digital Probabilities at t1 via Monte Carlo
    % ------------------------------------------------------------------ %

    digital_mask_t1 = (S1 <= Strike);       

    % ------------------------------------------------------------------ %
    %% 4. Early Redemption Logic (Autocall / Trigger)
    % ------------------------------------------------------------------ %

    is_auto  = digital_mask_t1;    % boolean mask for autocalled paths
    survived = ~is_auto;           % boolean mask for surviving paths

    % ------------------------------------------------------------------ %
    %% 5. Monte Carlo: Simulate S at t2 (Year 2 reset date) - Survived paths only
    % ------------------------------------------------------------------ %
    %
    if strcmpi(model, 'NIG')
        dt2     = delta_NIG(2);
        mu_IG2  = dt2;
        lam_IG2 = dt2^2 / kappa;
    
        % Simulate for ALL paths (vectorized); we will only use survived ones
        G2  = random('InverseGaussian', mu_IG2, lam_IG2, [numSim, 1]);
    end
    if strcmpi(model, 'VG')
        dt2     = delta_NIG(2);
        a = dt2/kappa;
        b = kappa;
    
        % Simulate for ALL paths (vectorized); we will only use survived ones
        G2  = random('Gamma', a, b, [numSim, 1]);
    end

    Z2  = randn(numSim, 1);

    dX2 = omega + (-(0.5 + eta) * sigma^2) .* G2*dt2 + sigma * sqrt(G2*dt2) .* Z2;
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
    
    % All quarterly discounted cash flows 
    CF_float_all =  1 * ((P_start - P_end) + Spread .* delta_libor .* P_end);
    
    PV_float_1Y = sum(CF_float_all(1:4)); % 4 quarters if autocalled at t1
    PV_float_2Y = sum(CF_float_all(1:8)); % 8 quarters if contract runs to t2
    PV_float_3Y = sum(CF_float_all); % all quarters if contract runs up to t3
    
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
        
        % Scenario 2: Digital triggered at t2
        PV_paths_coupon(idx_digital_t2) =  1 * (...
            0.00 * delta_coupon(1) * DF_coupon(1) + ...
            0.06 * delta_coupon(2) * DF_coupon(2) ); 
        
        % Scenario 3: Vanilla survival at t2
        PV_paths_coupon(idx_vanilla_t2) =  1 * (...
            0.00 * delta_coupon(1) * DF_coupon(1) + ...
            0.00 * delta_coupon(2) * DF_coupon(2) + ...
            0.02 * delta_coupon(3) * DF_coupon(3));
    else
        % Point (a): No digital at t2, just standard fixed coupon
        PV_paths_coupon(survived) =  1 * (...
            0.00 * delta_coupon(1) * DF_coupon(1) + ...
            0.02 * delta_coupon(2) * DF_coupon(2) );
    end

    % ------------------------------------------------------------------ %
    %% 8. Upfront Premium X & Monte Carlo Standard Error
    % ------------------------------------------------------------------ %
    
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
    if  BoolControl
        fprintf('--- Monte Carlo Diagnostics %s model 3 years maturity---\n', model);
        fprintf('  P^Q(S1 <= K=3200)              : %.4f%%\n', prob_auto_t1 * 100);
        fprintf('  P^Q(S2 <= K | survived to t1)  : %.4f%%\n', prob_digital_t2_conditional * 100);
    else
        fprintf('--- Monte Carlo Diagnostics %s model 2 years maturity ---\n', model);
        fprintf('  P^Q(S1 <= K=3200)              : %.4f%%\n', prob_auto_t1 * 100);
    end
    fprintf('  PV Coupon Leg (Avg)            : %.4f EUR\n', PV_CouponLeg);
    fprintf('  PV Floating Leg (Avg)          : %.4f EUR\n', PV_FloatingLeg);
    fprintf('  UPFRONT PREMIUM (X)            : %.4f EUR\n', Upfront_X);
    fprintf('  Upfront as %% of  1      : %.4f%%\n', (Upfront_X /  1) * 100);
    fprintf('  Monte Carlo Standard Error     : %.6f EUR\n', MC_StdErr);
    fprintf('-------------------------------\n');


end
