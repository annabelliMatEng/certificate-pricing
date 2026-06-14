function Upfront_X = UpfrontPricingBS(Fwd, DF_coupon, DF_libor, ...
                                       delta_coupon, delta_libor,delta_timesteps, strike, ...
                                       spread, vol, coupon_payment, varargin)
% UPFRONTPRICINGBS  Computes the upfront X (as fraction of notional) with
%                   optional smile/digital-risk correction (slope impact)
% =========================================================================
% INPUTS:
%   Fwd           - [F(t0,t1), F(t0,t2)]  Forward prices
%   DF_coupon     - [B(t0,t1), B(t0,t2)]  Discount factors at annual dates
%   DF_libor      - [B(t0,t0),...,B(t0,tQ8)]  9 quarterly discount factors
%   delta_coupon  - [dt1, dt2]  30/360 day count fractions
%   delta_libor   - [dQ1,...,dQ8]  ACT/360 day count fractions
%   strike        - Barrier level K 
%   spread        - Spread over Euribor (0.013)
%   vol           - Implied vol at (K, T1)  [scalar]
%   coupon_payment- First year and second year coupons
%
% OPTIONAL name-value pairs:
%   'digital_risk'  true/false   Enable smile slope correction (default: false)
%   'vol_up'        scalar       Implied vol at K+eps  (required if digital_risk=true)
%   'vol_down'      scalar       Implied vol at K-eps  (required if digital_risk=true)
%   'eps_strike'    scalar       Strike bump in index points (default: 1.0)
%
% OUTPUT:
%   Upfront_X  - Upfront premium as fraction of notional

    %% 0. Parse optional arguments
    ip = inputParser();
    addParameter(ip, 'digital_risk', false, @islogical);
    addParameter(ip, 'vol_up',       NaN,   @isnumeric);
    addParameter(ip, 'vol_down',     NaN,   @isnumeric);
    addParameter(ip, 'eps_strike',   1.0,   @isnumeric);
    parse(ip, varargin{:});

    digital_risk = ip.Results.digital_risk;
    vol_up       = ip.Results.vol_up;
    vol_down     = ip.Results.vol_down;
    eps_strike   = ip.Results.eps_strike;

    % Validate: if digital_risk=true, bumped vols must be provided
    if digital_risk
        assert(~isnan(vol_up),   'vol_up required when digital_risk=true');
        assert(~isnan(vol_down), 'vol_down required when digital_risk=true');
    end

    %% 1. Force all inputs to row vectors (robust against orientation)
    Fwd            = Fwd(:)';
    DF_coupon      = DF_coupon(:)';
    DF_libor       = DF_libor(:)';
    delta_coupon   = delta_coupon(:)';
    delta_libor    = delta_libor(:)';
    coupon_payment = coupon_payment(:)';

    %% 2. Split quarterly schedule
    delta_libor_1y = delta_libor(1:4);   % ACT/360 fracs, quarters 1-4
    delta_libor_2y = delta_libor(5:8);   % ACT/360 fracs, quarters 5-8
    DF_libor_1Y    = DF_libor(2:5);      % B(t0,tQ1)..B(t0,tQ4)
    DF_libor_2Y    = DF_libor(6:9);      % B(t0,tQ5)..B(t0,tQ8)

    %% 3. Digital put probability at t1:  p = P^Q(S_t1 <= K)
    T1 = delta_timesteps(1);
    F1 = Fwd(1);

    % Black-76 d1 and d2
    d1 = (log(F1 / strike) + 0.5 * vol^2 * T1) / (vol * sqrt(T1));
    d2 = d1 - vol * sqrt(T1);

    % Plain Black term: N(-d2)
    p_black = normcdf(-d2);

    if ~digital_risk
        % Ignore smile slope
        p                = p_black;
        dsigma_dK        = 0;
        vega_black       = 0;
        slope_correction = 0;

    else
        % Smile slope: central finite difference
        dsigma_dK = (vol_up - vol_down) / (2 * eps_strike);

        % Black-76 Vega (undiscounted)
        vega_black = F1 * sqrt(T1) * normpdf(d1);

        % Slope impact correction
        slope_correction = dsigma_dK * vega_black;

        % Smile-adjusted probability
        p = p_black + slope_correction;

        % Clamp to [0,1] for safety
        p = max(0, min(1, p));
    end

    %% 4. PV Floating Leg
    % Year 1: always paid (Euribor telescopes: 1 - B(t0,t1))
    NPV_floating_1y = (1 - DF_coupon(1)) + spread * dot(delta_libor_1y, DF_libor_1Y);

    % Year 2: paid only if NOT autocalled -> weight (1-p)
    NPV_floating_2y = (DF_coupon(1) - DF_coupon(2) + ...
                       spread * dot(delta_libor_2y, DF_libor_2Y)) * (1 - p);

    %% 5. PV Coupon Leg
    % coupon_payment are fixed period rates (not rate*delta)
    NPV_coupon = coupon_payment(1) * DF_coupon(1) * p       + ...
                 coupon_payment(2) * DF_coupon(2) * (1 - p);

    %% 6. Upfront X (fraction of notional)
    Upfront_X = NPV_floating_1y + NPV_floating_2y - NPV_coupon;
end