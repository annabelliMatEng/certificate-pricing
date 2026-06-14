function [Libor_leg, Coupon_leg, X_upfront] = UpfrontPricingLEWIS(params_NIG, DiscountFactorLibor,...
    spread,delta_libor,DiscountFactorCoupon,delta_coupon)

% LEWIS FORMULA: Calculates the Net Present Value (NPV) of a swap-linked certificate.
%
% INPUTS:
%   - params_NIG: Structure for the NIG model containing:
%       .kappa: Kurtosis/variance of the subordinator.
%       .eta: Asymmetry (skewness) parameter.
%       .sigma: Volatility parameter of the NIG process.
%       .moneyness: Log-moneyness defined as log(Forward/Strike).
%       .time_to_maturity: Time to the coupon reset date (T).
%   - DiscountFactorLibor: Vector of discount factors for the Libor leg (Party A).
%   - spread: The fixed spread added to Euribor 3m.
%   - delta_libor: Year fractions for the quarterly Libor payment periods (ACT/360).
%   - DiscountFactorCoupon: Vector of discount factors for the annual coupon payments.
%   - delta_coupon: Year fractions for the annual coupon periods (30/360).
%
% OUTPUTS:
%   - Libor_leg: NPV of the floating payments (Euribor + Spread) considering early redemption.
%   - Coupon_leg: NPV of the conditional coupons paid by Bank XX.
%   - X_upfront: The balancing payment (X%) received by Bank XX at Start Date.



kappa = params_NIG.kappa;
eta = params_NIG.eta;
sigma = params_NIG.sigma;
moneyness = params_NIG.moneyness;
t = params_NIG.time_to_maturity;

L_base = @(u) exp( (1 - sqrt(1 + 2 * kappa * u)) / kappa );
omega  = -log(L_base( eta * sigma^2 * t ));
u_arg  = @(xi) (0.5 * xi.^2 + 1i * (0.5 + eta) * xi) * sigma^2 * t;
phi    = @(xi) exp(1i * xi * omega) .* L_base(u_arg(xi));

integrand = @(u) real((exp(1i.*u.*moneyness)).*phi(u)./(1i.*u));
p = (1/2) + (1/pi)*quadgk(integrand,0,inf);
digital_put_prob = 1-p;

    



Libor_payments = (DiscountFactorLibor(1:end-1) - DiscountFactorLibor(2:end)) + ...
    spread*delta_libor.*DiscountFactorLibor(2:end);

NPV_libor_1 = sum(Libor_payments(1:4));
NPV_libor_2 = sum(Libor_payments(5:8)) * (1-digital_put_prob);
Libor_leg = NPV_libor_1+NPV_libor_2;

Coupon_leg = (0.06)*(DiscountFactorCoupon(1))*delta_coupon(1)*digital_put_prob + ...
    (0.02)*(DiscountFactorCoupon(2))*delta_coupon(2)*(1-digital_put_prob);

X_upfront = Libor_leg - Coupon_leg;


%% --- DISPLAY ALL INTERMEDIATE AND FINAL VARIABLES ---
fprintf('\n================================================\n');
fprintf('         LEWIS PRODUCT VALUATION         \n');
fprintf('================================================\n');

fprintf('Model Used:        NIG (Normal Inverse Gaussian)\n');
fprintf('Digital Put Prob:  %.4f%%\n', digital_put_prob * 100);
fprintf('NPV Libor Year 1:  %.6f\n', NPV_libor_1);
fprintf('NPV Libor Year 2:  %.6f (Adjusted for Early Red.)\n', NPV_libor_2);
fprintf('TOTAL LIBOR LEG:   %.6f\n', Libor_leg);
fprintf('TOTAL COUPON LEG:  %.6f\n', Coupon_leg);
fprintf('------------------------------------------------\n');
fprintf('X%% UPFRONT:        %.4f%%\n', X_upfront * 100);
fprintf('================================================\n\n');

end