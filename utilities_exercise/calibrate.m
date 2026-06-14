function [sigma, k,eta]=calibrate(F0,B,strikes,surface,dt,alpha)

r_mat = -log(B) / dt;
moneyness = log(F0 ./ strikes);
%% Characteristic function of the log-return
% Laplace exponent of the subordinator
if alpha>0
    lnL = @(omega, k, sigma) ...
   (dt/k) * ((1-alpha)/alpha) * (1 - (1 + (omega * k * sigma^2) / (1-alpha)).^alpha);
elseif alpha==0
    lnL = @(omega, k, sigma) -(dt/k)*log(1+k*omega*sigma^2);
else 
  error('Error: alpha can not be negative (alpha = %f).', alpha);
end

% Characteristic function with martingale correction
char_func = @(csi, sigma, kappa, eta) ...
    exp(-1i * csi .* real(lnL(eta, kappa, sigma))) .* ...
    exp(lnL((csi.^2 + 1i*(1 + 2*eta)*csi) / 2, kappa, sigma));

%% Market reference prices via Black model
[BlackPrices, ~] = blkprice(F0, strikes, r_mat, dt, surface);

%% FFT parameters
dz = 0.25 * 1e-2;   % moneyness grid step
M  = 15;             % FFT grid size exponent (N = 2^M)

%% Calibration: global least-squares with constant weights
objective = @(P) sum(abs( ...
    FFT_obj(moneyness, @(u) char_func(u, P(1), P(2), P(3)), B, F0, M, dz, 'dz') ...
    - BlackPrices).^2);

P0 = [0.20, 1.0, 0.5];      % initial guess: [sigma, kappa, eta]
lb = [0.01, 0.01, -50];     % lower bounds
ub = [1.00, 1.00,  50];     % upper bounds

% Nonlinear feasibility constraint on eta
condition = @(p) conditionOnEta(p, alpha);

[P_opt, ~] = fmincon(objective, P0, [], [], [], [], lb, ub, condition);

% Display calibrated parameters
fprintf('\n--- Calibration Results ---\n');
fprintf('sigma : %.4f\n', P_opt(1));
fprintf('kappa : %.4f\n', P_opt(2));
fprintf('eta   : %.4f\n', P_opt(3));
fprintf('\n---------------------------\n');

sigma=P_opt(1);
k=P_opt(2);
eta=P_opt(3);
end



