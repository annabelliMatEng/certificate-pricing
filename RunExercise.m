% Exercise 1 - Certificate Pricing

clc; clear; close all;

addpath(fullfile(pwd, 'utilities'))
addpath(fullfile(pwd, 'utilities_ex1'))

%% 1. Data Loading & Bootstrap
% Load market data and volatility surface (Assignment 5 Data)
temp = load('eurostoxx_Poli.mat');
S = temp.cSelect; 

% Market Curve Bootstrap (Assignment 2 Data)
formatDate = 'dd/mm/yyyy'; 
if ispc
    fprintf('Operating System: Windows detected. Loading data...\n');
    [datesSet, ratesSet] = readExcelData_windows('MktData_CurveBootstrap.xls', formatDate);
    
elseif ismac
    fprintf('Operating System: macOS detected. Loading data...\n');
    [datesSet, ratesSet] = readExcelData_mac('MktData_CurveBootstrap.xls', formatDate);
    
else
    error('Unsupported Operating System. Please use Windows or macOS.');
end
[BootStrapDates, ~, BootStrapZeroRates] = bootstrap(datesSet, ratesSet);

%% 2. Parameters & Financial Setup
% Market Data extraction
strikes   = S.strikes;
surface   = S.surface;      % Implied Volatility vector
reference = S.reference;    % Spot Price (S0)
q         = S.dividends;    % Continuous Dividend Yield

% Contract Details
coupon_payment=[0.06; 0.02];
Notional = 100e6; % 100 Million EUR

% Time Setup (using datetime for robust calendar math)
startdate_dt = datetime('19-Feb-2008');
enddate_dt   = startdate_dt + calyears(3); % 3-year maturity

% Generate schedule arrays (includes t0)
libor_dates_dt  = startdate_dt : calmonths(3) : enddate_dt; % Quarterly
coupon_dates_dt = startdate_dt : calyears(1) : enddate_dt;  % Annual

% Convert datetime objects to standard MATLAB convertDates for built-in functions
startdate    = ConvertDates(startdate_dt);
enddate      = ConvertDates(enddate_dt);
libor_dates  = ConvertDates(libor_dates_dt); % first one is t0
coupon_dates = ConvertDates(coupon_dates_dt); % first one is t0

% Reset dates are 2 business days prior to the payment dates
% NaT ensures standard weekends (Sat/Sun) are skipped.
reset_dates = busdate(coupon_dates(2:end), -1, NaT); 
reset_dates = busdate(reset_dates, -1, NaT); 

%% 3. Day Count Fractions (yearfrac)
% delta(t_i, t_i+1) for the Coupon Leg -> 30/360 EURO convention
delta_coupon = yearfrac(coupon_dates(1:end-1), coupon_dates(2:end), 6); 

% delta(t_i, t_i+1) for the Euribor Leg -> ACT/360 convention
delta_libor = yearfrac(libor_dates(1:end-1), libor_dates(2:end), 2); 

% delta(t_i, t_i+1) for NIG Simulation -> ACT/365 convention 
delta_timesteps = yearfrac(reset_dates(1:end-1), reset_dates(2:end), 3); 

%% 4. Market Variable Calculations
% Discount Factors from Bootstrapped Zero Rates
DiscountFactor_coupon = fromdatetodiscount(startdate, BootStrapDates, BootStrapZeroRates, coupon_dates);
DiscountFactor_coupon = DiscountFactor_coupon(2:end); % [B(t0,1y), B(t0,2y)]
DiscountFactor_libor  = fromdatetodiscount(startdate, BootStrapDates, BootStrapZeroRates, libor_dates);


% Forward Price Calculation
Forward_1Y = (reference / DiscountFactor_coupon(1)) * exp(-q * delta_timesteps(1)); 
Forward_2y = (reference / DiscountFactor_coupon(2)) * exp(-q * (delta_timesteps(1)+delta_timesteps(2) ) );
Fwd = [Forward_1Y, Forward_2y];

%% MonteCarlo Valuation NIG

Spread = 1.30*1e-2;
Strike = 3200;

[params_NIG.sigma, params_NIG.kappa, params_NIG.eta]=...
    calibrate(Forward_1Y,DiscountFactor_coupon(1), strikes, surface, delta_timesteps(1), 1/2);
params_NIG.moneyness = log(Forward_1Y / Strike);
params_NIG.time_to_maturity = delta_timesteps(1);

numSim = 1e7;
BoolControl = 0; % 0 for 2 years maturity , 1 for 3 years maturity
model = 'NIG';


[Upfront_X_NIG, MC_StdErr_NIG] = UpfrontPricingMC(numSim, Fwd, DiscountFactor_coupon, DiscountFactor_libor, ...
                                      delta_timesteps, delta_libor,delta_coupon, params_NIG, Strike, ...
                                      Spread, BoolControl, model);

%% Monte Carlo VG


[params_VG.sigma, params_VG.kappa, params_VG.eta]=...
    calibrate(Forward_1Y,DiscountFactor_coupon(1), strikes, surface, delta_timesteps(1), 0);
params_VG.time_to_maturity = delta_timesteps(1);
params_VG.moneyness = log(Forward_1Y / Strike);

numSim = 1e7;
BoolControl = 0; % 0 for 2 years maturity , 1 for 3 years maturity
model = 'VG';

[Upfront_X_VG, MC_StdErr_VG] = UpfrontPricingMC(numSim, Fwd, DiscountFactor_coupon, DiscountFactor_libor, ...
                                      delta_timesteps, delta_libor,delta_coupon, params_VG, Strike, ...
                                      Spread, BoolControl,model);

%% Lewis Valuation

[Libor_leg_lewis, Coupon_leg_lewis, X_upfront_lewis]= UpfrontPricingLEWIS(params_NIG, DiscountFactor_libor,...
    Spread,delta_libor,DiscountFactor_coupon,delta_coupon);


%% Black Valuation

eps_strike = 1.0;             % 1 index point bump
vol_K      = interp1(strikes, surface, Strike,               'spline');
vol_K_up   = interp1(strikes, surface, Strike + eps_strike,  'spline');
vol_K_down = interp1(strikes, surface, Strike - eps_strike,  'spline');
vol = surface(4);

Upfront_X_black_smile = UpfrontPricingBS(Fwd, DiscountFactor_coupon, DiscountFactor_libor, ...
                                       delta_coupon, delta_libor,delta_timesteps, Strike, ...
                                       Spread, vol, coupon_payment, ...
                                         'digital_risk', true,      ...
                                         'vol_up',       vol_K_up,  ...
                                         'vol_down',     vol_K_down, ...
                                         'eps_strike',   eps_strike);

Upfront_X_black = UpfrontPricingBS(Fwd, DiscountFactor_coupon, DiscountFactor_libor, ...
                                       delta_coupon, delta_libor,delta_timesteps, Strike, ...
                                       Spread, vol, coupon_payment);

fprintf('\n================================================\n');
fprintf('  Black Model Comparison Summary\n');
fprintf('================================================\n');
fprintf('  %-35s : %8.4f%%\n', 'BS plain ',   Upfront_X_black          * 100);
fprintf('  %-35s : %8.4f%%\n', 'BS smile-adjusted',   Upfront_X_black_smile    * 100);
fprintf('  %-35s : %8.4f%%\n', 'Diff: BS_smile - BS_plain',    (Upfront_X_black_smile - Upfront_X_black) * 100);
fprintf('================================================\n');

%% Upfront Valuation in case of 3 years maturity

BoolControl = 1;

model = 'NIG';
[Upfront_X_NIG_3y, MC_StdErr_NIG_3y] = UpfrontPricingMC(numSim, Fwd, DiscountFactor_coupon, DiscountFactor_libor, ...
                                      delta_timesteps, delta_libor,delta_coupon, params_NIG, Strike, ...
                                      Spread, BoolControl, model);
model = 'VG';
[Upfront_X_VG_3y, MC_StdErr_VG_3y] = UpfrontPricingMC(numSim, Fwd, DiscountFactor_coupon, DiscountFactor_libor, ...
                                      delta_timesteps, delta_libor,delta_coupon, params_VG, Strike, ...
                                      Spread, BoolControl,model);
