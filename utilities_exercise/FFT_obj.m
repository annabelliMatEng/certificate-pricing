function C_FFT = FFT_obj(x_target, cf, discount_factor, F0, M, value, type)
% European Call pricing via Lewis formula using FFT
%
% Inputs:
%   x_target        - Vector of log-moneyness target points [ln(F0/K)]
%   cf              - Function handle for the characteristic function @(u)
%   discount_factor - Discount factor B(t0, T)
%   F0              - Forward price
%   M               - Degree of freedom for grid size (N = 2^M)
%   value           - the numerical value of the parameter
%   type            - a string, either 'dz' or 'x1'
%
% Output:
%   C_FFT           - Vector of call prices corresponding to x_target

% 1. Parameters & Relations (Your exact notation)
N = 2^M;

% Grid Logic based on selection
if strcmpi(type, 'dz')
    % User selected dz
    dz = value;
    dx = (2*pi) / (N * dz);
    x1 = -dx * (N-1) / 2;
elseif strcmpi(type, 'x1')
    % User selected x1 (integration extremum)
    x1 = value;
    dx = -2 * x1 / (N-1);
    dz = (2*pi) / (N * dx);
else
    error('Invalid type. Choose either ''dz'' or ''x1''');
end

% Sanity check
[x_target, idx] = sort(x_target);

% Symmetry of the grids 
x1 = -dx * (N-1) / 2;
xN = -x1;
x_grid = x1 : dx : xN;

z1 = -dz * (N-1) / 2;
zN = -z1;   
z_grid = z1 : dz : zN;

% 2. Definition of Integrand
f_x = (1/(2*pi)) * cf(-x_grid - 1i/2) ./ (x_grid.^2 + 1/4);

% 3. FFT
j_minus_1 = 0:N-1;
input_fft = f_x .* exp(-1i * z1 * dx * j_minus_1);

% Fast Fourier Transform
Y = fft(input_fft);
    
% Integral Reconstruction via Prefactor
integral_values = dx * exp(-1i * x1 * z_grid) .* Y;

% 4. Interpolation to Moneyness
integral_interp = interp1(z_grid, real(integral_values), x_target, 'spline');
prices = discount_factor * F0 * (1 - exp(-x_target/2) .* integral_interp);

% Prices in the original order
C_FFT(idx) = prices;
C_FFT = C_FFT';

end