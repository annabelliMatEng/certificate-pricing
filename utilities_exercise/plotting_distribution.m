%% Comparison of Subordinators: Inverse Gaussian (NIG) vs Gamma (VG)
clear; clc;

% 1. Parameters
dt1 = 1;          % Time step (1 year)
kappa_NIG = 0.4525;
kappa_VG =0.6674 ;
numSim = 1e6;     % High number of simulations for smooth density

% 2. Generate Random Variables
% Inverse Gaussian (NIG)
mu_IG1  = dt1;
lam_IG1 = dt1^2 / kappa_NIG;
G_IG = random('InverseGaussian', mu_IG1, lam_IG1, [numSim, 1]);

% Gamma (VG)
a = dt1 / kappa_VG;
b = kappa_VG;
G_Gamma = random('Gamma', a, b, [numSim, 1]);

% 3. Plotting
figure('Color', 'w', 'Name', 'Subordinator Comparison');
hold on;

% Use ksdensity for smooth PDF estimation
[f_ig, x_ig] = ksdensity(G_IG);
[f_ga, x_ga] = ksdensity(G_Gamma);

plot(x_ig, f_ig, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]); % Orange
plot(x_ga, f_ga, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);     % Blue

% Graphics refinements
grid on;
xlim([0, 4]); % Focus on the main mass
title('Density of Subordinators: IG vs Gamma');
xlabel('Time Increment (G)');
ylabel('Probability Density');
legend('Inverse Gaussian (NIG)', 'Gamma (VG)', 'Location', 'NorthEast');

