% Monte Carlo Uncertainty Quantification for Japan's Population Projection
% Governing Equation: dP/dt = a*P^2 + b*P
% Numerical Method: Euler's Method (h = 0.125)
clear; clc; close all;
%% 0. Computing the k vector from cvs file
x = [ 2017  2020   2023   2026]; %year

P = [127.065 126.305 124.371 122.428]; %population


for i = 1: length(P)-1
k(i) =(1/3)* log(P(i+1)/P(i));
end
k

%% 0. Solving Matrix for Two linear Splines
k=[ -0.001999716565635  -0.005143527565717  -0.005248644023522];
  
x=[126.305 124.371 122.428];

M = [x(1)^1 x(1)^0 0 0 ; % row 1
%-------------------------------
     x(2)^1 x(2)^0 0 0 ; %row 2 
%-------------------------------
     0 0 x(2)^1 x(2)^0 ; % row 3
%-------------------------------
     0 0 x(3)^1 x(3)^0 ]; % row 4

      
 b= [k(1); k(2); k(2); k(3)];
 
 y = inv(M)*b % COEFFICIENTS


%% 1. Nominal Parameters and Time Discretization
a_nominal  =  0.000054100081217;
b_nominal  = -0.011872008766780;
P0_nominal = 126.305; % Population in millions (2020)
t_start = 2020;
t_end   = 2100;
h       = 0.125; % Time step
t_vec   = t_start:h:t_end;
num_steps = length(t_vec);
%% 2. Monte Carlo Simulation Setup
num_simulations = 10000;
% Standard deviations based on percentage uncertainties
sigma_a  = 0.05 * abs(a_nominal);  % 5% uncertainty
sigma_b  = 0.05 * abs(b_nominal);  % 5% uncertainty
sigma_P0 = 0.01 * P0_nominal;      % 1% uncertainty
% Generate Gaussian parameter samples
rng(42); % Set seed for reproducibility
a_samples  = a_nominal + sigma_a * randn(num_simulations, 1);
b_samples  = b_nominal + sigma_b * randn(num_simulations, 1);
P0_samples = P0_nominal + sigma_P0 * randn(num_simulations, 1);
%% 3. Euler's Method Simulation Loop (Vectorized)
P_matrix = zeros(num_simulations, num_steps);
P_matrix(:, 1) = P0_samples;
for i = 1:(num_steps - 1)
    % Derivative: dP/dt = a*P^2 + b*P
    dPdt = a_samples .* (P_matrix(:, i).^2) + b_samples .* P_matrix(:, i);
    % Euler update step: P(t + h) = P(t) + h * dP/dt
    P_matrix(:, i+1) = P_matrix(:, i) + h * dPdt;
end
%% 4. Statistical Calculations (Percentiles)
p_median  = median(P_matrix, 1);
p_80_low  = prctile(P_matrix, 10, 1);
p_80_high = prctile(P_matrix, 90, 1);
p_95_low  = prctile(P_matrix, 2.5, 1);
p_95_high = prctile(P_matrix, 97.5, 1);
%% 5. United Nations Projection Data & Historical Data Used
un_years = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100];
un_pop   = [126.305, 119.584, 112.158, 105.123, 98.103, 90.611, 84.729, 80.634, 76.845];
% Historical population data used
data_years = [2020, 2023, 2026];
data_pop   = [126.305, 124.371, 122.428];
%% 6. Visualization
fig = figure('Color', 'w', 'Position', [100, 100, 900, 550]);
hold on;
% Shaded region for 95% Confidence Interval (Light Blue)
fill_x  = [t_vec, fliplr(t_vec)];
fill_95 = [p_95_high, fliplr(p_95_low)];
fill(fill_x, fill_95, [0.75, 0.88, 0.98], 'EdgeColor', 'none', ...
    'DisplayName', '95% CI');
% Shaded region for 80% Confidence Interval (Blue)
fill_80 = [p_80_high, fliplr(p_80_low)];
fill(fill_x, fill_80, [0.40, 0.65, 0.90], 'EdgeColor', 'none', ...
    'DisplayName', '80% CI');
% Median trajectory (Red)
plot(t_vec, p_median, 'r-', 'LineWidth', 3.0, 'DisplayName', 'Calculated (Median)');
% UN Population Projections (Black square markers)
plot(un_years,un_pop, 'o-', ...
     'Color',[0.10 0.10 0.10], ...
     'LineWidth',2.0, ...
     'MarkerFaceColor','w', 'DisplayName', 'UN Projection','MarkerSize',8);
% Data Used (Black dots)
plot(data_years, data_pop, 'ko', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'k', 'LineStyle', 'none', 'DisplayName', 'Data used');
% Chart Formatting
xlabel('Year', 'FontSize', 18, 'FontWeight', 'normal');
ylabel('Population (Millions)', 'FontSize', 18, 'FontWeight', 'normal');
title('Example 3- Japan -Linear Splines - Monte Carlo UQ', ...
    'FontSize', 22, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 15);
set(gca, 'FontSize', 18)
grid on;
box on;
xlim([2020, 2100]);
hold off;

%% 7. Export in PDF 
% Configuration for PDF 
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperPositionMode', 'manual');

% Fitting to page
paperWidth = 11;
paperHeight = 11 * (550 / 900);
set(fig, 'PaperSize', [paperWidth, paperHeight]);
set(fig, 'PaperPosition', [0, 0, paperWidth, paperHeight]);

% Exportation to PDF
pdfFilename = 'C:\Users\Pc03\Desktop\Example3_Japan_MC.pdf';

if exist('exportgraphics', 'file') == 2
    % Para MATLAB R2020a y posteriores (ajuste perfecto sin márgenes)
    exportgraphics(fig, pdfFilename, 'ContentType', 'vector');
else
    % for other MATLAB version or GNU Octave
    print(fig, pdfFilename, '-dpdf', '-r300', '-fillpage');
end