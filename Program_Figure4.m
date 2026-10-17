% ================================================================
% LINEAR SPLINES MONTE CARLO UQ WORLD 2024 COMPARED TO LOGISTIC (2024–2100)
% 1. Quadratic Model with Monte Carlo Uncertainty Quantification
% 2. Logistic Model Fitted via Optimization (fminsearch + ode45)
% ================================================================

clear;
clc;
close all;

% ================================================================
% DATA-SET FORM 2024 UN
x = [2017   2019   2021  2023];

P = [7.645  7.811  7.954 8.091];
%length(P)
% ================================================================
% COMPUTING k VECTOR
for i = 1: length(P)-1
k(i) =(1/2)* log(P(i+1)/P(i));
end
k
%% ---------------------------------------------------------------
% 0. SOLVING MATRIX FOR LINEAR SPLINES
k=[ 0.010740578623567   0.009070975078669   0.008538693020102];
%P=[5810000000000000   6144000000000000   6471000000000000   6801000000000000   7141000000000000   7492000000000000];
x=[7.811  7.954 8.091];
%plot(x,P)
%grid

M = [x(1)^1 x(1)^0 0 0 ; %row 1
%-------------------------------
     x(2)^1 x(2)^0 0 0 ; % row2
%-------------------------------
     0 0 x(2)^1 x(2)^0 ; %row 3
%-------------------------------
     0 0 x(3)^1 x(3)^0 ]; % row 4

      
 b= [k(1); k(2); k(2); k(3)];
 
 y = inv(M)*b %Formula
%% ---------------------------------------------------------------
% 1. Observed Data and Common Setup
% ---------------------------------------------------------------
% Observed data used for calibration
data_used_years = [2019 2021 2023];
data_used_population = [7.811 7.888 8.091];

% UN Population projections (UN 2024 report)
UN_years = [2019 2021 2023 2030 2040 2050 2060 2070 2080 2100];
UN_population = [7.811 7.888 8.091 8.57 9.20 9.66 10.27 10.40 10.42 10.19];

% Time domain settings
t0 = 2019;
tf = 2100;
P0 = data_used_population(1);

%% ---------------------------------------------------------------
% 2. Model 1: Quadratic Differential Equation + Monte Carlo
%      dP/dt = a*P^2 + b*P
% ---------------------------------------------------------------
a0 = -0.003885270500489;
b0 =  0.039974416639559;
K=-b0/a0
% Numerical parameters
h = 0.125;                     % Time step (years)
t_quad = t0:h:tf;
nt_quad = length(t_quad);
N = 1000;                      % Number of Monte Carlo simulations

% Uncertainty standard deviations
sigma_a  = 0.03;               % 3% relative uncertainty
sigma_b  = 0.03;               % 3% relative uncertainty
sigma_P0 = 0.01;               % 1% relative uncertainty

% Deterministic Euler Solution
P_quad_det = zeros(1, nt_quad);
P_quad_det(1) = P0;
for k = 1:nt_quad - 1
    dPdt = a0 * P_quad_det(k)^2 + b0 * P_quad_det(k);
    P_quad_det(k+1) = P_quad_det(k) + h * dPdt;
end

% Monte Carlo Simulation
P_MC = zeros(N, nt_quad);
randn('seed', 12345);            % Reproducibility seed

for j = 1:N
    a_rand = a0 * (1 + sigma_a * randn());
    b_rand = b0 * (1 + sigma_b * randn());
    P_init_rand = P0 * (1 + sigma_P0 * randn());
    
    P_MC(j, 1) = P_init_rand;
    for k = 1:nt_quad - 1
        dPdt = a_rand * P_MC(j, k)^2 + b_rand * P_MC(j, k);
        P_MC(j, k+1) = P_MC(j, k) + h * dPdt;
    end
end

% Monte Carlo Percentiles and Confidence Intervals
CI95 = prctile(P_MC, [2.5 97.5], 1);
CI80 = prctile(P_MC, [10 90], 1);
P_quad_median = prctile(P_MC, 50, 1);

P95_low  = CI95(1, :);
P95_high = CI95(2, :);
P80_low  = CI80(1, :);
P80_high = CI80(2, :);

%% ---------------------------------------------------------------
% 3. Model 2: Logistic Differential Equation
%      dP/dt = r*P*(1 - P/K)
% ---------------------------------------------------------------
% Objective function for parameter estimation (Least Squares)
logistic_func = @(theta, t) theta(2) ./ (1 + ((theta(2) - P0) / P0) .* exp(-theta(1) .* (t - t0)));
objective = @(theta) sum((logistic_func(theta, data_used_years) - data_used_population).^2);

% Initial parameter guesses: theta0 = [r0, K0]
theta0 = [0.02, 12];
options = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-10);

% Optimization via fminsearch
theta_opt = fminsearch(objective, theta0, options);
r_opt = theta_opt(1);
K_opt = theta_opt(2);

% Numerical Solution via ode45
tspan_log = [t0 tf];
logistic_ode = @(t, P) r_opt * P * (1 - P / K_opt);
[t_log, P_log] = ode45(logistic_ode, tspan_log, P0);

% Performance metric (RMSE on observed points)
P_log_calc = logistic_func(theta_opt, data_used_years);
RMSE_log = sqrt(mean((P_log_calc - data_used_population).^2));

%% ---------------------------------------------------------------
% 4. Visualization
% ---------------------------------------------------------------
fig = figure('Color', 'w');
hold on;
box on;
grid on;

% 95% Confidence Interval (Quadratic Model)
fill([t_quad fliplr(t_quad)], [P95_low fliplr(P95_high)], ...
     [0.75 0.85 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.55);

% 80% Confidence Interval (Quadratic Model)
fill([t_quad fliplr(t_quad)], [P80_low fliplr(P80_high)], ...
     [0.35 0.60 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.45);

% Quadratic Model Median Curve
plot(t_quad, P_quad_median, 'r-', 'LineWidth', 2.5);

% Logistic Model Curve (ode45)
plot(t_log, P_log, 'b--', 'LineWidth', 1.8);

% UN Population Projection Benchmark
%plot(UN_years, UN_population, 'k--', 'LineWidth', 1.5);

plot(UN_years,UN_population, ...
    'o-', ...
     'Color',[0.10 0.10 0.10], ...
     'LineWidth',1.2, ...
     'MarkerFaceColor','w', 'DisplayName', 'UN Projection','MarkerSize',5);


% Observed Calibration Data Points
plot(data_used_years, data_used_population, 'ko', ...
     'MarkerFaceColor', 'k', 'MarkerSize', 7);

% Labels and Formatting
xlabel('Year', 'FontSize', 13);
ylabel('Population (billions)', 'FontSize', 13);
title('Example 1 - World 2024 UN - Linear Splines - Monte Carlo UQ', 'FontSize', 14);
legend('95% CI', ...
       '80% CI', ...
       'Calculated (Median)', ...
       'Logistic Model', ...
       '2024 UN Projection', ...
       'Data used', ...
       'Location', 'southeast');
xlim([2019 2100]);
ylim([7 12]);
set(gca, 'FontSize', 10);

% Export PDF figure
set(fig, 'Units', 'Centimeters');
pos = get(fig, 'Position');
set(fig, 'PaperUnits', 'Centimeters');
set(fig, 'PaperSize', [pos(3) pos(4)]);
set(fig, 'PaperPositionMode', 'manual');
set(fig, 'PaperPosition', [0 0 pos(3) pos(4)]);
exportgraphics(fig, 'C:\Users\Pc03\Desktop\Example1_WP2024_MC.pdf', 'ContentType', 'vector');

%% ---------------------------------------------------------------
% 5. Print Output Summary
% ---------------------------------------------------------------
fprintf('\n====================================================\n');
fprintf('        COMBINED WORLD POPULATION MODEL SUMMARY      \n');
fprintf('====================================================\n');
fprintf('Time Period          : %d to %d\n', t0, tf);
fprintf('Observed Data Points : %s\n', mat2str(data_used_years));
fprintf('----------------------------------------------------\n');
fprintf('QUADRATIC MODEL (Euler + Monte Carlo):\n');
fprintf('  a  = %.15f\n', a0);
fprintf('  b  = %.15f\n', b0);
fprintf('  2100 Population (Median) : %.4f billion\n', P_quad_median(end));
fprintf('  2100 Population (80%% CI) : [%.4f, %.4f] billion\n', P80_low(end), P80_high(end));
fprintf('  2100 Population (95%% CI) : [%.4f, %.4f] billion\n', P95_low(end), P95_high(end));
fprintf('----------------------------------------------------\n');
fprintf('LOGISTIC MODEL (fminsearch + ode45):\n');
fprintf('  r  = %.12f (growth rate)\n', r_opt);
fprintf('  K  = %.12f billion (carrying capacity)\n', K_opt);
fprintf('  RMSE (Calibration)        : %.12f billion\n', RMSE_log);
fprintf('  2100 Population Projection : %.4f billion\n', P_log(end));
fprintf('====================================================\n');