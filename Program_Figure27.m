% =========================================================================
% UNCERTAINTY QUANTIFICATION (MONTE CARLO) FOR OUR QUADRATIC TUMOR GROWTH MODEL
% COMPARISON WITH LOGISTIC AND GOMPERTZ MODELS
% =========================================================================
clear; clc; close all;

%% 0. Experimental Data in cvs and calculation of k vector
x = [0 7 12 14];

P = [1   13.57 47.27 68.75];

k(1) =(1/7)* log(P(2)/P(1));
k(2) =(1/5)* log(P(3)/P(2));
k(3) =(1/2)* log(P(4)/P(3));
k

%% 0. Solving Matrix for two linear Splines
k=[0.372551639120968   0.249602874191673   0.187300445870768];

x=[13.57 47.27 68.75];

M = [x(1)^1 x(1)^0 0 0 ; % row 1 
%-------------------------------
     x(2)^1 x(2)^0 0 0 ; %row 2 
%-------------------------------
     0 0 x(2)^1 x(2)^0 ; % row 3
%-------------------------------
     0 0 x(3)^1 x(3)^0 ]; % row 4

      
 b= [k(1); k(2); k(2); k(3)];
 
 y = inv(M)*b %COEFFICIENTS
%% 1. Experimental Data
% Initial training data (Days 0 to 14)
t_data = [0, 7, 12, 14];           % Time vector (days)
V_data = [1, 13.57, 47.27, 68.75]; % Tumor volume vector
V0 = V_data(1);                    % Initial volume V(0) = 1
t0 = 0;                            % Initial time
tf = 20;                           % Final simulation time (days)

% Unused experimental data points
t_unused_data = [11, 13];
V_unused_data = [46.30, 58.72];

% Final experimental test data points
t_final_data = [15, 18, 19];
V_final_data = [73.33, 109.77, 122.15];

%% 2. Model 1: Quadratic Model with Monte Carlo Uncertainty Quantification
% Differential Equation: dP/dt = a*P^2 + b*P
a_quad = -0.002900485489800;
b_quad =  0.386708823294521;
h = 0.125;                         % Time step size for Euler's method
t_euler = t0:h:tf;
N_steps = length(t_euler);
K = -b_quad/a_quad
% Monte Carlo Parameters
N_mc = 1000;                       % Number of Monte Carlo iterations
std_ratio = 0.05;                  % 5% perturbation standard deviation
rng(42); % Set random seed for reproducibility

% Perturb parameters a, b, and P(0) using Gaussian distribution
a_mc  = a_quad + std_ratio * abs(a_quad) * randn(N_mc, 1);
b_mc  = b_quad + std_ratio * abs(b_quad) * randn(N_mc, 1);
V0_mc = V0     + std_ratio * V0          * randn(N_mc, 1);

% Matrix to store simulation results: rows = MC samples, cols = time steps
V_mc = zeros(N_mc, N_steps);

for j = 1:N_mc
    a_j = a_mc(j);
    b_j = b_mc(j);
    V_mc(j, 1) = V0_mc(j);
    
    for i = 1:(N_steps - 1)
        f_val = a_j * V_mc(j, i)^2 + b_j * V_mc(j, i);
        V_mc(j, i+1) = V_mc(j, i) + h * f_val;
    end
end

% Compute Quantiles across simulations for Confidence Intervals
median_quad = prctile(V_mc, 50, 1);     % Median (50th percentile)
ci80_low    = prctile(V_mc, 10, 1);     % 80% CI lower bound (10th percentile)
ci80_high   = prctile(V_mc, 90, 1);     % 80% CI upper bound (90th percentile)
ci95_low    = prctile(V_mc, 2.5, 1);    % 95% CI lower bound (2.5th percentile)
ci95_high   = prctile(V_mc, 97.5, 1);   % 95% CI upper bound (97.5th percentile)

%% 3. Model 2: Logistic Growth Model
logistic_model = @(p, t) p(1) ./ (1 + ((p(1) - V0) / V0) * exp(-p(2) * t));
fun_error_log  = @(p) sum((V_data - logistic_model(p, t_data)).^2);
p_log_init     = [100, 0.3];
p_log_opt      = fminsearch(fun_error_log, p_log_init);

K_est = p_log_opt(1);
r_est = p_log_opt(2);

ode_logistic = @(t, V) r_est * V * (1 - V / K_est);
[t_log, V_log] = ode45(ode_logistic, [t0 tf], V0);

%% 4. Model 3: Gompertz Growth Model
gompertz_model = @(p, t) V0 * exp( (p(1)/p(2)) * (1 - exp(-p(2)*t)) );
p_gomp_init    = [0.5, 0.1];

if exist('lsqcurvefit', 'file') == 2
    options = optimoptions('lsqcurvefit', 'Display', 'off');
    p_gomp_opt = lsqcurvefit(gompertz_model, p_gomp_init, t_data, V_data, [0, 0], [], options);
else
    gomp_error = @(p) sum((V_data - gompertz_model(p, t_data)).^2);
    p_gomp_opt = fminsearch(gomp_error, p_gomp_init);
end

t_gomp = linspace(t0, tf, 200);
V_gomp = gompertz_model(p_gomp_opt, t_gomp);

%% 5. Unified Visualization
fig = figure('Name', 'Example 5 - Total Volume Fraction -Linear Spline Monte Carlo UQ', ...
             'Color', 'w', 'Position', [100, 100, 1600, 1300]);
hold on;

% --- 95% Confidence Interval (Light Blue Shading) ---
fill([t_euler, fliplr(t_euler)], [ci95_high, fliplr(ci95_low)], ...
     [0.70, 0.88, 1.00], 'EdgeColor', 'none', ...
     'DisplayName', '95% CI');

% --- 80% Confidence Interval (Blue Shading) ---
fill([t_euler, fliplr(t_euler)], [ci80_high, fliplr(ci80_low)], ...
     [0.20, 0.45, 0.85], 'EdgeColor', 'none', ...
     'DisplayName', '80% CI');

% --- Median Curve (Red Solid Line) ---
plot(t_euler, median_quad, 'r-', 'LineWidth', 2.7, ...
     'DisplayName', 'Calculated (Median)');

% --- Logistic Model (Black Dashed Line) ---
plot(t_log, V_log, 'k--', 'LineWidth', 2.3, ...
     'DisplayName', 'Logistic Model');

% --- Gompertz Model (Green Dashed Line) ---
plot(t_gomp, V_gomp, 'g--', 'LineWidth', 2.3, ...
     'DisplayName', 'Gompertz Model');

% --- Initial Experimental Data (Black Circles) ---
plot(t_data, V_data, 'ko', 'MarkerSize', 11, 'MarkerFaceColor', 'k', ...
     'DisplayName', 'Data used');

% --- Unused Experimental Data Points (Gray Circles) ---
plot(t_unused_data, V_unused_data, 'o', 'MarkerSize', 11, ...
     'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], ...
     'DisplayName', 'Unused Data');

% --- Final Experimental Data Points (Red Circles) ---
plot(t_final_data, V_final_data, 'ro', 'MarkerSize', 12, ...
     'MarkerFaceColor', 'r', 'DisplayName', 'Final Data Points (15, 18, 19)');

% --- Vertical Line at Day 14 and Text Label ---
xline(14, 'k:', 'LineWidth', 2.0, 'HandleVisibility', 'off');
text(14.2, 20, 'end of data used', 'FontSize', 16, 'FontAngle', 'italic', 'Color', 'k');

% --- Formatting ---
grid on;
xlabel('Time (days)', 'FontSize', 22, 'FontWeight', 'normal');
ylabel('Tumor Volume V(t)', 'FontSize', 22, 'FontWeight', 'normal');
title('Example 5 - Total Volume Fraction -Linear Spline Monte Carlo UQ', ...
      'FontSize', 19, 'FontWeight', 'bold');
set(gca, 'FontSize', 20)
legend('Location', 'northwest', 'FontSize', 20);
xlim([t0, tf]);
ylim([0, max([ci95_high, V_log', V_gomp, V_final_data]) * 1.05]);
hold off;

%% 6. Export Figure to PDF (Full Page)
% Configuración de tamaño de papel sin bordes/márgenes excedentes
set(fig, 'Units', 'Inches');
pos = get(fig, 'Position');
set(fig, 'PaperUnits', 'Inches');
set(fig, 'PaperSize', [pos(3), pos(4)]);
set(fig, 'PaperPositionMode', 'manual');
set(fig, 'PaperPosition', [0, 0, pos(3), pos(4)]);

% Exportación vectorial limpia (MATLAB R2020a o superior)
exportgraphics(fig, 'C:\Users\Pc03\Desktop\Example5_Tumor.pdf', 'ContentType', 'vector');

