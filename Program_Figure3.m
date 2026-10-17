% ================================================================
% WORLD POPULATION 2019 UN, 2 LINEAR SPLINES
% 
clear all, clc
% ================================================================
% CONSTANTS
%Data from CVS files
x = [2012 2014 2016 2018]; %year

P = [7.125 7.295 7.464 7.631]; %population
for i = 1: length(P)-1
k(i) =(1/2)* log(P(i+1)/P(i));
end
k %displays the k vector
% ================================================================
% SOLVING MATRIX FOR TWO SPLINES
%Matix

k=[  0.011789727910017   0.011451140785147   0.011063717381227]; %Constants vector

x=[7.295 7.464 7.631]; %Population vector
M = [x(1)^1 x(1)^0 0 0 ; %  row 1
%-------------------------------
     x(2)^1 x(2)^0 0 0 ; %2 row 2
%-------------------------------
     0 0 x(2)^1 x(2)^0 ; %  row 3
%-------------------------------
     0 0 x(3)^1 x(3)^0 ]; % row 4

      
 b= [k(1); k(2); k(2); k(3)]; 
 
 y = inv(M)*b %Formula
 
%-------------------------------
 %x=[7.295 7.464 7.631];
 x1 = [7.295:0.001:7.464]; %Evaluating spline 1
 S1 =  y(1)*x1+ y(2);
 x2 = [7.464:0.001:7.631]; %Evaluating spline 2
 S2 =  y(3)*x2+ y(4);
% ================================================================
% SOLVING EQUATION 
%
%       dP/dt = a*P^2 + b*P
%
% Numerical Method Explicit Euler
% Step: h = 0.125 anos
% Monte Carlo: N = 1000
%
% Population in billions
%
% United Nations Data:
% 2014  7.295
% 2016  7.464
% 2018  7.631
% 2020  7.795
% 2030  8.549
% 2040  9.199
% 2050  9.735
% 2060 10.175
% 2070 10.595
% 2080 10.956
% 2090 11.241
% 2100 11.465
% ================================================================

%% ---------------------------------------------------------------
% 1. MODEL PARAMETERS
% ---------------------------------------------------------------
% From Program called "Program_Figure1_Splines" we obtain the model parameters
a0 = -0.002319900622275; 
b0 =  0.028766879029811;
P0 =  7.295;
K=-b0/a0
%  Euler step
h = 0.125;
% Time interval
t0 = 2014;
tf = 2100;
% Number Monte Carlo simulations
N = 1000;
%% ---------------------------------------------------------------
% 2. UNITED NATIONS DATA
% ---------------------------------------------------------------
t_UN = [2014 2016 2018 2020 2030 2040 2050 2060 ...
        2070 2080 2090 2100];
P_UN = [7.295 7.464 7.631 7.795 8.549 9.199 9.735 10.175 ...
        10.595 10.956 11.241 11.465];
% Data used for calibration
t_obs = [2014 2016 2018];
P_obs = [7.295 7.464 7.631];
%% ---------------------------------------------------------------
% 3. TIME MESH
% ---------------------------------------------------------------
t = t0:h:tf;
nt = length(t);
%% ---------------------------------------------------------------
% 4. NUMERICAL SOLUTION THROUGH EULER
% ---------------------------------------------------------------
P_det = zeros(1,nt);
P_det(1) = P0;
for k = 1:nt-1
    dPdt = a0*P_det(k)^2 + b0*P_det(k);
    P_det(k+1) = P_det(k) + h*dPdt;
end
%% ---------------------------------------------------------------
% 5. MONTE CARLO
%
% We consider independent uncertainties for:
%
%       a
%       b
%       P0
%
% LasPerturbations are generated through normal distributions.
%
%
% ---------------------------------------------------------------
% Relative uncertainties
sigma_a_rel  = 0.03;     % 3 %
sigma_b_rel  = 0.03;     % 3 %
sigma_P0_rel = 0.01;     % 1 %
% Matrix for storing simulations
P_MC = zeros(N,nt);
% Seed for reproducibility
rng(12345); % Updated randn("seed", ...) to modern MATLAB rng syntax
for j = 1:N
    % -----------------------------------------------------------
    % Random parameter generation
    % -----------------------------------------------------------
    a_MC = a0 * (1 + sigma_a_rel*randn());
    b_MC = b0 * (1 + sigma_b_rel*randn());
    P0_MC = P0 * (1 + sigma_P0_rel*randn());
    % Avoid non physical values for P0
    if P0_MC <= 0
        P0_MC = P0;
    end
    % -----------------------------------------------------------
    % Euler for Monte Carlo
    % -----------------------------------------------------------
    P_MC(j,1) = P0_MC;
    for k = 1:nt-1
        dPdt_MC = a_MC*P_MC(j,k)^2 + ...
                  b_MC*P_MC(j,k);
        P_MC(j,k+1) = P_MC(j,k) + h*dPdt_MC;
        % Protection against non physical values
        if P_MC(j,k+1) < 0
            P_MC(j,k+1) = NaN;
        end
    end
end
%% ---------------------------------------------------------------
% 6. UNCERTAINTY QUANTIFICATION
%
% CI 80 %:
%       percentils 10 y 90
%
% CI 95 %:
%       percentills 2.5 y 97.5
%
% Median:
%       percentil 50
% ---------------------------------------------------------------
P_80 = prctile(P_MC,[10 90],1);
P_95 = prctile(P_MC,[2.5 97.5],1);
P_med = prctile(P_MC,50,1);
% Limits
P80_inf = P_80(1,:);
P80_sup = P_80(2,:);
P95_inf = P_95(1,:);
P95_sup = P_95(2,:);
%% ---------------------------------------------------------------
% 7. GRAPHICS
% ---------------------------------------------------------------
fig = figure('Color','w');
hold on;
grid on;
box on;
%% ---------------------------------------------------------------
% CI 95 %
% blue
% ---------------------------------------------------------------
fill([t fliplr(t)], ...
     [P95_sup fliplr(P95_inf)], ...
     [0.80 0.90 1.00], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.60);
%% ---------------------------------------------------------------
% CI 80 %
% light blue
% ---------------------------------------------------------------
fill([t fliplr(t)], ...
     [P80_sup fliplr(P80_inf)], ...
     [0.30 0.55 0.90], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.45);
%% ---------------------------------------------------------------
%  Monte Carlo Median
% Red
% ---------------------------------------------------------------
plot(t,P_med, ...
     'r-', ...
     'LineWidth',2.5);
%% ---------------------------------------------------------------
% Deterministic solution for our Spline
% Black discountinous line
% ---------------------------------------------------------------
%plot(t,P_det, ...
 %    'k--', ...
  %   'LineWidth',1.2);
%% ---------------------------------------------------------------
% United Nations Data
% ---------------------------------------------------------------
plot(t_UN,P_UN, ...
     'o-', ...
     'Color',[0.10 0.10 0.10], ...
     'LineWidth',1.2, ...
     'MarkerFaceColor','w', ...
     'MarkerSize',5);
%% ---------------------------------------------------------------
% Used data 2014, 2016 y 2018
% Black
% ---------------------------------------------------------------
plot(t_obs,P_obs, ...
     'ko', ...
     'MarkerFaceColor','k', ...
     'MarkerSize',7, ...
     'LineWidth',1.5);
%% ---------------------------------------------------------------
% Configuration of the graphics
% ---------------------------------------------------------------
xlabel('Year','FontSize',11);
ylabel('Population (billons)','FontSize',11);
title('Example1 - World 2019 UN - Linear Splines - Monte Carlo UQ', ...
      'FontSize',12);
xlim([2014 2100]);
ylim([6.5 12]);
xticks(2014:10:2100);
legend({'CI 95 %','CI 80 %', ...
        'Calculated (Median)', ...
        '2019 UN', ...
        'Data used'}, ...
        'Location','southeast');
set(gca,'FontSize',11);
hold off;

%% ---------------------------------------------------------------
% EXPORT TO PDF (FULL PAGE / ZERO MARGINS)
% ---------------------------------------------------------------
% Se ajustan las propiedades del papel de la figura para que coincidan 
% exactamente con la proporción de la pantalla y ocupe toda la página PDF.
set(fig, 'PaperPositionMode', 'auto');
fig_pos = get(fig, 'PaperPosition');
set(fig, 'PaperSize', [fig_pos(3) fig_pos(4)]);

% Exportación en PDF
exportgraphics(fig, 'C:\Users\Pc03\Desktop\Example1_WP2019_MC.pdf', 'ContentType', 'vector');

%% ---------------------------------------------------------------
% 8. Show parameters on screen
% ---------------------------------------------------------------
fprintf('\n');
fprintf('===============================================\n');
fprintf(' WORLD POPULATION MODEL\n');
fprintf('===============================================\n');
fprintf('a  = %.15f\n',a0);
fprintf('b  = %.15f\n',b0);
fprintf('P0 = %.3f billions\n',P0);
fprintf('h  = %.3f years\n',h);
fprintf('N  = %d Monte Carlo simulation\n',N);
fprintf('\n');
fprintf('Uncertainty used:\n');
fprintf('a  = %.1f %%\n',100*sigma_a_rel);
fprintf('b  = %.1f %%\n',100*sigma_b_rel);
fprintf('P0 = %.1f %%\n',100*sigma_P0_rel);
fprintf('===============================================\n');
%% ---------------------------------------------------------------
% 9. MEDIAN VALUES FOR SOME YEARS
% ---------------------------------------------------------------
anios_consulta = [2014 2016 2018 2020 2030 2040 ...
                  2050 2060 2070 2080 2090 2100];
fprintf('\n');
fprintf('Año       MC Median       UN\n');
fprintf('-----------------------------------\n');
for i = 1:length(anios_consulta)
    [~,idx] = min(abs(t-anios_consulta(i)));
    fprintf('%4d      %8.4f      %8.4f\n', ...
            anios_consulta(i), ...
            P_med(idx), ...
            P_UN(i));
end
fprintf('-----------------------------------\n');