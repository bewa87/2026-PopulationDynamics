%% ================================================================
%  AFRICA MONTE CARLO UQ LINEAR SPLINES
%  dP/dt = a*P^2 + b*P
%
%  Euler Method
%  Uncertainty Quantification via Monte Carlo
%  1000 particles
%
%  Population in billions of inhabitants
% ================================================================
clear;
clc;
close all;

%  FROM DATA-SET CVS FILE WE COMPUTE k VECTOR
% ================================================================
x = [2010   2012   2014  2016]; %year

P = [1.065   1.123   1.186   1.251]; %ppulation in billions


for i = 1: length(P)-1
k(i) =(1/2)* log(P(i+1)/P(i));
end
k

%  SOLVING THE MATRIX FOR LINEAR SPLINES
% ================================================================
k=[0.026514438297459   0.027291312409614   0.026678465454620];

x=[1.123   1.186   1.251];


M = [x(1)^1 x(1)^0 0 0 ; %row 1 
%-------------------------------
     x(2)^1 x(2)^0 0 0 ; %row 2 
%-------------------------------
     0 0 x(2)^1 x(2)^0 ; %row 3
%-------------------------------
     0 0 x(3)^1 x(3)^0 ]; %row 4

      
 b= [k(1); k(2); k(2); k(3)];
 
 y = inv(M)*b %COEFFICIENTS
 
%% ------------------------------------------------
% 1. MODEL PARAMETERS
% -------------------------------------------------
a0 = -0.009428414692215;
b0 =  0.038473412234581;
P0 =  1.123;                     % billions
h = 0.125;                       % Euler step size
Nmc = 1000;                      % number of particles
anio_inicio = 2010;
anio_final  = 2100;
%% ------------------------------------------------
% 2. PARAMETER UNCERTAINTY
% -------------------------------------------------
% A normal distribution with CV = 5% is assumed.
CV_a  = 0.05;
CV_b  = 0.05;
CV_P0 = 0.05;
sigma_a  = abs(a0)*CV_a;
sigma_b  = abs(b0)*CV_b;
sigma_P0 = P0*CV_P0;
% Reset to default generator to avoid legacy error
rng('default');
rng(1);       % Seed for reproducibility
% Generation of 1000 particles
a_MC  = a0  + sigma_a*randn(Nmc,1);
b_MC  = b0  + sigma_b*randn(Nmc,1);
P0_MC = P0  + sigma_P0*randn(Nmc,1);
%% ------------------------------------------------
% 3. TIME MESH
% -------------------------------------------------
t = anio_inicio:h:anio_final;
Nt = length(t);
P_MC = zeros(Nmc,Nt);
%% ------------------------------------------------
% 4. INITIAL CONDITION
% -------------------------------------------------
indice_2012 = round((2012-anio_inicio)/h) + 1;
indice_2010 = 1;
P_MC(:,indice_2012) = P0_MC;
%% ------------------------------------------------
% 5. BACKWARD EULER: 2012 -> 2010
% -------------------------------------------------
for j = indice_2012:-1:(indice_2010+1)
    P_actual = P_MC(:,j);
    P_MC(:,j-1) = P_actual - h .* ...
        (a_MC.*P_actual.^2 + b_MC.*P_actual);
end
%% ------------------------------------------------
% 6. FORWARD EULER: 2012 -> 2100
% -------------------------------------------------
for j = indice_2012:Nt-1
    P_actual = P_MC(:,j);
    dP = a_MC.*P_actual.^2 + b_MC.*P_actual;
    P_MC(:,j+1) = P_actual + h*dP;
end
%% ------------------------------------------------
% 7. MONTE CARLO STATISTICS
% -------------------------------------------------
P_media = mean(P_MC,1);
P_75_inf = prctile(P_MC,12.5,1);
P_75_sup = prctile(P_MC,87.5,1);
P_95_inf = prctile(P_MC,2.5,1);
P_95_sup = prctile(P_MC,97.5,1);
%% ------------------------------------------------
% 8. UNITED NATIONS DATA
% -------------------------------------------------
anios_UN = [ 2020 2030 2040 2050 ...
            2060 2070 2080 2090 2100];
%deleted (2012,1.123)
Africa = [ 1.360 1.727 2.096 2.467 ...
          2.823 3.160 3.450 3.690 3.814];
%% ------------------------------------------------
% 9. OBSERVED DATA POINTS
% -------------------------------------------------
anios_obs = [2012 2014 2016];
P_obs = [1.123 1.184 1.247];
%% ------------------------------------------------
% 10. PREDICTIONS 2017-2025
% -------------------------------------------------
%anios_pred = 2017:2025;
%P_pred = interp1(t,P_media,anios_pred,'linear');
%% ------------------------------------------------
% 11. PLOT (MEJORAS VISUALES APLICADAS)
% -------------------------------------------------
fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);
hold on;
grid on;
box on;
%% 95% Confidence Interval (BLUE)
fill([t fliplr(t)], ...
     [P_95_inf fliplr(P_95_sup)], ...
     [0.65 0.80 0.98], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.65);
%% 75% Confidence Interval (DEEP BLUE )
fill([t fliplr(t)], ...
     [P_75_inf fliplr(P_75_sup)], ...
     [0.15 0.45 0.90], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.65);
%% Monte Carlo Mean (Red Line)
plot(t,P_media, ...
     'r-', ...
     'LineWidth',4.5);
%% Observed data: black dots
plot(anios_obs,P_obs, ...
     'ko', ...
     'MarkerFaceColor','k', ...
     'MarkerSize',11, ...
     'LineWidth',1.8);
%% United Nations Curve
plot(anios_UN,Africa,...
    'o-', ...
     'Color',[0.10 0.10 0.10], ...
     'LineWidth',2.2, ...
     'MarkerFaceColor','w', ...
     'DisplayName', 'UN Projection', ...
     'MarkerSize',9);
%% ------------------------------------------------
% 12. PLOT FORMATTING
% -------------------------------------------------
xlabel('Year','FontSize',30,'FontWeight','normal');
ylabel('Population (billions)','FontSize',30,'FontWeight','normal');
title(['Example 2 - Africa - Linear Spline - Monte Carlo UQ'], ...
       'FontSize',34,'FontWeight','bold');
xlim([2010 2100]);
legend({'95% CI','75% CI', ...
        'Calculated (Median)', ...
        'Data used', ...
        '2024 UN Projection'}, ...
        'Location','northwest', ...
        'FontSize',20);
set(gca,'FontSize',24,'LineWidth',1.5);
hold off;
%% ------------------------------------------------
% 13. FULL PAGE PDF SETUP AND EXPORT 
% -------------------------------------------------
% Page Configuration  (21.0 cm -> 25.2 cm)
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperOrientation', 'landscape');
set(fig, 'PaperSize', [29.7 25.2]);
set(fig, 'PaperPosition', [0 0 29.7 25.2]);
% Full exportation
exportgraphics(fig, 'C:\Users\Pc03\Desktop\Example2_Africa_MC.pdf', 'ContentType', 'vector');
%% ------------------------------------------------
% 14. DISPLAY NUMERICAL RESULTS
% -------------------------------------------------
fprintf('\n============================================\n');
fprintf('     MONTE CARLO RESULTS\n');
fprintf('============================================\n');
fprintf('Number of particles = %d\n',Nmc);
fprintf('Euler step size h   = %.3f years\n',h);
fprintf('\nNominal parameters:\n');
fprintf('a  = %.15f\n',a0);
fprintf('b  = %.15f\n',b0);
fprintf('P0 = %.3f billions\n',P0);
fprintf('\nUncertainties used:\n');
fprintf('CV(a)  = %.1f %%\n',100*CV_a);
fprintf('CV(b)  = %.1f %%\n',100*CV_b);
fprintf('CV(P0) = %.1f %%\n',100*CV_P0);
%% ------------------------------------------------
% 15. RESULTS TABLE FOR SELECTED YEARS
% -------------------------------------------------
anios_tabla = [2012 2014 2016 2017 2020 2025 ...
               2030 2040 2050 2060 2070 2080 ...
               2090 2100];
media_tabla = interp1(t,P_media,anios_tabla);
inf75_tabla = interp1(t,P_75_inf,anios_tabla);
sup75_tabla = interp1(t,P_75_sup,anios_tabla);
inf95_tabla = interp1(t,P_95_inf,anios_tabla);
sup95_tabla = interp1(t,P_95_sup,anios_tabla);
Tabla = table(anios_tabla', ...
              media_tabla', ...
              inf75_tabla', ...
              sup75_tabla', ...
              inf95_tabla', ...
              sup95_tabla', ...
              'VariableNames', ...
              {'Year','Mean','CI75_Inf','CI75_Sup',...
               'CI95_Inf','CI95_Sup'});
disp(Tabla);