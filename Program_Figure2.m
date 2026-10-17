% WORLD UNITED NATIONS 2019 Splines 
% constant k vs Population P
% Two linear splines t=2 de 2014 al 2018
format long  
clear all, clc

% From files CVS for World Population UN 2019
year = [2012 2014 2016 2018];
Po = [7.125 7.295 7.464 7.631];

% Calculate the constants 
for i = 1: length(Po)-1
    k(i) = (1/2) * log(Po(i+1)/Po(i));
end
k
k = [0.011789727910017   0.011451140785147   0.011063717381227]; % Constants vector
P = [7.295 7.464 7.631]; % Population vector

% Matrix
M = [P(1)^1 P(1)^0 0 0 ; % row 1
     P(2)^1 P(2)^0 0 0 ; % row 2
     0 0 P(2)^1 P(2)^0 ; % row 3
     0 0 P(3)^1 P(3)^0 ]; % row 4

b = [k(1); k(2); k(2); k(3)]; 
y = inv(M)*b; % Formula

x1 = [7.295:0.001:7.464]; % Evaluating spline 1
S1 = y(1)*x1 + y(2);

x2 = [7.464:0.001:7.631]; % Evaluating spline 2
S2 = y(3)*x2 + y(4);

y

% ------------------- CODE FOR FIGURE AND EXPORTATION TO PDF -------------------
% Create figure 
fig = figure('Visible', 'off');

% Plot S1 (red, line) y S2 (black, line)
%plot(P, k, 'bo', 'MarkerFaceColor', 'b', 'LineWidth', 2);% hold on;
plot(x1, S1, 'r-', 'LineWidth', 2); hold on;
plot(x2, S2, 'k-', 'LineWidth', 2);
plot(P, k, 'bo', 'MarkerFaceColor', 'b', 'LineWidth', 2);% hold on;
% Title
title('Example 1 - Dependence of k on N -Linear spline', 'FontSize', 14);

% Axis labels
xlabel('World Population (billions)', 'FontSize', 12);
ylabel('k(N)', 'FontSize', 12);

% Inset / legend
legend({'Spline_1', 'Spline_2','Data used'}, 'FontSize', 12, 'Location', 'best');

% GCA font
set(gca, 'FontSize', 12);
grid on;
hold off;

% --- Configuration for page in PDF ---
% 1. Define size (Landscape): 11 x 8.5 pulgadas
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [11, 8.5]);

% 2. Fit size
set(fig, 'PaperPosition', [0, 0, 11, 8.5]);
set(fig, 'PaperPositionMode', 'manual');

% 3. Save / Export figure in PDF
exportgraphics(fig, 'C:\Users\Pc03\Desktop\Example1_WP_Splines.pdf', 'ContentType', 'vector');

% Figure created
close(fig);

disp('PDF has been created.');