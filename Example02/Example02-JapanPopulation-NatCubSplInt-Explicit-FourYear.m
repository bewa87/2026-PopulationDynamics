% Step 0: Clearing

clear all;
close all;
clc;

% Step 1: Load Japan data

pkg load statistics;

fid = fopen('Example02-Japan-Population.csv','r');
C   = textscan(fid,'%f32 %f32', 'delimiter', ';');
fclose(fid);

japan_year  = cell2mat(C(1));
japan_pop   = cell2mat(C(2));

start_ind   = find(japan_year==1992);
end_ind     = find(japan_year==2016);

% Step 2: Definition of Japan's population data points

time_int = 4;
year     = [japan_year(start_ind):time_int:japan_year(end_ind)]';
x        = [0:time_int:(japan_year(end_ind)-japan_year(start_ind))]';
y        = japan_pop(start_ind:time_int:end_ind);

% Step 3: Computation of k in dependence of y

k    = (1/time_int).*(log(y(2:1:end)./y(1:1:(end-1)))); %(1./x(2:1:end)).*log(y(2:1:end)/y(1));

figure(1)
plot(y(2:1:end),k,'marker','*')

% Step 4: Spline interpolation

y_spl = y(2:1:end);
k_spl = k;
n = (length(y_spl) - 1);
h = zeros(n,1);
h = diff(y_spl);

M = cubic_spline_natural(y_spl,k_spl,n,h);

cub_spl = cell(n,1);

for j = 1:1:n
  cub_spl{j} = @(z) (M(j)*(y_spl(j+1) - z).^3/(6*h(j))) + (M(j+1)*(z - y_spl(j)).^3/(6*h(j))) + (k_spl(j) - M(j)*h(j)^2/6)*((y_spl(j+1) - z)/h(j)) + (k_spl(j+1) - M(j+1)*h(j)^2/6)*((z - y_spl(j))/h(j));
endfor

m = 50;
figure(2)
for j = 1:1:n
  plot(y_spl(j):h(j)/m:y_spl(j+1),cub_spl{j}(y_spl(j):h(j)/m:y_spl(j+1)),'color','blue','linewidth',1)
  hold on
  plot(y_spl(j),k_spl(j),'marker','*','color','blue')
  hold on
  title('Example 2 - Dependence of k on P (Cubic spline)','fontsize',14);
  xlabel('Population P','fontsize',10);
  ylabel('k(P)','fontsize',10);
endfor
  plot(y_spl(end),k_spl(end),'marker','*','color','blue')

% Step 5: Explicit Eulerian Time-Stepping Method

h_eul    = 0.125;
x_eul    = [0:h_eul:x(end-1)+12]';
y_eul    = zeros(length(x_eul),1);
y_eul(1) = y(2);

m = 1;
for l = 1:1:(length(x_eul)-1)
    if (x_eul(l) < x(m+1))
      y_eul(l+1) = y_eul(l) + h_eul*(cub_spl{m}(y_eul(l)))*y_eul(l);
    elseif ((x_eul(l) >= x(m+1)) & (m < length(cub_spl)))%(m < length(cub_spl)))
      m          = m+1;
      y_eul(l+1) = y_eul(l) + h_eul*(cub_spl{m}(y_eul(l)))*y_eul(l);
    else
      y_eul(l+1) = y_eul(l) + h_eul*(cub_spl{length(cub_spl)}(y_eul(l)))*y_eul(l);
    endif
endfor

% Step 6: Rigid Logistic Growth Model

% Define the Logistic Function
% p(1)=K (carrying capacity), p(2)=r (growth rate), p(3)=P0 (initial)
logistic = @(t, p) (p(1)*p(3)* exp(p(2)*t))./(p(1)+p(3)*(exp(p(2)*t)-1))

% Step 7: Resulting plots

figure(3)
plot(x-4,y,'marker','+','color','black','linewidth',1);
hold on
plot(x_eul,y_eul,'color','red','linewidth',1);
hold on
plot([0:0.1:28]',logistic([0:0.1:28]', [1.278*10^8,0.3,1.2598*10^8]'),'linestyle','--','linewidth',1)
hold on
title('Example 2 - Total Population of Japan (Cubic spline)','fontsize',14);
legend({'Data used','Calculation','Logistic Model'},'location','northeast','fontsize',10)
xlabel('Time (Year)','fontsize',10);
ylabel('Total Population of Japan (In Billions)','fontsize',10);
xticks([-4:4:28]);
xticklabels({'1992','1996','2000','2004','2008','2012','2016','2020','2024'});

% Step 8: Relative Error Plot

x_new   = [0:4:32]';
y_new   = [y; japan_pop(end-5); japan_pop(end-1)];

err_rel = abs(y_new(2:1:end)-y_eul(1:4/h_eul:end-4/h_eul))./y_new(2:1:end);

figure(4)
plot([0:4:28]',err_rel,'marker','*','color','blue','linewidth',1)
hold on
title('Example 2 - Relative Error (Cubic spline)','fontsize',14);
xlabel('Time (Year)','fontsize',10);
ylabel('Relative Error','fontsize',10);
xticks([0:4:28]);
xticklabels({'1996','2000','2004','2008','2012','2016','2020','2024'});
yticks([0:0.004:0.02]);
yticklabels({'0%','0.4%','0.8%','1.2%','1.6%','2.0'});
