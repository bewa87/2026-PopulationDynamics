% Step 0: Clearing

clear all;
close all;
clc;

% Step 1: Load world data

pkg load statistics;

fid = fopen('Example03-Germany-Population.csv', 'r');
C   = textscan(fid,'%f32 %f32', 'delimiter', ';');
fclose(fid);

world_year  = cell2mat(C(1));
world_pop   = cell2mat(C(2));

start_ind   = find(world_year==1992);
end_ind     = find(world_year==2016);

% Step 2: Definition of world population data points

time_int = 4;
year     = [world_year(start_ind):time_int:world_year(end_ind)]';
x        = [0:time_int:(world_year(end_ind)-world_year(start_ind))]';
y        = world_pop(start_ind:time_int:end_ind);  %[5.47 5.81 6.144 6.471 6.801 7.141 7.492]';

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

cub_spl      = cell(n,1);
cub_spl_diff = cell(n,1);

for j = 1:1:n
  cub_spl{j}      = @(z) (M(j)*(y_spl(j+1) - z).^3/(6*h(j)))   + (M(j+1)*(z - y_spl(j)).^3/(6*h(j))) + (k_spl(j) - M(j)*h(j)^2/6)*((y_spl(j+1) - z)/h(j)) + (k_spl(j+1) - M(j+1)*h(j)^2/6)*((z - y_spl(j))/h(j));
  cub_spl_diff{j} = @(z) - (M(j)*(y_spl(j+1) - z).^2/(2*h(j))) + (M(j+1)*(z - y_spl(j)).^2/(2*h(j))) - (k_spl(j) - M(j)*h(j)^2/6)*(1/h(j))                + (k_spl(j+1) - M(j+1)*h(j)^2/6)*(1/h(j));
endfor

m = 50;
figure(2)
for j = 1:1:n
  plot(y_spl(j):h(j)/m:y_spl(j+1),cub_spl{j}(y_spl(j):h(j)/m:y_spl(j+1)),'color','black','linewidth',1)
  hold on
  plot(y_spl(j),k_spl(j),'marker','*','color','black')
  hold on
  title('Example 1 - Dependence of k on P (Cubic spline)','fontsize',14);
  xlabel('Population P','fontsize',10);
  ylabel('k(P)','fontsize',10);
endfor

% Step 5: Implicit Eulerian Time-Stepping Method

h_eul    = 2/3;
x_eul    = [0:h_eul:x(end-1)+2*time_int]';
y_eul    = zeros(length(x_eul),1);
y_eul(1) = y(2);

N          = 20;
y_zeros    = zeros(N+1,1);
y_zeros(1) = y_eul(1);

m = 1;
for l = 1:1:(length(x_eul)-1)
    if (x_eul(l) < x(m+1))
      for p = 1:1:N
        y_zeros(p+1) = y_zeros(p) - (h_eul*cub_spl{m}(y_zeros(p))*y_zeros(p)+y_eul(l)-y_zeros(p))/(h_eul*(cub_spl{m}(y_zeros(p))+y_zeros(p)*cub_spl_diff{m}(y_zeros(p)))-1);
      endfor
      y_eul(l+1) = y_zeros(end);
      y_zeros(1) = y_eul(l+1);
    elseif ((x_eul(l) >= x(m+1)) & (m < length(cub_spl)))%(m < length(cub_spl)))
      m          = m+1;
      for p = 1:1:N
        y_zeros(p+1) = y_zeros(p) - (h_eul*cub_spl{m}(y_zeros(p))*y_zeros(p)+y_eul(l)-y_zeros(p))/(h_eul*(cub_spl{m}(y_zeros(p))+y_zeros(p)*cub_spl_diff{m}(y_zeros(p)))-1);
      endfor
      y_eul(l+1) = y_zeros(end);
      y_zeros(1) = y_eul(l+1);
    else
      for p = 1:1:N
        y_zeros(p+1) = y_zeros(p) - (h_eul*cub_spl{length(cub_spl)}(y_zeros(p))*y_zeros(p)+y_eul(l)-y_zeros(p))/(h_eul*(cub_spl{length(cub_spl)}(y_zeros(p))+y_zeros(p)*cub_spl_diff{length(cub_spl_diff)}(y_zeros(p)))-1);
      endfor
      y_eul(l+1) = y_zeros(end);
      y_zeros(1) = y_eul(l+1);
    endif
endfor

figure(3)
plot(x(2:1:end)-time_int,y(2:1:end),'marker','+','color','black','linewidth',1)
hold on
plot(x_eul,y_eul,'color','red','linewidth',1)
hold on
title('Example 1 - Total World Population (Cubic spline)','fontsize',14);
legend({'Data used','Calculation'},'location','northwest','fontsize',10);
xlabel('Time (Year)','fontsize',10);
ylabel('Total World Population (In Billions)','fontsize',10);
xticks([0:time_int:28]);
xticklabels({'1996','2000','2004','2008','2012','2016','2020','2024'});

% Step 6: Relative Error Plot

y_pop = world_pop(start_ind+time_int:time_int:end);

err_rel = abs(y_pop-y_eul(1:time_int/h_eul:end))./y_pop;

figure(4)
plot([0:time_int:28]',err_rel,'marker','*','color','blue','linewidth',1);
hold on
title('Example 1 - Relative Error (Cubic spline)','fontsize',14);
xlabel('Time (Year)','fontsize',10);
ylabel('Relative Error','fontsize',10);
xticks([0:time_int:28]);
xticklabels({'1996','2000','2004','2008','2012','2016','2020','2024'});
yticks([0:0.005:0.025]);
yticklabels({'0%','0.5%','1%','1.5%','2%','2.5%'});
