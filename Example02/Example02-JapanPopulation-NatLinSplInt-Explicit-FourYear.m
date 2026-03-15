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

%year = [1992:4:2016]'; % 2020, 2024
%x    = [0:4:24]';      % 28, 32
%y    = [124.3837 125.9753 127.0278 127.8158 128.1520 127.9283 127.1600]'; % 126.3045 123.7530

% Step 3: Computation of k in dependence of y

k    = (1/time_int).*(log(y(2:1:end)./y(1:1:(end-1)))); %(1./x(2:1:end)).*log(y(2:1:end)/y(1));

figure(1)
plot(y(2:1:end),k,'marker','*')

% Step 4: Spline interpolation

y_spl = y(2:1:end);
k_spl = k(1:1:end);
n = (length(y_spl) - 1);
h = zeros(n,1);
h = diff(y_spl);

[A B] = linear_spline_natural(y_spl,k_spl,n,h);

lin_spl      = cell(n,1);
lin_spl_diff = cell(n,1);

for j = 1:1:n
  lin_spl{j}      = @(z) A(j)*z + B(j);
  lin_spl_diff{j} = @(z) A(j);
endfor

m = 50;
figure(2)
for j = 1:1:n
  plot(y_spl(j):h(j)/m:y_spl(j+1),lin_spl{j}(y_spl(j):h(j)/m:y_spl(j+1)))
  hold on
  plot(y_spl(j),k_spl(j),'marker','*')
  hold on
endfor
  plot(y_spl(end),k_spl(end),'marker','*')

% Step 5: Explicit Eulerian Time-Stepping Method

h_eul    = 0.125;
x_eul    = [0:h_eul:x(end-1)+12]';
y_eul    = zeros(length(x_eul),1);
y_eul(1) = y(2);

m = 1;
for l = 1:1:(length(x_eul)-1)
    if (x_eul(l) < x(m+1))
      y_eul(l+1) = y_eul(l) + h_eul*(lin_spl{m}(y_eul(l)))*y_eul(l);
    elseif ((x_eul(l) >= x(m+1)) & (m < length(lin_spl)))%(m < length(cub_spl)))
      m          = m+1;
      y_eul(l+1) = y_eul(l) + h_eul*(lin_spl{m}(y_eul(l)))*y_eul(l);
    else
      y_eul(l+1) = y_eul(l) + h_eul*(lin_spl{length(lin_spl)}(y_eul(l)))*y_eul(l);
    endif
endfor

figure(3)
plot(x-4,y,'marker','+','color','black','linewidth',1);
hold on
plot(x_eul,y_eul,'marker','*','color','red','linewidth',1);
hold on
title('Example 2 - Total Population of Japan (Linear spline)','fontsize',14);
legend({'Data used','Calculation'},'location','northeast','fontsize',10)
xlabel('Time (Year)','fontsize',10);
ylabel('Total Population of Japan','fontsize',10);
xticks([-4:4:28]);
xticklabels({'1992','1996','2000','2004','2008','2012','2016','2020','2024'});

% Step 6: Relative Error Plot

x_new   = [0:4:32]';
y_new   = [y; japan_pop(end-5); japan_pop(end-1)];

err_rel = abs(y_new(2:1:end)-y_eul(1:4/h_eul:end-4/h_eul))./y_new(2:1:end);

figure(4)
plot([0:4:28]',err_rel,'marker','*','color','blue','linewidth',1)
hold on
title('Example 2 - Relative Error (Linear spline)','fontsize',14);
xlabel('Time (Year)','fontsize',10);
ylabel('Relative Error','fontsize',10);
xticks([0:4:28]);
xticklabels({'1996','2000','2004','2008','2012','2016','2020','2024'});
yticks([0:0.05:0.30]);
yticklabels({'0%','5%','10%','15%','20%','25%','30%'});
