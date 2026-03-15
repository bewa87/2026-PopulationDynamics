% Step 0: Clearing

clear all;
close all;
clc;

% Step 1: Load world data

pkg load statistics;

fid = fopen('Example01-World-Population.csv', 'r');
C   = textscan(fid,'%f32 %f32', 'delimiter', ';');
fclose(fid);

world_year  = cell2mat(C(1));
world_pop   = cell2mat(C(2));

start_ind   = find(world_year==1994);
end_ind     = find(world_year==2016);

% Step 2: Definition of World's population data points

time_int = 2;
year     = [world_year(start_ind):time_int:world_year(end_ind)]';
x        = [0:time_int:(world_year(end_ind)-world_year(start_ind))]';
y        = world_pop(start_ind:time_int:end_ind);  %[5.47 5.81 6.144 6.471 6.801 7.141 7.492]';

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
  plot(y_spl(j):h(j)/m:y_spl(j+1),lin_spl{j}(y_spl(j):h(j)/m:y_spl(j+1)),'color','blue','linewidth',1)
  hold on
  plot(y_spl(j),k_spl(j),'marker','*','color','blue')
  hold on
  number = 1994 + time_int*j;
  title('Example 1 - Dependence of k on P (Linear spline)','fontsize',14);
  xlabel('Population P','fontsize',10);
  ylabel('k(P)','fontsize',10);
endfor
  plot(y_spl(end),k_spl(end),'marker','*','color','blue')

% Step 5: Implicit Eulerian Time-Stepping Method

h_eul    = 0.25;
x_eul    = [0:h_eul:x(end-1)+3*time_int]';
y_eul    = zeros(length(x_eul),1);
y_eul(1) = y(3);

N          = 20;
y_zeros    = zeros(N+1,1);
y_zeros(1) = y_eul(1);

m = 1;
for l = 1:1:(length(x_eul)-1)
    if (x_eul(l) < x(m+1))
      for p = 1:1:N
        y_zeros(p+1) = y_zeros(p) - (h_eul*lin_spl{m}(y_zeros(p))*y_zeros(p)+y_eul(l)-y_zeros(p))/(h_eul*(lin_spl{m}(y_zeros(p))+y_zeros(p)*lin_spl_diff{m}(y_zeros(p)))-1);
      endfor
      y_eul(l+1) = y_zeros(end);
      y_zeros(1) = y_eul(l+1);
    elseif ((x_eul(l) >= x(m+1)) & (m < length(lin_spl)))%(m < length(cub_spl)))
      m          = m+1;
      for p = 1:1:N
        y_zeros(p+1) = y_zeros(p) - (h_eul*lin_spl{m}(y_zeros(p))*y_zeros(p)+y_eul(l)-y_zeros(p))/(h_eul*(lin_spl{m}(y_zeros(p))+y_zeros(p)*lin_spl_diff{m}(y_zeros(p)))-1);
      endfor
      y_eul(l+1) = y_zeros(end);
      y_zeros(1) = y_eul(l+1);
    else
      for p = 1:1:N
        y_zeros(p+1) = y_zeros(p) - (h_eul*lin_spl{length(lin_spl)}(y_zeros(p))*y_zeros(p)+y_eul(l)-y_zeros(p))/(h_eul*(lin_spl{length(lin_spl)}(y_zeros(p))+y_zeros(p)*lin_spl_diff{length(lin_spl_diff)}(y_zeros(p)))-1);
      endfor
      y_eul(l+1) = y_zeros(end);
      y_zeros(1) = y_eul(l+1);
    endif
endfor

figure(3)
plot(x(3:1:end)-2*time_int,y(3:1:end),'marker','+','color','black','linewidth',1)
hold on
plot(x_eul,y_eul,'color','red','linewidth',1)
hold on
% Start: Local Sensitivity
plot([x(3)-2*time_int:0.1:x(3)-1*time_int],(y(3)-5e+07).*exp(min(lin_spl{1}([y(3)-1e+08:1e+07:y(3)+1e+08])).*[x(3)-2*time_int:0.1:x(3)-1*time_int]),'color','black','linestyle','--','linewidth',0.5)
hold on
plot([x(3)-2*time_int:0.1:x(3)-1*time_int],(y(3)+5e+07).*exp(max(lin_spl{1}([y(3)-1e+08:1e+07:y(3)+1e+08])).*[x(3)-2*time_int:0.1:x(3)-1*time_int]),'color','black','linestyle','--','linewidth',0.5)
hold on
A1 = ((y(3)-5e+07)*exp(min(lin_spl{1}([y(3)-1e+08:1e+07:y(3)+1e+08])).*(x(3)-1*time_int)));
B1 = ((y(3)+5e+07)*exp(max(lin_spl{1}([y(3)-1e+08:1e+07:y(3)+1e+08])).*(x(3)-1*time_int)));
for k = 4:1:11
plot([x(k)-2*time_int:0.1:x(k)-1*time_int],A1.*exp(min(lin_spl{k-2}([y(k)-1e+08:1e+07:y(k)+1e+08])).*[x(3)-2*time_int:0.1:x(3)-1*time_int]),'color','black','linestyle','--','linewidth',0.5)
hold on
plot([x(k)-2*time_int:0.1:x(k)-1*time_int],B1.*exp(max(lin_spl{k-2}([y(k)-1e+08:1e+07:y(k)+1e+08])).*[x(3)-2*time_int:0.1:x(3)-1*time_int]),'color','black','linestyle','--','linewidth',0.5)
hold on
A1 = A1*exp(min(lin_spl{k-2}([y(k)-1e+08:1e+07:y(k)+1e+08]))*(x(3)-1*time_int));
B1 = B1*exp(max(lin_spl{k-2}([y(k)-1e+08:1e+07:y(k)+1e+08]))*(x(3)-1*time_int));
endfor
% End: Local Sensitivity
title('Example 1 - Total World Population (Linear spline)','fontsize',14);
legend({'Data used','Calculation','Estimated Sensitivity'},'location','northwest','fontsize',10);
xlabel('Time (Year)','fontsize',10);
ylabel('Total World Population','fontsize',10);
xticks([0:time_int:28]);
xticklabels({'1996','1998','2000','2002','2004','2006','2008','2010','2012','2014','2016','2018','2020','2022','2024'});

% Step 6: Relative Error Plot

y_pop = world_pop(start_ind+2*time_int:time_int:end);

err_rel = abs(y_pop-y_eul(1:time_int/h_eul:end))./y_pop;

figure(4)
plot([0:time_int:x(end)+2*time_int]',err_rel,'marker','*','color','blue','linewidth',1);
hold on
title('Example 1 - Relative Error (Linear spline)','fontsize',14);
xlabel('Time (Year)','fontsize',10);
ylabel('Relative Error','fontsize',10);
xticks([0:time_int:28]);
xticklabels({'1996','1998','2000','2002','2004','2006','2008','2010','2012','2014','2016','2018','2020','2022','2024'});
yticks([0:0.0025:0.0125]);
yticklabels({'0%','0.25%','0.5%','0.75%','1%','1.25%'});

% Step 7: Figure of Allee Threshold/Carrying Capacity

at_cc_thresholds = -B./A;

figure(5)
for j = 1:1:length(at_cc_thresholds)
  plot([time_int*j-time_int:0.1:time_int*j]',at_cc_thresholds(j)*ones(length([time_int*j-time_int:0.1:time_int*j]),1),'color','blue','linewidth',1);
  hold on
endfor
title('Example 1 - Time-varying thresholds (Linear spline)','fontsize',14);
xlabel('Time (Year)','fontsize',10);
ylabel('Threshold size','fontsize',10);
xticks([0:time_int:20]);
xticklabels({'1996','1998','2000','2002','2004','2006','2008','2010','2012','2014','2016'});

##figure(3)
##plot(x_new,y_new,'marker','+','color','black','linewidth',1);
##hold on
##plot(x_eul,y_eul,'marker','*','color','red','linewidth',1);
##hold on
##title('Example 2 - Total Population of Japan (Linear spline)','fontsize',14);
##legend({'Data','Calculation'},'location','northeast','fontsize',10)
##xlabel('Time (Year)','fontsize',10);
##ylabel('Total Population of Japan (In Billions)','fontsize',10);
##xticks([0:4:28]);
##xticklabels({'1996','2000','2004','2008','2012','2016','2020','2024'});
##
##% Step 5: Relative Error Plot
##
##y_pop   = [125.178 125.757 126.400 126.843 127.445 127.761 127.854 128.063 128.070 127.629 127.276 127.076 126.811 126.261 125.125 123.975]';
##
##err_rel = abs(y_pop-y_eul(1:4/h_eul:end))./y_pop;
##
##figure(4)
##plot([0:2:32]',err_rel,'marker','*','color','blue','linewidth',1)
##hold on
##title('Example 2 - Relative Error (Linear spline)','fontsize',14);
##xlabel('Time (Year)','fontsize',10);
##ylabel('Relative Error','fontsize',10);
##xticks([0:4:32]);
##xticklabels({'1996','2000','2004','2008','2012','2016','2020','2024'});
##yticks([0:0.03:0.24]);
##yticklabels({'0%','3%','6%','9%','12%','15%','18%','21%','24%'});
##
##figure(5)
##plot([123:1:130]',lin_spl{4}([123:1:130]'))
