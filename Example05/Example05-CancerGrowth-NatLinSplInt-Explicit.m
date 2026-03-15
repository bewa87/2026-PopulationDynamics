% Step 0: Clearing

clear all;
close all;
clc;

% Step 1: Load Cancer Growth Data

pkg load statistics;

fid = fopen('Example05-CancerGrowth.csv','r');
C   = textscan(fid,'%f32 %f32', 'delimiter', ';');
fclose(fid);

cancer_day  = cell2mat(C(1));
cancer_vol  = cell2mat(C(2));

start_ind   = find(cancer_day==0);
end_ind     = find(cancer_day==19);

x           = cancer_day(start_ind:1:end_ind);
y           = cancer_vol(start_ind:1:end_ind);

% Step 2: Computation of k in dependence of y

k    = (1./diff(x)).*(log(y(2:1:end)./y(1:1:(end-1))));

figure(1)
plot(y(2:1:end),k,'marker','*')

% Step 3: Spline interpolation

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
endfor
  plot(y_spl(end),k_spl(end),'marker','*','color','blue')
  hold on
  xlabel('Total volume fraction','fontsize',10)
  ylabel('k','fontsize',10)
  title('Example 5: Spline interpolation (Linear spline)','fontsize',14)

% Step 5: Explicit Eulerian Time-Stepping Method

h_eul    = (1/25);
x_eul    = [x(2):h_eul:x(end)+5]';
y_eul    = zeros(length(x_eul),1);
y_eul(1) = y(2);

m = 1;
for l = 1:1:(length(x_eul)-1)
    if (x_eul(l) < x(m+1))
      y_eul(l+1) = y_eul(l) + h_eul*(lin_spl{m}(y_eul(l)))*y_eul(l);
    elseif ((x_eul(l) >= x(m+1)) & (m < length(lin_spl)))
      m          = m+1;
      y_eul(l+1) = y_eul(l) + h_eul*(lin_spl{m}(y_eul(l)))*y_eul(l);
    else
      y_eul(l+1) = y_eul(l) + h_eul*(lin_spl{length(lin_spl)}(y_eul(l)))*y_eul(l);
    endif
endfor

figure(3)
plot(x,y,'marker','+','color','black','linewidth',1);
hold on
plot(x_eul,y_eul,'color','red','linewidth',1);
hold on
title('Example 5 - Total Volume Fraction (Linear spline)','fontsize',14);
legend({'Data used','Calculation'},'location','northwest','fontsize',10)
xlabel('Time (Days)','fontsize',10);
ylabel('Total Volume Fraction','fontsize',10);
xticks([0:4:32]);
xticklabels({'0','4','8','12','16','20','24','28','32'});

% Step 6: Figure of Allee Threshold/Carrying Capacity

at_cc_thresholds = -B./A;

figure(4)
for j = 1:1:length(y)-2
  plot([x(j):0.1:x(j+1)]',at_cc_thresholds(j)*ones(length([x(j):0.1:x(j+1)]),1),'color','blue','linewidth',1)
  hold on
  plot(x(j),at_cc_thresholds(j),'marker','*','color','blue')
  hold on
endfor
plot(x(length(y)-1),at_cc_thresholds(end),'marker','*','color','blue')
title('Example 5 - Time-varying thresholds (Linear spline)','fontsize',14);
xlabel('Time (Days)','fontsize',10);
ylabel('Threshold size','fontsize',10);
xticks([0:4:20]);
ylim([-10 170])
xticklabels({'0','4','8','12','16','20'});
