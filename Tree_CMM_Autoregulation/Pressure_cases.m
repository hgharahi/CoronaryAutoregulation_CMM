clear; clear -global;
clc;
close all;

global p_in p_out M a f_in

k_act = 1.0;
kton = 2;

K = [10.3389    1.0333    6.9158    0.9795];
K_est = [K(1) 0 K(2) K(3) 0 K(4)];

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);

subendo = load([newdir,'\OptimizationResults\subendo_2020']);
subepi = load([newdir,'\OptimizationResults\subepi_2020']);

subendo = tree_param_set(subendo);

subepi = tree_param_set(subepi);

PF = xlsread('Autoreg','meta','A1:D8');
% 
P = PF(:,1)*133.32;
F = PF(:,2)/0.578411232; % normalization by flow @ P = 100;
Fmax = PF(:,3)/0.578411232 - F; % normalization by flow @ P = 100;
Fmin = F - PF(:,4)/0.578411232; % normalization by flow @ P = 100;
f_in = F;
M = 70;
p_out = 20*133.32;
p_in  = [59, 38]*133.32;

a = 1;
[q_network1, A1, R1, SMAX1, TAU1, PMID1, S1] = Eval_K(K_est(1:3), subendo);
[q_network2, A2, R2, SMAX2, TAU2, PMID2, S2] = Eval_K(K_est(4:6), subepi);

%% Plots

[~,sheet_name]=xlsfinfo('DiameterChangeAutoregulation.xlsx');
for k=1:numel(sheet_name)
  data{k}=xlsread('DiameterChangeAutoregulation.xlsx',sheet_name{k},'A2:B100');
  eval([sheet_name{k},'=data{',num2str(k),'};']);
end

h1 = figure(1);
plot(2*subendo.Radius*1e6,(R1(:,1)'-subendo.Radius)./subendo.Radius*100,'--','LineWidth',1.5);hold on
plot(2*subepi.Radius*1e6,(R2(:,1)'-subepi.Radius)./subepi.Radius*100,'LineWidth',1.5);hold on
scatter(P59(:,1),P59(:,2),'^k','filled');
plot([0,1000],[0,0],':k');
axis([0 500 -40 40]);set(gca,'XDir','reverse');
xlabel('D_h (\mum)','FontSize',12); ylabel('% \DeltaD','FontSize',12);
legend('Subendo.','Subepi.','[Kanatsuka et al. 1990]','FontSize',12,'Location','NorthWest');
title('Mild pressure reduction (p_{in} = 59 mmHg)','FontSize',12);
grid on; box on;
print(h1, '-dpdf', ['Delta_Mild','.pdf']);
movefile('Delta_Mild.pdf','Figs\Delta_Mild.pdf');



h2 = figure(2);
plot(2*subendo.Radius*1e6,(R1(:,2)'-subendo.Radius)./subendo.Radius*100,'--','LineWidth',1.5);hold on
plot(2*subepi.Radius*1e6,(R2(:,2)'-subepi.Radius)./subepi.Radius*100,'LineWidth',1.5);hold on
scatter(P38(:,1),P38(:,2),'^k','filled');
plot([0,1000],[0,0],':k');
axis([0 500 -40 40]);set(gca,'XDir','reverse');
xlabel('D_h (\mum)','FontSize',12); ylabel('% \DeltaD','FontSize',12);
legend('Subendo.','Subepi.','[Kanatsuka et al. 1990]','FontSize',12,'Location','NorthWest');
title('Severe pressure reduction (p_{in} = 38 mmHg)','FontSize',12);
grid on; box on;
print(h2, '-dpdf', ['Delta_Severe','.pdf']);
movefile('Delta_Severe.pdf','Figs\Delta_Severe.pdf');
