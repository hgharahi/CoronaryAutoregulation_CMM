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

M = 70;
p_out = 20*133.32;
p_in  = 100*133.32;

a = 1;
[q_network1, A1, R1, SMAX1, TAU1, PMID1, S1] = Adenosine(K_est(1:3), subendo);
[q_network2, A2, R2, SMAX2, TAU2, PMID2, S2] = Adenosine(K_est(4:6), subepi);





%% Plots

[~,sheet_name]=xlsfinfo('LNAMEandADE.xlsx');
for k=1:numel(sheet_name)
  data{k}=xlsread('LNAMEandADE.xlsx',sheet_name{k},'A2:B100');
  eval([sheet_name{k},'=data{',num2str(k),'};']);
end

h1 = figure(1);
plot(2*subendo.Radius*1e6,(R1'-subendo.Radius)./subendo.Radius*100,'--','LineWidth',1.5);hold on
plot(2*subepi.Radius*1e6,(R2'-subepi.Radius)./subepi.Radius*100,'LineWidth',1.5);hold on
scatter(ADE(:,1),ADE(:,2),'^k','filled');
plot([0,1000],[0,0],':k');
axis([0 500 -40 40]);
xlabel('D_{h} (\mum)','FontSize',12); ylabel('% \DeltaD','FontSize',12);
legend('Subendo.','Subepi.','[Jones et al. 1995]','FontSize',12,'Location','NorthWest');
title('Adenosine','FontSize',12);
set(gca,'XDir','reverse');
grid on; box on;
% print(h1, '-dpdf', ['ADE','.pdf']);
% movefile('ADE.pdf','Figs\ADE.pdf');

[q_network1, A1, R1, SMAX1, TAU1, PMID1, S1] = L_NAME(K_est(1:3), subendo);
[q_network2, A2, R2, SMAX2, TAU2, PMID2, S2] = L_NAME(K_est(4:6), subepi);

h2 = figure(2);hold on
% rectangle('Position',[0 -40 50 80],'FaceColor',[0.8 1 .1]);
% alpha(0.1);
% rectangle('Position',[50 -40 50 80],'FaceColor',[1 .8 1]);
% alpha(0.1);
% rectangle('Position',[100 -40 90 80],'FaceColor',[1 1 .8]);
% alpha(0.1);
% rectangle('Position',[190 -40 310 80],'FaceColor',[0.8 .8 1]);
% alpha(0.1);

plot(2*subendo.Radius*1e6,(R1'-subendo.Radius)./subendo.Radius*100,'--','LineWidth',1.5)
plot(2*subepi.Radius*1e6,(R2'-subepi.Radius)./subepi.Radius*100,'LineWidth',1.5);hold on
scatter(LNAME(:,1),LNAME(:,2),'^k','filled');
plot([0,1000],[0,0],':k');
axis([0 500 -45 40]);

set(gca,'XDir','reverse');
xlabel('D_{h} (\mum)','FontSize',12); ylabel('% \DeltaD','FontSize',12);
title('L-NAME','FontSize',12);
legend('Subendo.','Subepi.','[Jones et al. 1995]','FontSize',12,'Location','NorthWest');
grid on; box on;
% print(h2, '-dpdf', ['NOinhibit','.pdf']);
% movefile('NOinhibit.pdf','Figs\NOinhibit.pdf');


