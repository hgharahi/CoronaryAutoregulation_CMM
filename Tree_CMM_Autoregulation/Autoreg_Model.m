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
p_in  = [20:5:160]*133.32;

a = 1;
[q_network1, A1, R1, SMAX1, TAU1, PMID1, S1, SIGMA_H1] = Eval_K(K_est(1:3), subendo);
[q_network2, A2, R2, SMAX2, TAU2, PMID2, S2, SIGMA_H2] = Eval_K(K_est(4:6), subepi);

a = 0;
[q_network1p, A1p, R1p, SMAX1p, TAU1p, PMID1p, S1p] = Eval_K(K_est(1:3), subendo);
[q_network2p, A2p, R2p, SMAX2p, TAU2p, PMID2p, S2p] = Eval_K(K_est(4:6), subepi);

% a = 2;
% [q_network1c, A1c, R1c, SMAX1c, TAU1c, PMID1c, S1c] = Eval_K(K_est(4:6), subendo);
% [q_network2c, A2c, R2c, SMAX2c, TAU2c, PMID2c, S2c] = Eval_K(K_est(4:6), subepi);
%% Plots

h1 = figure;

p1 = errorbar(P/133.32, F, Fmin, Fmax,'ob');hold on
p2 = plot(p_in/133.32, (q_network2 + q_network1)./(subepi.q(end)+subendo.q(end)),'-k','LineWidth',1.5);
p3 = plot(p_in/133.32, (q_network2p + q_network1p)./(subepi.q(end)+subendo.q(end)),'--k','LineWidth',1.5);
% p4 = plot(p_in/133.32, (q_network2c + q_network1c)./(subepi.q(end)+subendo.q(end)),'--k','LineWidth',1.5);

xlabel('p_{in} (mmHg)','FontSize',12); ylabel('$q/q_{h}$','Interpreter','latex','FontSize',12);
legend([p2, p3, p1],'Model fit', 'Fully dilated trees','Swine data [Dick et al. 2018]','Location','Best','FontSize',12);
grid on;
box on;
print(h1, '-dmeta', ['AutoReg','.emf']);
movefile('AutoReg.emf','Figs\AutoReg.emf');

prompt = 'Should I continue with plots? (1:yes, else:no)';
Y = input(prompt);

if Y==1
    ParamEst_plots(q_network1, A1, R1, SMAX1, TAU1, PMID1, S1, p_in, subendo);
    
    ParamEst_plots(q_network2, A2, R2, SMAX2, TAU2, PMID2, S2, p_in, subepi);
    
    h2 = figure;
    plot(p_in/133.32, q_network1./q_network2,'-k','LineWidth',1.5);
    xlabel('p_{in} (mmHg)','FontSize',12); ylabel('ENDO/EPI','Interpreter','latex','FontSize',12);
    grid on;
    box on;
    print(h2, '-dmeta', ['ENDOEPI','.emf']);
    movefile('ENDOEPI.emf','Figs\ENDOEPI.emf');
end