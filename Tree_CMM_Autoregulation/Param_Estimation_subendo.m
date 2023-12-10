clear; clear -global;
clc;
close all;
%% Load Parameters from the results of the homeostatic optimization!
global p_in p_out M a f_in

k_act = 1.0;
kton = 2;

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


p_out = 20*133.32;
p_in  = PF(:,1)*133.32;
f_in = F;%interp1(P,F,p_in,'spline'); 

M = 70;

% K0 = [1, 1, 1];
K0 = rand(1,4)*5;
% K0(2) = 0.0001;
% K0(5) = 0.0001;
% Eval_K(K0, subendo, subepi);

f = @(x)Obj_K(x, subendo, subepi);

options = optimset('MaxIter',250,'TolFun',1e-5,'TolX',1e-5,'Display','iter');

[K,fval,exitflag,output] = fminsearch(f, K0, options);
% 
K

K_est = [K(1) 0 K(2) K(3) 0 K(4)];

a = 1;
[q_network1, A1, R1, SMAX1, TAU1, PMID1, S1] = Eval_K(K_est(1:3), subendo);
[q_network2, A2, R2, SMAX2, TAU2, PMID2, S2] = Eval_K(K_est(4:6), subepi);

a = 0;
[q_network1p, A1p, R1p, SMAX1p, TAU1p, PMID1p, S1p] = Eval_K([0, 0, 0], subendo);
[q_network2p, A2p, R2p, SMAX2p, TAU2p, PMID2p, S2p] = Eval_K([0, 0, 0], subepi);
%% Plots


figure;

errorbar(P/133.32, F, Fmin, Fmax,'ob');hold on
plot(p_in/133.32, (q_network2 + q_network1)./(subepi.q(end)+subendo.q(end)),'--k','LineWidth',1.5);
plot(p_in/133.32, (q_network2p + q_network1p)./(subepi.q(end)+subendo.q(end)),'--k','LineWidth',1.5);


ParamEst_plots(q_network1, A1, R1, SMAX1, TAU1, PMID1, S1, p_in, subendo);

ParamEst_plots(q_network2, A2, R2, SMAX2, TAU2, PMID2, S2, p_in, subepi);

