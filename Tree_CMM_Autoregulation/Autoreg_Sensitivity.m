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
p_in  = PF(:,1)*133.32;

f = @(K)Obj_K_Sensitivity(K, subendo, subepi);

Ep = f(K);

for j = 1:length(K)
    
    Kp = K;
    Kp(j) = K(j) + 0.5*K(j);
    
    E = f(Kp);
    Theta(j,1) = 10/(Ep)*abs(E-Ep);

    
    Km = K;
    Km(j) = K(j) - 0.1*K(j);
    
    E = f(Km);
    Theta(j,2) = 10/(Ep)*abs(E-Ep);
    
end



