clear; clear -global;
clc;
close all;
%% Load Parameters from the results of the homeostatic optimization!
global tree p_in p_out M a f_in

k_act = 1.0;
kton = 2;

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);

% subepi = load([newdir,'\OptimizationResults\subepi']);
tree = load([newdir,'\OptimizationResults\subendo_Umich']);
tree.R = tree.Radius;
tree.r_act = tree.R;
tree.r_h = tree.R;
tree.dt = 0.5*6.9444e-04;
tree.t = 0;
for i=1:tree.N_gen
    if tree.R(i) < 2.5e-5
        tree.ID(i) = 'E';
    elseif tree.R(i) >= 2.5e-5 && tree.R(i) < 5.0e-5
        tree.ID(i) = 'D';
    elseif tree.R(i) >= 5.0e-5 && tree.R(i) < 9.5e-5
        tree.ID(i) = 'C';
    elseif tree.R(i) >= 9.5e-5
        tree.ID(i) = 'B';
    else
        disp('error in parameter assignment');
    end
end
tree.Cap_res = (tree.Pout - 20*133.32)/tree.q(end);

a = 1;

tree.tau_h = test_optimization(tree);




PF = xlsread('../Autoreg','meta','A1:D8');
% 
P = PF(:,1)*133.32;
F = PF(:,2)/0.578411232; % normalization by flow @ P = 100;
Fmax = PF(:,3)/0.578411232 - F; % normalization by flow @ P = 100;
Fmin = F - PF(:,4)/0.578411232; % normalization by flow @ P = 100;
% K = [17.1454    2.3641    0.4292];
% Obj_K(K)


p_out = 20*133.32;
p_in  = PF(:,1)*133.32;
f_in = F;%interp1(P,F,p_in,'spline'); 

M = 70;

K0 = [15., 1.0, 0.3];
options = optimset('MaxIter',250,'TolFun',1e-5,'TolX',1e-5,'Display','iter');

[K,fval,exitflag,output] = fminsearch(@Obj_K, K0, options);


