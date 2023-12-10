clear; clear -global;
clc;close all;
%% Load Parameters from the results of the homeostatic optimization!
global k_act kton

k_act = 0.1;
kton = 2;

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
close all;
subepi = load([newdir,'\OptimizationResults\subepi']);
subendo = load([newdir,'\OptimizationResults\subendo']);


mu = subendo.mu;

subendo.tau_h = test_optimization(subendo);

% return
% semilogx(R,p_mid/133.32,'.'); hold on;
%% autoregulation of coronary arteries
% subendo.q(1) = subendo.q(1)*2.5;
subendo.R = subendo.Radius;
subendo.r_act = subendo.R;
subendo.r_h = subendo.R;
subendo.dt = 0.002;
subendo.t = 0;
%    Pin = subendo.Pin*1;
p_out = 55*133.32;%subendo.Pout;%60*133.32;
p_in = p_out:subendo.Pin*0.02:subendo.Pin*1.8;
% for j = 1:length(Pin)
%     PF = myogenic_control(Pin(j) - subendo.Pim , subendo.Pmid(1) - subendo.Pim);
%     scatter(Pin(j)/133.32,PF,'.c');hold on;
% end


for i=1:subendo.N_gen
    if subendo.R(i) < 2.5e-5
        subendo.ID(i) = 'E';
    elseif subendo.R(i) >= 2.5e-5 && subendo.R(i) < 5.0e-5
        subendo.ID(i) = 'D';
    elseif subendo.R(i) >= 5.0e-5 && subendo.R(i) < 9.5e-5
        subendo.ID(i) = 'C';
    elseif subendo.R(i) >= 9.5e-5
        subendo.ID(i) = 'B';
    else
        disp('error in parameter assignment');
    end
end

% control mechanisms activation (0 deactivated, 1 active, 2 maximum constriction): a = total control, af = flow-dependent, ap = myogenic, am = metabloic.
M = 240;
[subendo.q_incr, subendo.gol] = metabolic_control(M);


af = 1;
ap = 1;
am = 1;
a = 1;

marker = '-.';

Autoregulation(subendo, p_in, p_out, af, ap, am, a, marker)

af = 1;
ap = 1;
am = 0;
a = 1;

marker = ':';
Autoregulation(subendo, p_in, p_out, af, ap, am, a, marker)

af = 0;
ap = 1.5;%%% Please check the activation levels and make sure!!!
am = 0;
a = 1;

marker = '--';
Autoregulation(subendo, p_in, p_out, af, ap, am, a, marker)

af = 1;
ap = 1;
am = 1;
a = 0;

marker = '-';
Autoregulation(subendo, p_in, p_out, af, ap, am, a, marker)

legend('All active','myo+flow','max constriction','passive','Location','NorthWest');

xlabel('p (mmHg)');
ylabel('q/q_0');









