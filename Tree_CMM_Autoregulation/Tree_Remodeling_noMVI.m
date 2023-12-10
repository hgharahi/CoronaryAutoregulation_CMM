clear; clear -global;
clc;
close all;
%% Load Parameters from the results of the homeostatic optimization!
global k_act kton

k_act = 1.0;
kton = 2;

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);

subepi = load([newdir,'\OptimizationResults\subepi']);
subendo = load([newdir,'\OptimizationResults\subendo']);

% subepi.Act_lvl = subepi.Act_lvl;
% subendo.Sbasal = subendo.Sbasal*1/1.5;
% subepi.Sbasal = subepi.Sbasal*1/1.5;
% subepi.l_max(:) = 1.4;
% subendo.l_max(:) = 1.4;
% subepi.l_min(:) = 0.4;
% subendo.l_min(:) = .4;
% mu = subendo.mu;
subendo.Cap_res = (subendo.Pout - 20*133.32)/subendo.q(end);

subepi.Cap_res = (subepi.Pout - 20*133.32)/subepi.q(end);
% subendo.tau_h = test_optimization(subendo);
% subepi.tau_h = test_optimization(subepi);
% return
% semilogx(R,p_mid/133.32,'.'); hold on;
%% autoregulation of coronary arteries
% subendo.q(1) = subendo.q(1)*2.5;
subendo.R = subendo.Radius;
subendo.r_act = subendo.R;
subendo.r_h = subendo.R;
subendo.dt = 0.5*6.9444e-04;
subendo.t = 0;

subepi.R = subepi.Radius;
subepi.r_act = subepi.R;
subepi.r_h = subepi.R;
subepi.dt = 0.5*6.9444e-04;
subepi.t = 0;

p_out = 20*133.32;
p_in  = p_out:1/4*p_out:subendo.Pin*1.5;
% p_in = 120.*133.32;
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

for i=1:subepi.N_gen
    if subepi.R(i) < 2.5e-5
        subepi.ID(i) = 'E';
    elseif subepi.R(i) >= 2.5e-5 && subepi.R(i) < 5.0e-5
        subepi.ID(i) = 'D';
    elseif subepi.R(i) >= 5.0e-5 && subepi.R(i) < 9.5e-5
        subepi.ID(i) = 'C';
    elseif subepi.R(i) >= 9.5e-5
        subepi.ID(i) = 'B';
    else
        disp('error in parameter assignment');
    end
end


subendo.tau_h = test_optimization(subendo);
subepi.tau_h = test_optimization(subepi);
% control mechanisms activation (a = 0 deactivated, 1 active, 2 maximum constriction): a = total control, af = flow-dependent, ap = myogenic, am = metabloic.

%% Subendo
% 0.9944    0.1896    0.3333

% 
am = 1.0500;
af = 0.30;
ap = 0.3159;

a = 1;
M = 70;
% 
% % for i = 1:e]
marker1 = '-k';
% 
[LL, q_network1, Pmid1, tau1, A1, R1] = Remodeling4(subendo, p_in, p_out, M, af, ap, am, a);
disp('done!')
figure(2);plot(p_in/133.32, q_network1./subendo.q(1),marker1);hold on


am = 1.0500;
af = 0.30;
ap = 0.3159;

a = 0;
M = 70;
% 
% % for i = 1:e]
marker1 = '--k';
% 
[LL, q_network1, Pmid1, tau1, A1, R1] = Remodeling4(subendo, p_in, p_out, M, af, ap, am, a);
disp('done!')
figure(2);plot(p_in/133.32, q_network1./subendo.q(1),marker1);hold on


am = 1.0500;
af = 0.30;
ap = 0.3159;
subendo.Pim = 0;

a = 0;
M = 70;
% 
% % for i = 1:e]
marker1 = '.-k';
% 
[LL, q_network1, Pmid1, tau1, A1, R1] = Remodeling4(subendo, p_in, p_out, M, af, ap, am, a);
disp('done!')
figure(2);plot(p_in/133.32, q_network1./subendo.q(1),marker1);hold on

%% Subepicardial
marker2 = ':k';

% figure(2);plot(p_in/133.32, q_network2/subepi.q(1)),marker2);hold on
% q_network2./subepi.q(1)
% 
% am = 3.6606;
% af = 0.2327;
% ap = 1.1326;

am = 2.8900;
af = 1.0100;
ap = 0.2266;

a = 1;
M = 70;

[LL, q_network2, Pmid2, tau2, A2, R2] = Remodeling4(subepi, p_in, p_out, M, af, ap, am, a);
disp('done!')

figure(2);plot(p_in/133.32, q_network2./subepi.q(1),marker2);hold on




xlabel('p (mmHg)'); ylabel('q mm^3/s')
legend('Subendo. Autoregulation, MVI','Subendo Fully dilated, MVI','Subendo Fully dilated noMVI','Subepi. Autoregulation','Location','Best');