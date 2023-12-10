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
subendo.Cap_res = (subendo.Pout - 30*133.32)/subendo.q(end);

subepi.Cap_res = (subepi.Pout - 30*133.32)/subepi.q(end);
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

p_out = 30*133.32;
% p_in  = p_out:1/6*p_out:subendo.Pin*2.0;
p_in = 100.*133.32;
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
am = 5.6972;
af = 2.6496;
ap = 2.632;

a = 0;
M = 70;

% for i = 1:e]
marker1 = '-';
% figure(1);
[LL, q_network1, Pmid1, tau, A, R1] = Remodeling4(subendo, p_in, p_out, M, af, ap, am, a);
disp('done!')
% figure(2);plot(p_in/133.32, q_network1./subendo.q(1),marker1);hold on


am = 5.6972;
af = 2.6496;
ap = 2.632;

a = 1;
M = 70;

% for i = 1:e]
marker1 = '-k';
% figure(1);
[LL, q_network2, Pmid2, tau, A, R2] = Remodeling4(subendo, p_in, p_out, M, af, ap, am, a);
disp('done!')
% figure(2);plot(p_in/133.32, q_network1./subendo.q(1),marker1);hold on
% q_network1./subendo.q(1)

am = 5.6972;
af = 2.6496;
ap = 2.632;

a = 2;
M = 70;

% for i = 1:e]
marker3 = '-';
% figure(1);
[LL, q_network3, Pmid3, tau, A, R3] = Remodeling4(subendo, p_in, p_out, M, af, ap, am, a);
disp('done!')
% figure(2);plot(p_in/133.32, q_network3./subendo.q(1),marker3);hold on

Rendo = [R1;R2;R3]./subendo.Radius;
Pmid_endo = [Pmid1';Pmid2';Pmid3'];

figure(3); hold on
for i=1:1:subendo.N_gen
    plot(Rendo(:,i),Pmid_endo(:,i)/133.32,marker1)
    scatter(Rendo(1,i),Pmid_endo(1,i)/133.32,'*');
    scatter(Rendo(2,i),Pmid_endo(2,i)/133.32,'s');
    scatter(Rendo(3,i),Pmid_endo(3,i)/133.32,'d');    
end
figure(4);plot([1,2,3],[q_network1,q_network2,q_network3]*1e9,marker1); hold on
names = {'Dilated','Baseline','Constricted'};
set(gca,'xtick',[1:3],'xticklabel',names);
ylabel('q (mL/min)');
%% Subepicardial

am = 15.6972;
af = 2.6496;
ap = 2.632;
a = 0;
M = 70;

% for i = 1:e]
marker2 = '-';
% figure(1);
[LL, q_network4, Pmid4, tau, A, R4] = Remodeling4(subepi, p_in, p_out, M, af, ap, am, a);
disp('done!')
% figure(2);plot(p_in/133.32, q_network2./subepi.q(1),marker2);hold on
% 



am = 5.6972;
af = 2.6496;
ap = 2.632;
a = 1;
M = 70;

% for i = 1:e]
marker2 = '-';
% figure(1);
[LL, q_network5, Pmid5, tau, A, R5] = Remodeling4(subepi, p_in, p_out, M, af, ap, am, a);
disp('done!')
% figure(2);plot(p_in/133.32, q_network2./subepi.q(1),marker2);hold on
% q_network2./subepi.q(1)

am = 5.6972;
af = 2.6496;
ap = 2.632;

a = 2;
M = 70;

% for i = 1:e]
marker3 = '--k';
% figure(1);
[LL, q_network6, Pmid6, tau, A, R6] = Remodeling4(subepi, p_in, p_out, M, af, ap, am, a);
disp('done!')



% figure(2);plot(p_in/133.32, q_network4./subepi.q(1),marker3);hold on
Repi = [R4;R5;R6]./subepi.Radius;
Pmid_epi = [Pmid4';Pmid5';Pmid6'];

figure(3); hold on
for i=1:1:subendo.N_gen
    plot(Repi(:,i),Pmid_epi(:,i)/133.32,marker3)
    scatter(Repi(1,i),Pmid_epi(1,i)/133.32,'*');
    scatter(Repi(2,i),Pmid_epi(2,i)/133.32,'s');
    scatter(Repi(3,i),Pmid_epi(3,i)/133.32,'d');
end
figure(3);xlabel('R/R_{basal}');ylabel('Intraluminal Pressure (mmHg)');
figure(4);plot([1,2,3],[q_network4,q_network5,q_network6]*1e9,marker3); 
ylabel('q_{terminal} (mm^3/s)');

figure(5);plot([1,2,3],[q_network1/q_network4,q_network2/q_network5,q_network3/q_network6],marker3); 
set(gca,'xtick',[1:3],'xticklabel',names);
ylabel('ENDO/EPI');

%%
% PF = xlsread('../Autoreg','Rest','A1:B17');
% 
% P = PF(:,1);
% F = PF(:,2)/130; % normalization by flow @ P = 100;
% 
% figure(2);scatter(P, F,'+b');hold on
% 
% PF = xlsread('../Autoreg','Exercise','A1:B17');
% 
% P = PF(:,1);
% F = PF(:,2)/130; % normalization by flow @ P = 100;
% 
% figure(2);scatter(P, F,'+r');hold on
% % 
% PF = xlsread('../Autoreg','meta','A1:D8');
% % 
% P = PF(:,1);
% F = PF(:,2)/0.578411232; % normalization by flow @ P = 100;
% Fmax = PF(:,3)/0.578411232 - F; % normalization by flow @ P = 100;
% Fmin = F - PF(:,4)/0.578411232; % normalization by flow @ P = 100;
% 
% figure(2);errorbar(P, F, Fmin, Fmax,'ok');hold on
% 
% 
% % legend('K_m = 5.6972, K_\tau = 2.6496, K_p=2.63, M/M0=1','K_m = 5.6972, K_\tau = 2.6496, K_p=2.63, M/M0=1.33','Laird et al. 1981, M/M0=1','Laird et al. 1981, M/M0=1.33','Dick et al. 2018','Location','best');
% 
% xlabel('P (mmHg)'); ylabel('F/F_{100}');
% 
% 
