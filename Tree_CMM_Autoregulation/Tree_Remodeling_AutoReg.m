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
p_in  = p_out:1/5*p_out:subendo.Pin*1.5;
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
am = .500;
af = 1.30;
ap = .3159;

a = 1;
M = 70;
% 
% % for i = 1:e]
marker1 = '--k';
% 
[LL, q_network1, Pmid1, tau1, A1, R1] = Remodeling4(subendo, p_in, p_out, M, af, ap, am, a);
disp('done!')
% figure(2);plot(p_in/133.32, q_network1./subendo.q(end),marker1);hold on


%% Subepicardial
marker2 = '-.k';

% figure(2);plot(p_in/133.32, q_network2./subepi.q(1),marker2);hold on
% q_network2./subepi.q(1)
% 
% am = 3.6606;
% af = 0.2327;
% ap = 1.1326;

am = 1.5900;
af = 1.100;
ap = .2266;

a = 1;
M = 70;

[LL, q_network2, Pmid2, tau2, A2, R2] = Remodeling4(subepi, p_in, p_out, M, af, ap, am, a);
disp('done!')

% figure(2);plot(p_in/133.32, q_network2./subepi.q(end),marker2);hold on


figure(2);plot(p_in/133.32, (q_network2 + q_network1)./(subepi.q(end)+subendo.q(end)),'-k','LineWidth',1.5);hold on




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
% 
PF = xlsread('../Autoreg','meta','A1:D8');
% 
P = PF(:,1);
F = PF(:,2)/0.578411232; % normalization by flow @ P = 100;
Fmax = PF(:,3)/0.578411232 - F; % normalization by flow @ P = 100;
Fmin = F - PF(:,4)/0.578411232; % normalization by flow @ P = 100;

h1 = figure(2);errorbar(P, F, Fmin, Fmax,'o');hold on

% p_out = 20*133.32;
% p_in  = min(P):0.5*p_out:max(P);
% f_in = interp1(P,F,p_in','spline'); 
% legend('K_m = 5.6972, K_\tau = 2.6496, K_p=2.63, M/M0=1','K_m = 5.6972, K_\tau = 2.6496, K_p=2.63, M/M0=1.33','Laird et al. 1981, M/M0=1','Laird et al. 1981, M/M0=1.33','Dick et al. 2018','Location','best');

xlabel('p_{in} (mmHg)','FontSize',12); ylabel('$\bar{q}$','Interpreter','latex','FontSize',12);
legend('subendo. + subepi.','Meta analysis on swine data [Dick et al. 2018]','Location','Best','FontSize',12);
grid on;
box on;
%print(h1, '-dpdf', ['AutoReg','.pdf']);
%movefile('AutoReg.pdf','Figs\AutoReg.pdf');
% 
% 
%% some plots

h2 = figure(3);
tau_leg = {};
markers = {'-','--','-.','.-','.','+','x'};

c = 0;
for i=3:3:12
    c = c + 1;
    plot(p_in/133.32,tau1(i,:)/subendo.shear(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = ['D = ',num2str(round(2*subendo.Radius(i)*1e6)),' \mum'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('$\bar{\tau}$','Interpreter','latex','FontSize',16); 
title('subendo.','FontSize',12); 
grid on;
box on;
%print(h2, '-dpdf', ['tau_endo','.pdf']);
%movefile('tau_endo.pdf','Figs\tau_endo.pdf');



h3 = figure(4);
tau_leg = {};
c = 0;
for i=3:3:12
    c = c + 1;
    plot(p_in/133.32,tau2(i,:)/subepi.shear(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = ['D = ',num2str(round(2*subepi.Radius(i)*1e6)),' \mum'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('$\bar{\tau}$','Interpreter','latex','FontSize',16); 
title('subepi.','FontSize',12); 
grid on;
box on;
%print(h3, '-dpdf', ['tau_epi','.pdf']);
%movefile('tau_epi.pdf','Figs\tau_epi.pdf');


h4 = figure(5);
tau_leg = {};

c = 0;
for i=3:3:12
    c = c + 1;
    plot(p_in/133.32,A1(i,:),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = ['D = ',num2str(round(2*subendo.Radius(i)*1e6)),' \mum'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('A','FontSize',12); 
% title('subendo.','FontSize',12); 
grid on;
box on;
%print(h4, '-dpdf', ['A_endo','.pdf']);
%movefile('A_endo.pdf','Figs\A_endo.pdf');

h5 = figure(6);
tau_leg = {};
c = 0;
for i=3:3:12
    c = c + 1;
    plot(p_in/133.32,A2(i,:),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = ['D = ',num2str(round(2*subepi.Radius(i)*1e6)),' \mum'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)'); ylabel('A','FontSize',12); 
% title('epicardial','FontSize',12); 
grid on;
box on;
%print(h5, '-dpdf', ['A_epi','.pdf']);
%movefile('A_epi.pdf','Figs\A_epi.pdf');


h6 = figure(7);
tau_leg = {};
c = 0;
for i=3:3:12
    c = c + 1;
    plot(p_in/133.32,R1(i,:)/subendo.Radius(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = ['D = ',num2str(round(2*subendo.Radius(i)*1e6)),' \mum'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('$\bar{D}$','Interpreter','latex','FontSize',12); 
title('Subendo.','FontSize',12); 
grid on;
box on;
%print(h6, '-dpdf', ['D_endo','.pdf']);
%movefile('D_endo.pdf','Figs\D_endo.pdf');

h7 = figure(8);
tau_leg = {};
c = 0;
for i=3:3:12
    c = c + 1;
    plot(p_in/133.32,R2(i,:)/subepi.Radius(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = ['D = ',num2str(round(2*subepi.Radius(i)*1e6)),' \mum'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('$\bar{D}$','Interpreter','latex','FontSize',12); 
title('Subepi.','FontSize',12); 
grid on;
box on;
%print(h7, '-dpdf', ['D_epi','.pdf']);
%movefile('D_epi.pdf','Figs\D_epi.pdf');

% close all;
%%


h8 = figure(9);plot(p_in(2:end)/133.32, q_network1(2:end)./q_network2(2:end),'-k','LineWidth',1.5);hold on
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('ENDO/EPI','Interpreter','latex','FontSize',12); 
grid on;
ENDOtoEPI = [
24.45852973	0.309772398	0.590142061	0.524911573
27.95563425	0.438671731	0.750815348	0.584260474
31.08448531	0.5544946	0.683465379	0.811298739
36.86636016	0.825402557	0.715115751	1.154222313
49.0138585	0.941038123	0.754123083	1.24785747
83.66930758	1.028168851	0.828167692	1.241498383
];

ax2 = axes('Position',[.6 .2 .25 .25]);
plot(ENDOtoEPI(:,1),ENDOtoEPI(:,4),'^-k','LineWidth',1.5);
axis([10 90 0.5 1.3]);
%print(h8, '-dpdf', ['ENDOEPI','.pdf']);
%movefile('ENDOEPI.pdf','Figs\ENDOEPI.pdf');



