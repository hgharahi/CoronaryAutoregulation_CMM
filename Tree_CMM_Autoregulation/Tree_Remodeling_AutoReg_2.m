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

subepi = load([newdir,'\OptimizationResults\subepi_2020']);
subendo = load([newdir,'\OptimizationResults\subendo_2020']);


subendo.Cap_res = (subendo.Pout - 20*133.32)/subendo.q(end);

subepi.Cap_res = (subepi.Pout - 20*133.32)/subepi.q(end);

TYPEs = {'Small Artery';'Large Arteriole';'Intermediate Arteriole';'Small Arteriole'};
%% autoregulation of coronary arteries

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
% am = 1.1900;
% af = .30;
% ap = .3388;

am = 6.3;
af = .30;
ap = 0.75;

[Ks0, Kp0, Km0, phi_s, phi_p, phi_m, a] =  stimuli_params(subendo, af, ap, am);

% phi_s(:) = 0;
% phi_m(:) = 0;
% phi_p(:) = 0;
% a = 1;
M = 70;
% 
% % for i = 1:e]
marker1 = '--k';
% 
[LL, q_network1, Pmid1, tau1, A1, T1, Tmax1, R1] = Remodeling5(subendo, p_in, p_out, M, af, ap, am, a, Ks0, Kp0, Km0, phi_s, phi_p, phi_m);

% Full dilation
a = 0*a;
[LL, q_network1_2, Pmid1_2, tau1_2, A1_2, T1_2, Tmax1_2, R1_2] = Remodeling5(subepi, p_in, p_out, M, af, ap, am, a, Ks0, Kp0, Km0, phi_s, phi_p, phi_m);

disp('done!')
% figure(2);plot(p_in/133.32, q_network1./subendo.q(end),marker1);hold on


%% Subepicardial
marker2 = '-.k';

am = 5.3;
af = .33;
ap = 0.68;

% Autoregulation
[Ks0, Kp0, Km0, phi_s, phi_p, phi_m, a] =  stimuli_params(subepi, af, ap, am);

% phi_s(:) = 0;
% phi_m(:) = 0;
% phi_p(:) = 0;
% a = 1;
M = 70;

[LL, q_network2, Pmid2, tau2, A2, T2, Tmax2, R2] = Remodeling5(subepi, p_in, p_out, M, af, ap, am, a, Ks0, Kp0, Km0, phi_s, phi_p, phi_m);

% Full dilation
a = 0*a;
[LL, q_network2_2, Pmid2_2, tau2_2, A2_2, T2_2, Tmax2_2, R2_2] = Remodeling5(subepi, p_in, p_out, M, af, ap, am, a, Ks0, Kp0, Km0, phi_s, phi_p, phi_m);


disp('done!')

% figure(2);plot(p_in/133.32, q_network2./subepi.q(end),marker2);hold on


figure(2);plot(p_in/133.32, (q_network1)./(subendo.q(end)),'-k','LineWidth',1.5);hold on
plot(p_in/133.32, (q_network2_2 + q_network1_2)./(subepi.q(end)+subendo.q(end)),'--k','LineWidth',1.5);hold on



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

h1 = figure(2);errorbar(P, F, Fmin, Fmax,'ob');hold on

% Data from Kiel ... Tune

% CPP = [40 50 60 70 80 90 100 110 120 130 140];
% Q_tune = [0.28 0.38 0.47 0.50 0.56 0.58 0.63 0.67 0.7 0.83 0.95];
% plot(CPP, Q_tune/Q_tune(7));
% p_out = 20*133.32;
% p_in  = min(P):0.5*p_out:max(P);
% f_in = interp1(P,F,p_in','spline'); 
% legend('K_m = 5.6972, K_\tau = 2.6496, K_p=2.63, M/M0=1','K_m = 5.6972, K_\tau = 2.6496, K_p=2.63, M/M0=1.33','Laird et al. 1981, M/M0=1','Laird et al. 1981, M/M0=1.33','Dick et al. 2018','Location','best');

xlabel('p_{in} (mmHg)','FontSize',12); ylabel('$q/q_{h}$','Interpreter','latex','FontSize',12);
legend('Model fit', 'Fully dilated trees','Swine data [Dick et al. 2018]','Location','Best','FontSize',12);
grid on;
box on;
print(h1, '-dmeta', ['AutoReg_FF','.emf']);
%movefile('AutoReg.emf','Figs\AutoReg.emf');
% 
% 
%% some plots

h2 = figure(3);
tau_leg = {};
markers = {'-','--','-.','.-','.','+','x'};

c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,tau1(i,:)/subendo.shear(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = [TYPEs{c},' (',num2str(round(2*subendo.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('$\tau/\tau_{h}$','Interpreter','latex','FontSize',16); 
title('Subendocardial Tree','FontSize',12); 
grid on;
box on;
print(h2, '-dmeta', ['tau_endo','.emf']);
movefile('tau_endo.emf','Figs\tau_endo.emf');



h3 = figure(4);
tau_leg = {};
c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,tau2(i,:)/subepi.shear(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} =  [TYPEs{c},' (',num2str(round(2*subepi.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('$\tau/\tau_{h}$','Interpreter','latex','FontSize',16); 
title('Subepicardial Tree','FontSize',12); 
grid on;
box on;
print(h3, '-dmeta', ['tau_epi','.emf']);
movefile('tau_epi.emf','Figs\tau_epi.emf');


h4 = figure(5);
tau_leg = {};

c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,A1(i,:),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} =  [TYPEs{c},' (',num2str(round(2*subendo.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('A','FontSize',12); 
title('Subendocardial Tree','FontSize',12); 
grid on;
box on;
axis([20 160 0 1]);
print(h4, '-dmeta', ['A_endo','.emf']);
movefile('A_endo.emf','Figs\A_endo.emf');

h5 = figure(6);
tau_leg = {};
c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,A2(i,:),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = [TYPEs{c},' (',num2str(round(2*subepi.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)'); ylabel('A','FontSize',12); 
title('Subepicardial Tree','FontSize',12); 
grid on;
box on;
axis([20 160 0 1]);
print(h5, '-dmeta', ['A_epi','.emf']);
movefile('A_epi.emf','Figs\A_epi.emf');

%%
h4_2 = figure(52);
tau_leg = {};

c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,T1(i,:),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} =  [TYPEs{c},' (',num2str(round(2*subendo.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('S (MPa)','FontSize',12); 
title('Subendocardial Tree','FontSize',12); 
grid on;
box on;
axis([20 160 0 2.2]);
print(h4_2, '-dmeta', ['S_endo','.emf']);
movefile('S_endo.emf','Figs\S_endo.emf');

h5_2 = figure(62);
tau_leg = {};
c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,T2(i,:),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = [TYPEs{c},' (',num2str(round(2*subepi.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)'); ylabel('S (MPa)','FontSize',12); 
title('Subepicardial Tree','FontSize',12); 
grid on;
box on;
axis([20 160 0 2.2]);
print(h5_2, '-dmeta', ['S_epi','.emf']);
movefile('S_epi.emf','Figs\S_epi.emf');

%%

h4_3 = figure(53);
tau_leg = {};

c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,Tmax1(i,:),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} =  [TYPEs{c},' (',num2str(round(2*subendo.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('S_m (MPa)','FontSize',12); 
title('Subendocardial Tree','FontSize',12); 
grid on;
box on;
axis([20 160 0 2.2]);
print(h4_3, '-dmeta', ['Sm_endo','.emf']);
movefile('Sm_endo.emf','Figs\Sm_endo.emf');

h5_3 = figure(63);
tau_leg = {};
c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,Tmax2(i,:),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = [TYPEs{c},' (',num2str(round(2*subepi.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)'); ylabel('S_m (MPa)','FontSize',12); 
title('Subepicardial Tree','FontSize',12); 
grid on;
box on;
axis([20 160 0 2.2]);
print(h5_3, '-dmeta', ['Sm_epi','.emf']);
movefile('Sm_epi.emf','Figs\Sm_epi.emf');

%%
h6 = figure(7);
tau_leg = {};
c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,R1(i,:)/subendo.Radius(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = [TYPEs{c},' (',num2str(round(2*subendo.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('$D/D_{h}$','Interpreter','latex','FontSize',12); 
title('Subendocardial Tree','FontSize',12); 
grid on;
box on;
axis([20 160 0.5 1.3]);
print(h6, '-dmeta', ['D_endo','.emf']);
movefile('D_endo.emf','Figs\D_endo.emf');

h7 = figure(8);
tau_leg = {};
c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,R2(i,:)/subepi.Radius(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} =  [TYPEs{c},' (',num2str(round(2*subepi.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('$D/D_{h}$','Interpreter','latex','FontSize',12); 
title('Subepicardial Tree','FontSize',12); 
grid on;
box on;
axis([20 160 0.5 1.3]);
print(h7, '-dmeta', ['D_epi','.emf']);
movefile('D_epi.emf','Figs\D_epi.emf');

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

% ax2 = axes('Position',[.6 .2 .25 .25]);
% plot(ENDOtoEPI(:,1),ENDOtoEPI(:,4),'^-k','LineWidth',1.5);
axis([20 160 0.5 1.35]);
print(h8, '-dmeta', ['ENDOEPI','.emf']);
movefile('ENDOEPI.emf','Figs\ENDOEPI.emf');


figure(1000);
plot(1:11,T2(1:11,21)/T2(end,21)); hold on;
plot(1:11,A2(1:11,21)/A2(end,21));


%% 

h200 = figure(300);
tau_leg = {};
markers = {'-','--','-.','.-','.','+','x'};

c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,Pmid1(i,:)/subendo.p_mid(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = [TYPEs{c},' (',num2str(round(2*subendo.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('p/p_{h}','FontSize',12); 
title('Subendocardial Tree','FontSize',12); 
grid on;
box on;
print(h200, '-dmeta', ['p_endo','.emf']);
movefile('p_endo.emf','Figs\p_endo.emf');



h300 = figure(400);
tau_leg = {};
c = 0;
for i=1:3:subepi.N_gen
    c = c + 1;
    plot(p_in/133.32,Pmid2(i,:)/subepi.p_mid(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} =  [TYPEs{c},' (',num2str(round(2*subepi.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('p/p_{h}','FontSize',12);  
title('Subepicardial Tree','FontSize',12); 
grid on;
box on;
print(h300, '-dmeta', ['p_epi','.emf']);
movefile('p_epi.emf','Figs\p_epi.emf');