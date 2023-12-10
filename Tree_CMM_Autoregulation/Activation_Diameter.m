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

subepi = load([newdir,'\OptimizationResults\subepi_Umich']);
subendo = load([newdir,'\OptimizationResults\subendo_Umich']);

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
%%
p_out = 20*133.32;
p_in  = 59*133.32;
%% Subendo

% 
am = .500;
af = .30;
ap = .3159;

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
[LL, q_network1, Pmid1, tau1, A1, R1] = Remodeling5(subendo, p_in, p_out, M, af, ap, am, a, Ks0, Kp0, Km0, phi_s, phi_p, phi_m);
disp('done!')
% figure(2);plot(p_in/133.32, q_network1./subendo.q(end),marker1);hold on


%% Subepicardial
marker2 = '-.k';

am = 1.5900;
af = .100;
ap = .2266;

[Ks0, Kp0, Km0, phi_s, phi_p, phi_m, a] =  stimuli_params(subepi, af, ap, am);

% phi_s(:) = 0;
% phi_m(:) = 0;
% phi_p(:) = 0;
% a = 1;
M = 70;

[LL, q_network2, Pmid2, tau2, A2, R2] = Remodeling5(subepi, p_in, p_out, M, af, ap, am, a, Ks0, Kp0, Km0, phi_s, phi_p, phi_m);
disp('done!')




[~,sheet_name]=xlsfinfo('../DiameterChangeAutoregulation.xlsx');
for k=1:numel(sheet_name)
  data{k}=xlsread('../DiameterChangeAutoregulation.xlsx',sheet_name{k},'A2:B100');
  eval([sheet_name{k},'=data{',num2str(k),'};']);
end

h1 = figure(1);
plot(2*subendo.Radius*1e6,(R1'-subendo.Radius)./subendo.Radius*100,'--','LineWidth',1.5);hold on
plot(2*subepi.Radius*1e6,(R2'-subepi.Radius)./subepi.Radius*100,'LineWidth',1.5);hold on
scatter(P59(:,1),P59(:,2),'^k','filled');
plot([0,1000],[0,0],':k');
axis([0 500 -40 40]);
xlabel('D_{100} (\mum)','FontSize',12); ylabel('% \DeltaD','FontSize',12);
legend('Subendo.','Subepi.','[Kanatsuka et al. 1990]','FontSize',12);
title('p_{in}=59 mmHg','FontSize',12);
grid on; box on;
print(h1, '-dmeta', ['Delta_mild','.emf']);



h2 = figure(2);
plot(2*subendo.Radius*1e6,(R1'-subendo.Radius)./subendo.Radius*100,'--','LineWidth',1.5);hold on
plot(2*subepi.Radius*1e6,(R2'-subepi.Radius)./subepi.Radius*100,'LineWidth',1.5);hold on
scatter(P38(:,1),P38(:,2),'^k','filled');
plot([0,1000],[0,0],':k');
axis([0 500 -40 40]);
xlabel('D_{100} (\mum)','FontSize',12); ylabel('% \DeltaD','FontSize',12);
title('p_{in}=38 mmHg','FontSize',12);
legend('Subendo.','Subepi.','[Kanatsuka et al. 1990]','FontSize',12);
grid on; box on;
print(h2, '-dmeta', ['Delta_severe','.emf']);
movefile('Delta_severe.emf','Figs\Delta_severe.emf');


