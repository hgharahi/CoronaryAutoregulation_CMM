clear; clear -global;
% clc;close all;
%% Load Parameters from the results of the homeostatic optimization!
global k_act kton

k_act = 0.1;
kton = 2;

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
% close all;
% subepi = load([newdir,'\OptimizationResults\subepi']);
subendo = load([newdir,'\OptimizationResults\subendo']);


mu = subendo.mu;

R = subendo.Radius;
r_act = R;
r_h = R;

%    Pin = subendo.Pin*1;

i = 6;
subendo.ratio(i)
%% R100
x = [subendo.c1(i), subendo.c2(i), subendo.c3(i), subendo.c4(i), subendo.c5(i), ...
    subendo.Ge1(i), subendo.Ge2(i), subendo.Gm(i), subendo.Gc(i), subendo.angle(i), ...
    subendo.Act_lvl(i)*subendo.Smax(i), subendo.l_min(i), subendo.l_max(i)];

Me = subendo.Me(i);
Mm = subendo.Mt(i)*(subendo.SMCtoCOL(i)/(subendo.SMCtoCOL(i)+1));
Mc = 1.25*subendo.Mt(i)*(1/(subendo.SMCtoCOL(i)+1));
Mk = [0.1, 0.1, 0.4, 0.4]*Mc;


hm = 1/subendo.rho_w;% subendo.Thickness(i)*Mm/(Mm+Mc+Me)*3/10;

P = subendo.Pmid(i) - subendo.Pim;
R1 = NR_iterate_fit(Me,Mk,Mm,subendo.Radius(i),hm,P,x, subendo.Radius(i), subendo.Radius(i));
R1/R(i)
R100 = NR_iterate_fit(Me,Mk,Mm,subendo.Radius(i),hm,100*133.32,x, r_act(i),r_h(i));
LL = R100/R(i);

%%
if R(i) < 2.5e-5
    ID = 'E';
elseif R(i) >= 2.5e-5 && R(i) < 5.0e-5
    ID = 'D';
elseif R(i) >= 5.0e-5 && R(i) < 9.5e-5
    ID = 'C';
elseif R(i) >= 9.5e-5
    ID = 'B';
else
    disp('error in parameter assignment');
end
% ID
%% Data
data_active = xlsread('../myogenic_response','A2:A8')*0.73555912101486; % Converted from cmH2O to mmHg
data_active(:,2) = xlsread('../myogenic_response',[ID,'2:',ID,'8']).*R100;

R0 = subendo.Ru(i);
Data_1 = load('../TestDataFromAnnemieketalMicroCoron/ParamEstimation_Active/Passive1.3_micro.txt'); % lz = 1.3 van Andel et al. 2003
Data_1(:,2) = Data_1(:,2)./(110.4679/2)*(R0*1e+6)/2;
Data = Data_1; %% This is to scale to the desired radius!

Pin = subendo.Pim:subendo.Pin*0.05:subendo.Pin*1.7;
% Pin = 13.0*133.32 + subendo.Pim;
dt = 0.002;
t = 0;

% figure;
% for j = 1:length(Pin)
%     PF = myogenic_control_subepi(100000000000000000000 , subendo.Pmid(1) - subendo.Pim, ID);
%     scatter((Pin(j)-subendo.Pim)/133.32,PF,'.c');hold on;
% end

% figure;
% 

a = 0;
Pressure_Loop;

a = 1;
Pressure_Loop;
j
% PF = myogenic_control(Pin(end) - subendo.Pim , subendo.Pmid(i) - subendo.Pim, ID);
PF0 = myogenic_control(subendo.Pmid(i) - subendo.Pim , subendo.Pmid(i) - subendo.Pim, ID);

x0 = mechanical_params(subendo, i, PF0*0.51, ID);

RR = NR_iterate_fit(Me,Mk,Mm,subendo.Radius(i),hm,subendo.Pmid(i) - subendo.Pim,x0, subendo.Radius(i), subendo.Radius(i));
RR/subendo.Radius(i)

scatter(data_active(:,2)/subendo.Radius(i),(data_active(:,1)),'.k');
scatter(Data_1(:,2)*1e-6/subendo.Radius(i),(Data_1(:,1))/133.32,'.b');

ID




