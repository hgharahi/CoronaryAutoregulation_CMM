clear;close all;clc;

%% Geometry (initial value for root vessel optimization)
% and other reference values

% Olufsen 2012: 
%    R.pulm outlet diam=1.21cm=0.0121m -> r=0.00605m=6.05mm,L=0.0575cm
%    MPA outlet diam=2.6cm -> r=0.013m, L=0.045m
% R0 = 0.013; %0.00605; 
R0 = 0.0055; % [Olufsen 2012]RIA: 0.011/2; 
% Banks 1978: 
%   unstressed pulm artery r=6.96mm, h=0.4mm; estimated r=7mm, h=0.1r
% H0 = 0.1*R0;
H0 = 0.07*R0;% as disscussed with Baek;

% Reference values for check-test:
% Olufsen 2012: 
%    R.pulm L=0.0575cm, MPa L=0.045m
%    Rmin=0.005cm=0.00005mn - not cappilary but arterioles
Rmin = 0.00005;
% L0 = 0.045;
L0 = 0.0125; %[Olufsen 2012]RIA
% elastic modulus at pulmonary artery
%   estimated as 3C/4/6=0.093MPa from C=0.75MPa from Roccabiance ATA,
%   assuming poisson's ratio=0.5
%   estimated as 3C/4/5=0.1125MPa=112.5kPa from C=0.75MPa from Roccabiance ATA
%   assuming poisson's ratio=0.5
%   (<30yo) data taken from Hasket-2010
%   [Zambrano 2018] - E=0.5MPa too high, due to H/D=0.019 too low
E0 = 112500;   Ez=500000;
EhrK = 37.5*133.32; %Olufsen 2012, Eh/r0=3/4/Lambda, (Krenz 2003) Lambda=0.02/mmHg
EhrQ = 195*133.32; %Qureshi 2014
EhrY = 3/(4*0.012)*133.32;%62.5*133.32*0.01655/0.012; %Yen1990 (from Qureshi 2014) lambda=0.012 mmHg, HG: average lambda from Yen et al (1990) 0.0146 mmHg
EhZ=100; %Pa*m, from Zambrano h=0.0002m? approximately, MPA
EhrExperiment = [10530.92136,9496.887786,9311.097003,8679.807346,8432.782151]; % Ehr from porcine experiments
EhrEx_mean = mean(EhrExperiment);EhrEx_err = EhrExperiment - EhrEx_mean;

% Banks (1978): Moens-Korteweg PVW for 20-30yo humans c0=2.24m/s
% Milnor (1969): c0=1.68m/s;
PwvB = 2.24; PwvM = 1.68;

%% VF:  Flow and Pressure Parameters (by Vasilina)
% waveform data from Olufsen J. Fluid Mech. (2012),
% in main pulmonary artery, interpolated, scaled to m^3/s=10^6cm^3/s
% load FlowIntOl.dat;
% time = FlowIntOl(:,1);
% flow = FlowIntOl(:,2); 

load FlowBZ.dat; % flow from [Zambrano 2018]
time = FlowBZ(:,1);
flow = FlowBZ(:,2); 

%start from the 4th generation (after RIA in Olufsen 2012) of symm. tree
flow = flow./2^3; %4th gen from MPA
% flow = flow./2^5; % 6th gen from MPA

Nf = size(time,1);
T = time(Nf) - time(1);
Nt = 125;
% limit to number of first modes (filter higher modes 
% even if they have slightly larger amplitude)
NumModes = 10; 

% get flow at root vessel in frequency domain
% QnInp(1) - steady flow, QnInp(k>1) oscillatory flow in frequency domain
% Nf-1: don't consider the last point, it's identical to the first point
[QnInp] = FlowRateFrequencyDomain(flow,NumModes,Nf-1);

% define steady input flow
Q_parent =  mean(flow);

% %test Hamid's code
% Q_parent = (10+10)*0.000001;

% Capilary pressure: 
%   6-8mmHg, Ganter et al (2006), 10665.76dyn/cm2
%   Attinger 1963: 
%        approximetely 0.5-0.6 of the mean pulm. press.(15mmHg)
% p_terminal = 8*133.32; %6-8 mmHg
% Quereshi 2014: mean pressure= 10mmHg at artery-vein boudary set at
%   r=rmin=50micron. Note that for Olufsen2012 ptotal=0 at terminal which
%   is not correct physiologically
p_terminal = 10*133.32; % 11.5 mmHg as mean between arterial and capillary

% Terminal reflection coefficient
% 0 - matching impedance (Z0=ZT=Zinp), 
% 1 - closed end (ZT>>Z0, e.g. Z0=0, Zinp=0)
% -1 - open end (Z0>>ZT, e.g. ZT=0)
% McDonald book (6th ed): 
%   pulasatililty remais at capillary level, 5-10% of that in arteries
GammaVal = -1;  %open-end according to Hollander-2001

% Other refernce values for check-test:
%   normal mean arterial pressure 15, diastolic 10, systolic 25 mmHg;
Pm0=15*133.32;  Pd0=10*133.32; Ps0=25*133.32;

%%
% Tree properties
f=0.5; %0.75; % Reduction in elastin's mass fraction in each generation

%% ============= Parameter study on number of generations=================
% - Optimization of Mass, Geometry and Material parameters
% - VF: Steady solution for the symmetric tree
 
% preliminary max number of generations
N_gen = 2;

% preallocate material and geometry parameters
Thickness = zeros(1,N_gen);
YoungMod_tt = zeros(1,N_gen);
YoungMod_zz = zeros(1,N_gen);
nu_tz = zeros(1,N_gen);
nu_zt = zeros(1,N_gen);
StiffMatrix = zeros(2,2,N_gen);
ksi = zeros(1,N_gen-1);

% preallocate homeostatic values
ratio = zeros(1,N_gen);
sigma_h = zeros(1,N_gen);
shear = zeros(1,N_gen);
Pmid = zeros(1,N_gen);

% preallocate steady state solution
qSteady = zeros(1,N_gen);
pInpSteady = zeros(1,N_gen);
pTermSteady =zeros(1,N_gen);
HydRes = zeros(1,N_gen);


% assign given input mean flow and terminal mean pressure
qSteady(1) = Q_parent; 
pTermSteady(N_gen) = p_terminal; 

 %% Optimization for current number of generations
% [Mt,Me,Radius,Length,Table,p_mid] = TreeOptimization2(Q_parent,p_terminal,phi_e,N_gen); 
[Mt,Me,Radius,Length,Table,p_mid] = TreeCMM(Q_parent,p_terminal,N_gen); 



