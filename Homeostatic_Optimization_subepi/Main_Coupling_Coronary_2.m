clear; close all; clc;
% Homeostatic Optimization and Hemodynamics in arterial tree
% set some global parameters (NOT variable!):
global H0 rho_w beta gamma R0 Rmin alpha_t mu Pin Pout N_gen Pim A LA IA SA

% Metabolic costs
alpha_t = 1.0701e+03; %   .8(microMoleATP/(min.g)*1.06(g/cm3)/60(s)*31.8*1000(kJ/molATP)*100/42 (efficiency of OXPHOS), From Paul 1983 (Porcine coronaries)
beta = 160;         % W/m^3, Lindstrom et al (2014)
gamma = 0.00891267681;    % J*s/m^3, Lindstrom et al (2014)

%% Tree Parameters

% initial radius of root vessel [Olufsen 2012]RIA: 0.011/2; 
R0 = 0.1436*1.85E-03; %m
% assumed initial thickness of root vessel (as disscussed with Baek)
H0 = 0.1436*1.93E-04;

% Minimum vesel radius to define terminal vessel (or number of generations)
% [Olufsen 2012]: Rmin=50microns=0.005cm=0.00005mn 
%   - arterioles, not cappilaries 
Rmin = 0.00001; % m

% preliminary maximum number of generations
% >=19
N_gen = 10;



%% Hemodynamics Parameters
% Assumptions for Pulsatile hemodynamics:
%   - Find pulsatile solution from Womersley theory (longitudinally 
%     constrained wall)
%   - Given: input flowaveform and zero oscillatory terminal pressure at
%     the last generation (open-end);
%   - Steps:
%   1) Knowing geometry (R,L,symmetry)and wall properties (E,h or Stiffness 
%   matrix component) compute impedances in each segment;
%   2) Knowing impedance and steady solution, compute the total p and q for
%   each segment;
%   3) Postprocess results
%   4) Do that for the realistic pulmonary vessels parameters;
%   - Terminal reflection coefficient GammaVal
%       0   - matching impedance (Z0=ZT=Zinp), 
%       1   - closed end (ZT>>Z0, e.g. Z0=0, Zinp=0)
%       -1  - open end (Z0>>ZT, e.g. ZT=0)

% Blood properties
mu = 0.0035;    % Pa*s, dynamic viscosity
rho = 1060;     % kg/m^3, blood density
nu = mu/rho;    % m^2/s, dynamic viscosity of blood

% Terminal reflection coefficient
% Hollander-2001: pulmonary circulation is of open-end type reflector
% although in McDonald book (6th ed): pulsatililty remais at capillary 
% level, 5-10% of that in arteries
GammaVal = -1;

%% Flow and Pressure Parameters
% waveform data from [Zambrano 2018]

load PressureCoronary.dat
time = PressureCoronary(:,1);
pressure = PressureCoronary(:,2) * 133.32; 

inputs = ReadInputs('subendocardium');
Pin = inputs(1)*133.32;    Pout = inputs(2)*133.32;	Pim = inputs(3)*133.32;

A = ReadParams('B'); %small artery mechanical parameters   
LA = ReadParams('C'); %large arteriole mechanical parameters   
IA = ReadParams('D'); %intermediate arteriole mechanical parameters   
SA = ReadParams('E'); %small arteriole mechanical parameters   

%-------Step 0: Optimization for initial number of generations-------------
TreeOptimization4;



%% ============================ Post Processing ===========================

% preallocate material and geometry parameters
Thickness = zeros(1,N_gen);
YoungMod_tt = zeros(1,N_gen);
YoungMod_zz = zeros(1,N_gen);
nu_tz = zeros(1,N_gen);
nu_zt = zeros(1,N_gen);
StiffMatrix = zeros(2,2,N_gen);

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

% assign given input mean-flow and terminal mean-pressure
qSteady(1) = q(1); 
pTermSteady(N_gen) = Pout; 

%% Find unstressed radius
Ru = zeros(1,N_gen);
Hu = zeros(1,N_gen);

for kk=1:N_gen    
    
    Pmid(kk) = p_mid(kk);
    Thickness(kk) =(1/(SMCtoCOL(kk)+1)*(1.25)*Mt(kk) + ...
         SMCtoCOL(kk)/(SMCtoCOL(kk)+1)*Mt(kk) + Me(kk))/(0.3*rho_w);
     
    RzeroP(kk)=ZeroP(Me(kk),Mt(kk),Radius(kk),p_mid(kk),0);
    [Ru(kk), Hu(kk)] = UnloadedConfig(Radius(kk),Thickness(kk), Me(kk),Mt(kk), Pmid(kk)); %% Ru and Hu are radius and thickness in unloaded configuration
end


%% update parameters and get steady solution
for k=1:N_gen
    % steady solution
    % pInpSteady(k) = pTermSteady(k)+ HydRes(k)*qSteady(k);
    qSteady(k) = Table(k,2);
    pTermSteady(k) = Table(k,3);
    HydRes(k) = Table(k,4);

    % geometry and material parameters
    % use mid-pressure
    [YoungMod_tt(k),YoungMod_zz(k),nu_tz(k),nu_zt(k),StiffMatrix(:,:,k)] =...
           YoungMod_2(Radius(k), Me(k),Mt(k),1,1,p_mid(k),SMCtoCOL(k));
       
    % for output
    ratio(k) = Hu(k)/(2*Ru(k));
    sigma_h(k) = (Pmid(k))/(2*ratio(k));       % homeostatic stress
    shear(k) = 4*mu*qSteady(k)/(pi*Radius(k)^3);% homeostatic shear stress
end  

% get input mean pressure
pInpSteady(1) = pTermSteady(1)+ HydRes(1)*qSteady(1);
for k=2:N_gen
    pInpSteady(k) = pTermSteady(k-1);
end

% get Murray's law exponent for output
for k=1:N_gen-1
    ksi(k) = 1/(log2(Radius(k))-log2(Radius(k+1))); 
end

% RunPulsatile;
Plotting;

plot(ksi);

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
close all;
save([newdir,'\OptimizationResults\subepi']);