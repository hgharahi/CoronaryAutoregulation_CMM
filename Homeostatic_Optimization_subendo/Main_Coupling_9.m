clear; close all; clc;
% Homeostatic Optimization and Hemodynamics in arterial tree
%   - Extension of Murray's Law for homeostasis (steady flow)
%   - Symmetric bifurcating tree
%   - Pulsitile flow postprocessing as for Deformable wall Womersley theory 
%       (longitudinally constrained)
%   - Axisymmetric straight cylindrical vessels

% code by Hamid and Vasilina on August 2 2018

% Important Note: 
%   - System of Units (m,s,kg)
%   - mmHg = 1,333.2 dyn/cm2 = 1,33.32 Pa

%% Wall Tissue Parameters
% Assumptions: 
%   1) we have 4 collagen fiber families in 0, 90, 45, and -45 degrees. 
%   2) SMCs that are circumferentially oriented. 
%   3) isotropic elastin matrix.
%   4) The mass fraction in generations is known.
%   5) Ratio of collagen to SMC is constant throughout the arterial tree. 
%   6) The deposition stretch of SMC and Collagen fibers are constant.
%   8) Stiffness components are computed using Small on Large theory (SoL).
%   9) The metabolic cost of maintenance of collagen and SMCs is the same 
%   and constant. 
%   7) Blood viscosity and metabolic cost (of ?) is constant. 

% set some global parameters (NOT variable!):
global H0 rho_w beta gamma R0 Rmin alpha_t mu lmax lmin S ...
       Ghe1 Ghe2 Ghc Ghm fiber_ang

% Metabolic costs
alpha_t = 2*746;    % W/m^3, Taber 1998 (Biophys J.), porcine carotid artery
beta = 160;         % W/m^3, Lindstrom et al (2014)
gamma = 0.00891267681;    % J*s/m^3, Lindstrom et al (2014)

% Active wall
S =  2.0e+004;
lmax = 1.2;     %  nondim., maximumu tension in active tone
lmin = 0.7;     %  nondim., minimum tension in active tone

% Wall density
rho_w = 1060;   % kg/m3 

% Mechanical_properties (see function mechanical_properties_PA)
Ghe1 = 1.27;
Ghe2 = 1.27;
fiber_ang = 45*pi/180; % angle of helical collagen fibers
Ghc = 1.154; %HG changed from 1.034. June 22nd
Ghm = 1.21;

%% Tree Parameters
% Reduction in elastin's mass fraction in each generation
f = 0.5; 

% initial radius of root vessel [Olufsen 2012]RIA: 0.011/2; 
R0 = 0.0055; %m
% assumed initial thickness of root vessel (as disscussed with Baek)
H0 = 0.07*R0;

% Minimum vesel radius to define terminal vessel (or number of generations)
% [Olufsen 2012]: Rmin=50microns=0.005cm=0.00005mn 
%   - arterioles, not cappilaries 
Rmin = 0.00005; % m

% preliminary maximum number of generations
% >=19
N_gen = 19;

%% Reference Values for Results Comparisons
L0 = 0.0125; %[Olufsen 2012]RIA

% Stiffness parameter E11*H*Rp0 = 3/4/Lambda for isotropic incompressible
%   - Lanbda - distensibility parameter lambda
%   - Rp0 - radius at zero pressure
% (Krenz 2003) Lambda=0.02/mmHg
% (Yen1990 according to Qureshi 2014) lambda=0.012 mmHg
%   - HG: average lambda from Yen et al (1990) 0.0146 mmHg
EhrY = 3/(4*0.012)*133.32;

% Ehr from porcine experiments
EhrExperiment = [10530.92136,9496.887786,9311.097003,8679.807346,8432.782151]; 
EhrEx_mean = mean(EhrExperiment);EhrEx_err = EhrExperiment - EhrEx_mean;

% questionable assumptions
EhrK = 37.5*133.32; %Olufsen 2012, Eh/r0=3/4/Lambda, 
EhrQ = 195*133.32;  %Qureshi 2014 

% Elastic modulus at pulmonary artery
%   as 3C/4/5=0.1125MPa=112.5kPa (from C=0.75MPa from Roccabiance ATA)
%   assuming poisson's ratio=0.5
E0 = 112500;

% wave speed - pulse wave velocity
% Banks (1978): Moens-Korteweg PVW for 20-30yo humans c0=2.24m/s
% Milnor (1969): c0=1.68m/s;
PwvB = 2.24; PwvM = 1.68;

% Normal mean arterial pressure 15, diastolic 10, systolic 25 mmHg;
Pm0=15*133.32;  Pd0=10*133.32; Ps0=25*133.32;

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

% Terminal mean pressure
% Quereshi 2014: mean pressure = 10mmHg at artery-vein boundary set at
%   Rmin=50micron
p_terminal = 10*133.32; % Pa, or 10 mmHg

% Terminal reflection coefficient
% Hollander-2001: pulmonary circulation is of open-end type reflector
% although in McDonald book (6th ed): pulsatililty remais at capillary 
% level, 5-10% of that in arteries
GammaVal = -1;

%% Flow and Pressure Parameters
% waveform data from [Zambrano 2018]
load FlowBZ.dat;
time = FlowBZ(:,1);
flow = FlowBZ(:,2); 

% start from the 4th generation (after RIA in Olufsen 2012) of symm. tree
flow = flow./2^3; %4th gen from MPA

% define steady input flow
q_parent =  mean(flow);

% time steps
Nf = size(time,1);
T = time(Nf) - time(1);
Nt = 125;

% frequency modes, filter higher modes
NumModes = 10; 

%% =================== RUN HOMEOSTATIC OPTIMIZATION =======================
% - Optimization of Mass, Geometry and Material parameters
% - Steady solution for the symmetric tree
%   - Steady flow and mid-pressure for Mass Optimization
%   - R>Rmin, pTermSteady(Newgen) is given
%   - Mt is the mass of Collagen and SMC, Me is the mass of elastin 
%   - SMCtoCOL is ratio of content of smc to collagen
%   - Mtotal = Mc + Ms + Me; Mtotal = Mt*(5/4+SMCtoCOL)/(SMCtoCOL+1) + Me (VF?)

%-------Step 0: Optimization for initial number of generations-------------
disp('Start Homeostatic Optimization');
[Mt,Me,SMCtoCOL,Radius,Length,Table,p_mid]=...
                            TreeOptimization2(q_parent,p_terminal,N_gen); 

% output                        
disp(['    N_gen=',num2str(N_gen),': Radius=',num2str(Radius(N_gen)),...
    ', pTermSteady=',num2str(Table(N_gen,3))]);
% quick verification test
if Table(N_gen,3)~=p_terminal
    disp('ERROR: Step 0 - terminal pressure is not assigned correctly')
end

% find the number of generations related to Rmin
kk=1;
while (kk<=N_gen)&&(Radius(kk)>Rmin)
    kk=kk+1;
end
if (kk==N_gen+1) % initial number of gen. gets R>Rmin for all gen.
    Newgen=N_gen;
else
    Newgen=kk-1;
end

%--------Step 1: for new number of generations rerun optimization----------
% VF warning: may need iterative process to find number of generations
if N_gen~=Newgen
    [Mt,Me,SMCtoCOL,Radius,Length,Table,p_mid]=...
                            TreeOptimization2(q_parent,p_terminal,Newgen); 
    % output 
    disp(['Newgen=',num2str(Newgen),': Radius=',num2str(Radius(Newgen)),...
        ', pTermSteady=',num2str(Table(Newgen,3))]);
    % quick verification test
    if Table(Newgen,3)~=p_terminal
        disp('ERROR: Step 1 - terminal pressure is not assigned correctly')
    end
end
disp('End Homeostatic Optimization');

% preallocate material and geometry parameters
Thickness = zeros(1,Newgen);
YoungMod_tt = zeros(1,Newgen);
YoungMod_zz = zeros(1,Newgen);
nu_tz = zeros(1,Newgen);
nu_zt = zeros(1,Newgen);
StiffMatrix = zeros(2,2,Newgen);
ksi = zeros(1,Newgen-1);

% preallocate homeostatic values
ratio = zeros(1,Newgen);
sigma_h = zeros(1,Newgen);
shear = zeros(1,Newgen);
Pmid = zeros(1,Newgen);

% preallocate steady state solution
qSteady = zeros(1,Newgen);
pInpSteady = zeros(1,Newgen);
pTermSteady =zeros(1,Newgen);
HydRes = zeros(1,Newgen);

% assign given input mean-flow and terminal mean-pressure
qSteady(1) = q_parent; 
pTermSteady(Newgen) = p_terminal; 

%% Find unstressed radius
RzeroP = zeros(1,Newgen);
for kk=1:Newgen
    RzeroP(kk)=ZeroP(Me(kk),Mt(kk),Radius(kk),p_mid(kk),0);
end

%% update parameters and get steady solution
for k=1:Newgen
    % steady solution
    % pInpSteady(k) = pTermSteady(k)+ HydRes(k)*qSteady(k);
    qSteady(k) = Table(k,2);
    pTermSteady(k) = Table(k,3);
    HydRes(k) = Table(k,4);

    % geometry and material parameters
    % use mid-pressure
    [YoungMod_tt(k),YoungMod_zz(k),nu_tz(k),nu_zt(k),StiffMatrix(:,:,k)] =...
           YoungMod_2(Radius(k), Me(k),Mt(k),1,1,p_mid(k),SMCtoCOL(k));
       
    Thickness(k) =(1/(SMCtoCOL(k)+1)*(1.25)*Mt(k) + ...
         SMCtoCOL(k)/(SMCtoCOL(k)+1)*Mt(k) + Me(k))/(0.3*rho_w);

    % for output
    Pmid(k) = p_mid(k);
    ratio(k) = Thickness(k)/(2*Radius(k));
    sigma_h(k) = Pmid(k)/(2*ratio(k));       % homeostatic stress
    shear(k) = 4*mu*qSteady(k)/(pi*Radius(k)^3);% homeostatic shear stress
end  

% get input mean pressure
pInpSteady(1) = pTermSteady(1)+ HydRes(1)*qSteady(1);
for k=2:Newgen
    pInpSteady(k) = pTermSteady(k-1);
end

% get Murray's law exponent for output
for k=1:Newgen-1
    ksi(k) = 1/(log2(Radius(k))-log2(Radius(k+1))); 
end

%% ======================= RUN PULSATILE FLOW =============================
%   - deformable wall Womersley solution

disp('Start Pulsatile Solutions');

% get flow at root vessel in frequency domain
% QnInp(1) - steady flow, QnInp(k>1) oscillatory flow in frequency domain
% Nf-1: don't consider the last point, it's identical to the first point
[QnInp] = FlowRateFrequencyDomain(flow,NumModes,Nf-1);

% for longitudinally constrained vessel (tethered)
Ctt(1:Newgen)=StiffMatrix(1,1,:);

[pInpTime,pTermTime,qInpTime,qTermTime,cn,Alpha,zinp,zterm,zchar,InpImpedance] ...
    = WomersleySolutionSymmetricTree(Newgen,NumModes,Nt,T,rho,...
                             Radius,Length,Thickness,Ctt,GammaVal,QnInp,...
                             qSteady,pInpSteady,pTermSteady);

% real wave speed for 1st frequency (n=2)
c_R2 = zeros(1,Newgen);c_R10 = zeros(1,Newgen);
for k=1:Newgen
    c_R2(k) = 1.0/real(cn(2,k)^(-1));
    c_R10(k) = 1.0/real(cn(10,k)^(-1));
end

% Moens-Korteweg pulse wave velocity
c0_MK = zeros(1,Newgen); 
for k=1:Newgen
    c0_MK(k) = sqrt(Thickness(k)*YoungMod_tt(k)/(2*rho_w*Radius(k)));
end

% Ett pulse wave velocity
c0_tt = zeros(1,Newgen); 
for k=1:Newgen
    c0_tt(k) = sqrt(Thickness(k)*Ctt(k)/(2*rho_w*Radius(k)));
end
%% Estimate delta parameter (Vasilina's draft of CMM-Womersley paper)
% it checks if the longitudinal velocity is much smaller than the pulse
% wave velocity, delta<<1, needed to satisfy long wave approximation
% I estimated based on the maximum longitudinal velocity, computed from
% the input flow
delta1 = zeros(1,Newgen);delta9 = zeros(1,Newgen); deltaMK = zeros(1,Newgen);
for k=1:Newgen
    maxqpulse = max(qInpTime(:,k)) - qSteady(k);
    delta1(k) = maxqpulse/(pi*Radius(k)^2*c_R2(k));
    delta9(k) = maxqpulse/(pi*Radius(k)^2*c_R10(k));
    deltaMK(k) = maxqpulse/(pi*Radius(k)^2*c0_MK(k));
end

%% Get Root Impedance
[Zn,z] = GetImpedance (qInpTime(:,1),pInpTime(:,1),T,Nf,Nt,NumModes);
% Alternatively, can compute as:
% Zn(1)=InpImpedance(1,1) + pInpSteady(1)./qSteady(1)
% Zn(2:NumMode)=InpImpedance(2:NumMode,1)

% total for the root vessel zinp,zterm,zchar,InpImpedance
zinp_total=pInpSteady(1)./qSteady(1) + zinp(:,1);
zterm_total=pTermSteady(1)./qSteady(1) + zterm(:,1);

disp('End Pulsatile Solutions');

%% Output results
Plotting;

