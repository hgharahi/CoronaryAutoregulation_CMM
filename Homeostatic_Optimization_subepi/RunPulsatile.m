%% ======================= RUN PULSATILE FLOW =============================
%   - deformable wall Womersley solution

disp('Start Pulsatile Solutions');

load FlowCoronary.dat;
time = FlowCoronary(:,1);
scaleFlow = q(1)/mean(FlowCoronary(:,2));
flow = scaleFlow * FlowCoronary(:,2); 

% get flow at root vessel in frequency domain
% QnInp(1) - steady flow, QnInp(k>1) oscillatory flow in frequency domain
% Nf-1: don't consider the last point, it's identical to the first point
[QnInp] = FlowRateFrequencyDomain(flow,NumModes,Nf-1);

% for longitudinally constrained vessel (tethered)
Ctt(1:N_gen)=StiffMatrix(1,1,:);

[pInpTime,pTermTime,qInpTime,qTermTime,cn,Alpha,zinp,zterm,zchar,InpImpedance] ...
    = WomersleySolutionSymmetricTree(N_gen,NumModes,Nt,T,rho,...
                             Radius,Length,Thickness,Ctt,GammaVal,QnInp,...
                             qSteady,pInpSteady,pTermSteady);

% real wave speed for 1st frequency (n=2)
c_R2 = zeros(1,N_gen);c_R10 = zeros(1,N_gen);
for k=1:N_gen
    c_R2(k) = 1.0/real(cn(2,k)^(-1));
    c_R10(k) = 1.0/real(cn(10,k)^(-1));
end

% Moens-Korteweg pulse wave velocity
c0_MK = zeros(1,N_gen); 
for k=1:N_gen
    c0_MK(k) = sqrt(Thickness(k)*YoungMod_tt(k)/(2*rho_w*Radius(k)));
end

% Ett pulse wave velocity
c0_tt = zeros(1,N_gen); 
for k=1:N_gen
    c0_tt(k) = sqrt(Thickness(k)*Ctt(k)/(2*rho_w*Radius(k)));
end
%% Estimate delta parameter (Vasilina's draft of CMM-Womersley paper)
% it checks if the longitudinal velocity is much smaller than the pulse
% wave velocity, delta<<1, needed to satisfy long wave approximation
% I estimated based on the maximum longitudinal velocity, computed from
% the input flow
delta1 = zeros(1,N_gen);delta9 = zeros(1,N_gen); deltaMK = zeros(1,N_gen);
for k=1:N_gen
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
