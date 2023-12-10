clear;clc;close all;
%% Parameters
global kc Ghe1 Ghe2 Ghc Ghm H0 rho_w R0 S_basal Lmax Lmin tau_h Mm Mm_h Mc phi_f kton fiber_ang k_act

% Wall mechanical properties (mostly) from Zeinali et al. 2011 (IJSCS)
% (stated otherwise)

%% Simulation params
N_gen = 10;
time_step = 10.0;
dt = 1.0/time_step;
a_max = 350;
final_time = 1500;
num_pa = time_step * a_max+1;
kq_c = 0.1; kq_m = 0.1;
kg = 7; kg_sh = 52;
k_act = 0.1;



%% Hemodynamics params

mu = 0.0035;

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

% define steady input flow
Q_parent =  mean(flow);

p_terminal = 90*133.32; %5 mmHg




%% Wall composition
R0 =  0.0033; % meters! arbitrarily prescribed.
H0 = 0.1 * R0;

rho_w = 1060; % Wall density
phi_f = 0.7; 
Mtotal = (1-phi_f) * rho_w * H0;

nu_e = 0.4; nu_m = 0.3; nu_c = 1 - nu_e - nu_m;

Me = nu_e*Mtotal;
Mc = nu_c*Mtotal;
Mm = nu_m*Mtotal;
Mm_h = Mm;

nu_k = [0.2, 0.2, 0.3, 0.3];
Mk = Mc*nu_k;
%% Mechanical Properties

Ghe1 = 1.25;%
Ghe2 = 1.25;%

fiber_ang = pi/4; % angle of helical collagen fibers
Ghc = 1.07;
%
% c4 = 2.696518979858126e+01;% Pa/kg
% c5 = 8.5;%
Ghm = 1.2;
%
% % Active wall
S_basal = 0.862055*1.75; % N/m (active tone basal tone)
Lmax = 1.2;  % maximumu tension in active tone
Lmin = 0.7; % minimum tension in active tone

% % paremeters
tau_h = 1.5;
tau = tau_h;
kton = 0.1;
sigma_h = p_terminal*R0/H0;
c1 =  sigma_h / ((Ghe1^2-(Ghe1^2*Ghe2^2)^(-1))*(0.3*rho_w));  % material paramter for elastin
c3 = 25.0;    % for collagen fibers
c5 = 8.5;  % for passive smooth muscle
str_ch = 0.0;
ang = [0.0 90.0 45.0 135.0]*(pi/180);
c_frac0 = [.2 .2 .3 .3];
for i=1:4
    str_ch = str_ch + c_frac0(i)*Ghc^2*(Ghc^2-1)*...
        exp(c3*(Ghc^2-1.0)^2)*(sin(ang(i)))^2;
end
c2 = sigma_h/ (str_ch*(0.3*rho_w));
% Active muscle tone (From Biochemomechanical model)
L_act = 1; % no time evolution involved (see Eq 15, Getachew's draft)

[t_act, ~] = active_tone(R0, R0, R0);
hm = Mm/((1-phi_f)*rho_w);
c4 = (sigma_h-t_act/(hm))/(Ghm^2*(Ghm^2-1.0)*...
    exp(c5*(Ghm^2-1.0)^2)*(0.3*rho_w));
kc = [c1, c2, c3, c4, c5];


% age_distibution;

%% Mass production params
homeostatic_calculation;

rho_w = 1060; % Wall density
phi_f = 0.7;

m_basal_c = Mk./a_mean;
m_basal_m = Mm/a_mean;

Me = Me*ones(1,N_gen);
Mc = zeros(1,N_gen);
Mm = zeros(1,N_gen);
Mk = zeros(4,N_gen);
Mk_new = zeros(4,N_gen);
Mtotal = zeros(1,N_gen);

pk1_p = repmat(pk1_p,1,N_gen);
pk2_p = repmat(pk2_p,1,N_gen);
pk3_p = repmat(pk3_p,1,N_gen);
pm_p = repmat(pm_p,1,N_gen);
r_p = repmat(r_p,1,N_gen);

for j=1:N_gen
    
    R(j) = (R0^3/(2^(j-1)))^(1/3);
    r_p(:,j) = R(j)/R0*(r_p(:,j));
    H(j) = 0.1 * R(j);

    Mtotal(j) = (1-phi_f) * rho_w * H(j);
    
    nu_e = 0.4; nu_m = 0.3; nu_c = 1 - nu_e - nu_m;
    
    Me(j) = nu_e*Mtotal(j);
    Mc(j) = nu_c*Mtotal(j);
    Mm(j) = nu_m*Mtotal(j);
    Mm_h(j) = Mm(j);
    
    nu_k = [0.2, 0.2, 0.3, 0.3];
    Mk(:,j) = Mc(j)*nu_k(:);
end





for j=1:N_gen
[Mk(:,j), Mm(:,j), pk1_a(:,j), pk2_a(:,j), pk3_a(:,j), pm_a(:,j)] = age_distribution(pk1_p(:,j), pk2_p(:,j), pk3_p(:,j), pm_p(:,j), m_basal_c, m_basal_m, num_pa, dt, kq_c);
end

Mc = sum(Mk,1);

%% Hemodynamics initialization
% preallocate arrays

R = zeros(1,N_gen);
L = zeros(1,N_gen);
p_mid = zeros(1,N_gen); 
Res = zeros(1,N_gen);   

% initialization
% lhat = ones(1,N_gen);

% intitalize R,L,HydRes and Table
for j=1:N_gen
    R(j) = (R0^3/(2^(j-1)))^(1/3);
    L(j) = LengthSegmentk(R(j));
    Res(j) = 8*mu*(L(j))/(pi*R(j)^4);
end
Table = SymmetricArterialTree2(Q_parent,p_terminal,N_gen,Res);

q_seg(:) = Table(:,2);
Res(:)= Table(:,4);

p_mid(:) = Table(:,3) + (Res(:)/2).*q_seg(:);

%% mass iterations initialization
err = 1000;
it_mass = 1;
c = 1;
t(1) = 0;
Lz = 1.0;
r_act = R;
for j=1:N_gen
    [r(c,j), sigma_k, sigma_m] = NR_iterate(Me(1,j),pk1_a(:,j),pk2_a(:,j),pk3_a(:,j),pm_a(:,j),R(j),p_mid(j),dt,r_p(:,j),r_act(j), j);
    r_act(j) = r(1,j);
end


%% G&R with Hemodynamics update
er = 1000;
while er > 10^-6
    
    t(c+1) = t(c) + dt;
    c = c+1;
    
    for j=1:N_gen
        R(j) =r(c-1,j);
        L(j) = LengthSegmentk(R(j));
        Res(j) = 8*mu*(L(j))/(pi*R(j)^4);
    end
    Table = SymmetricArterialTree2(Q_parent,p_terminal,N_gen,Res);
    
    q_seg(:) = Table(:,2);
    Res(:)= Table(:,4);

    disp([num2str(t(c-1)),'    ',num2str(r(c-1,1)),'    ',num2str(r(c-1,2))]);

    for j=1:N_gen
        
%         j
%         
        tau = 4*mu*q_seg(j)/(pi*r(c-1,j)^3);
        
        p_mid(j) = Table(j,3)+ Res(j)/2*q_seg(j);
        
       
        err = 1000;
        num_it = 1;
        while err > 10^-6 || num_it < 5
            
            [r(c,j), sigma_k, sigma_m, T_act] = NR_iterate( Me(1,j),pk1_a(:,j),pk2_a(:,j),pk3_a(:,j),pm_a(:,j),r(1,j),p_mid(j),dt,r_p(:,j),r_act(j), j);
            Lt = r(c,j)/R0;
            
            for k=1:4
                mPk(k) = f_mP(sigma_k(k), sigma_h, kg, tau, tau_h, kg_sh, m_basal_c(k));
            end
            
            mPm = f_mP(sigma_m, sigma_h, kg, tau, tau_h, kg_sh, m_basal_m);
            
            [Mk_new(:,j), Mm_new(:,j), pk1_a(:,j), pk2_a(:,j), pk3_a(:,j), pm_a(:,j)] = ...
                age_distribution(pk1_p(:,j), pk2_p(:,j), pk3_p(:,j), pm_p(:,j), mPk, mPm, num_pa, dt, kq_c);

            err = sqrt(sum((Mk_new(:,j)-Mk(:,j)).^2)/sum(Mk(:,j).^2) + (Mm_new(:,j)-Mm(:,j))^2/Mm(:,j)^2);
            Mk(:,j) = Mk_new(:,j);
            Mm(:,j) = Mm_new(:,j);
            Mc(:,j) = sum(Mk(:,j));
            total_M = (Me(j)+Mm(j)+Mc(j));
            h(j) = total_M / ((1-phi_f)*rho_w*Lz*Lt);
            hm(j) = Mm(j)/((1-phi_f)*rho_w*Lt*Lz);
            
            num_it = num_it + 1;
        end
        
        for k=1:4
            mPk_p(k) = mPk(k);
            sigma_k_p(k) = sigma_k(k);
        end
        mPm_p = mPm;
        sigma_m_p = sigma_m;
        
        [pk1_p(:,j), pk2_p(:,j), pk3_p(:,j), pm_p(:,j), r_p(:,j)] = time_evolv(pk1_a(:,j), pk2_a(:,j), pk3_a(:,j), pm_a(:,j), r_p(:,j), r(c,j));
        

        
        r_act(j) = r_act(j) + k_act*(r(c,j)-r_act(j))*dt;
    end

    if mod(c-1,20) == 0
        
        for j=1:N_gen
        
        scatter(t,r(:,j),'.'); hold on; 
        pause(0.01);
        
        end
%                     figure(3);
%                     scatter(t(c),T_act/hm,'.r'); hold on;
%                     pause(0.1)
%                     figure(4);
%                     scatter(t(c),P*r(c)/h,'.r'); hold on;
%                     pause(0.1)
    end
    er = sqrt(sum((r(c,:)-r(c-1,:)).^2)/sum(r(c-1,:).^2))
end














