function [KC] = mechanical_properties_LCA(p,R)

 global Ghe1 Ghe2 Ghc Ghm S rho_w lmax lmin fiber_ang
 
 rho_w = 1060;
 lmax = 1.76;
 lmin = 0.35;
% Passive wall
% elastin
Ghe1 = 1.69;% nondim.
Ghe2 = 1.37;% nondim.
% Sigma_e = p/(HtoR(R));
% Sigma_c = p/(HtoR(R)*0.8);%p/(0.07); %120000;
% Sigma_m = p/(HtoR(R));%15*133.32/(0.05/20*k+0.1);
S = 89554;
% c1 =  Sigma_e / ((Ghe1^2-(Ghe1^2*Ghe2^2)^(-1))*(0.3*rho_w));
c1 = 116.76;
c3 = 5.30; % for collagen fibers
c5 = 1.72;  % for passive smooth muscle
Sigma_e = c1*((Ghe1^2-(Ghe1^2*Ghe2^2)^(-1))*(0.3*rho_w));
% fiber_ang = 45*pi/180; % angle of helical collagen fibers
Ghc = 1.18;% nondim. %HG changed from 1.034. June 22nd
str_ch = 0.0;
fiber_ang = 1.0690e+00;
ang = [0.0 90.0 1.0690e+00 -1.0690e+00];
c_frac0 = [.1 .1 .4 .4];
for i=1:4
    str_ch = str_ch + c_frac0(i)*Ghc^2*(Ghc^2-1)*...
        exp(c3*(Ghc^2-1.0)^2)*(sin(ang(i)))^2;
end
% c2 = Sigma_c/ (str_ch*(0.3*rho_w));
c2 = 496.2;
Sigma_c = c2*(str_ch*(0.3*rho_w));
Ghm = 1.36;% nondim.
VSM_act = dWmdx_act(1)*0.3*rho_w;
% c4  = (Sigma_m-VSM_act)/(Ghm^2*(Ghm^2-1.0)*exp(c5*(Ghm^2-1.0)^2)*(0.3*rho_w));
c4 = 68.49;%(Sigma_m-VSM_act)/(Ghm^2*(Ghm^2-1.0)*exp(c5*(Ghm^2-1.0)^2)*(0.3*rho_w));

Sigma_m = c4*(Ghm^2*(Ghm^2-1.0)*exp(c5*(Ghm^2-1.0)^2)*(0.3*rho_w)) + VSM_act;

KC = [c1, c2, c3, c4, c5];

end