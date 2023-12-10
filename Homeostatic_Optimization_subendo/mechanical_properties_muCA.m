function [KC] = mechanical_properties_muCA(p,R)

 global Ghe1 Ghe2 Ghc Ghm S rho_w lmax lmin fiber_ang
 
 rho_w = 1060;
 lmax = 1.76;
 lmin = 0.35;
% Passive wall
% elastin
Ghe1 = 1.4199e+00;% nondim.
Ghe2 = 1.4509e+00;% nondim.
S = 89554;
Ghc = 1.1066e+00;% nondim.
Ghm =  1.2183e+00;% nondim.
fiber_ang = 1.1066e+00;

c3 = 5.9782e+00; % for collagen fibers
c5 = 1.8475e+00;  % for passive smooth muscle


% c1 =  Sigma_e / ((Ghe1^2-(Ghe1^2*Ghe2^2)^(-1))*(0.3*rho_w));
c1 = 173.0185962;
Sigma_e = c1*((Ghe1^2-(Ghe1^2*Ghe2^2)^(-1))*(0.3*rho_w));
% fiber_ang = 45*pi/180; % angle of helical collagen fibers



% c2 = Sigma_c/ (str_ch*(0.3*rho_w));
c2 = 1380.521753;
str_ch = 0.0;
ang = [0.0 90.0 fiber_ang -fiber_ang];
c_frac0 = [.1 .1 .4 .4];
for i=1:4
    str_ch = str_ch + c_frac0(i)*Ghc^2*(Ghc^2-1)*...
        exp(c3*(Ghc^2-1.0)^2)*(sin(ang(i)))^2;
end
Sigma_c = c2*(str_ch*(0.3*rho_w));

% c4  = (Sigma_m-VSM_act)/(Ghm^2*(Ghm^2-1.0)*exp(c5*(Ghm^2-1.0)^2)*(0.3*rho_w));
c4 = 163.9467;
VSM_act = dWmdx_act(1)*0.3*rho_w;
Sigma_m = c4*(Ghm^2*(Ghm^2-1.0)*exp(c5*(Ghm^2-1.0)^2)*(0.3*rho_w)) + VSM_act;


KC = [c1, c2, c3, c4, c5];

end