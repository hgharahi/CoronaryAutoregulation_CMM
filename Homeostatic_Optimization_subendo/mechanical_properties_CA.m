function [KC, Act_lvl, Smax] = mechanical_properties_CA(p,R)

global Ghe1 Ghe2 Ghc Ghm S rho_w lmax lmin fiber_ang A LA IA SA

if R < 2.5e-5
    x = SA;
elseif R >= 2.5e-5 && R < 5.0e-5
    x = IA;
elseif R >= 5.0e-5 && R < 9.5e-5
    x = LA;
elseif R >= 9.5e-5
    x = A;
else
    disp('error in parameter assignment');
end

p = p/133.32;
rho_w = 1060;

corr1 = R/x.Rh;
corr2 = p/(x.Ph/133.32);

c1 = x.c1;
c2 = x.c2;
c3 = x.c3;
c4 = x.c4 ;
c5 = x.c5;

Ghe1 = x.Ghe1;
Ghe2 = x.Ghe2;
Ghm = x.Ghm;
Ghc = x.Ghc;

Act_lvl = x.Act_lvl0;
if R < 2.5e-5
    Act_lvl = 1.075*x.Act_lvl0*(corr2)^(0.0);
elseif R >= 2.5e-5 && R < 5.0e-5
    Act_lvl = 1.04*x.Act_lvl0*(corr2)^(1.0);
elseif R >= 5.0e-5 && R < 9.5e-5
    Act_lvl = 1.02*x.Act_lvl0*(corr2)^(2.5);
elseif R >= 9.5e-5
    Act_lvl = 1.03*x.Act_lvl0 * (corr2)^(3.0);
end


fiber_ang = x.fiber_ang;

% S = 0.5*x.S;
lmax = x.lmax;
lmin = x.lmin;

%% Elastin
Sigma_e = c1*((Ghe1^2-(Ghe1^2*Ghe2^2)^(-1))*(0.3*rho_w));

Smax = Pressure_Dependent_Tension(p*133.32, x.Smax, x.beta, x.p0);
S = Act_lvl * Smax;


%% Collagen
str_ch = 0.0;
ang = [0.0 90.0 fiber_ang -fiber_ang];
c_frac0 = [.1 .1 .4 .4];
for i=1:4
    str_ch = str_ch + c_frac0(i)*Ghc^2*(Ghc^2-1)*...
        exp(c3*(Ghc^2-1.0)^2)*(sin(ang(i)))^2;
end
Sigma_c = c2*(str_ch*(0.3*rho_w));

%% SMCs
VSM_act = dWmdx_act(1)*rho_w;
Sigma_m = c4*(Ghm^2*(Ghm^2-1.0)*exp(c5*(Ghm^2-1.0)^2)*(0.3*rho_w)) + 0.3*VSM_act;


KC = [c1, c2, c3, c4, c5];

end