function mu_vivo = viscosity(D)

D = D*1e6;

%% relative viscosities, based on Pries et al. 1994
mu_45 = ( 6*exp(-0.085*D) + 3.2 - 2.44*exp(-0.06*D^0.645)); 

Hd_lvl = 0.45;

C = (0.8+exp(-0.07*D))*(-1+1/(1+1e-11*D^12))+1/(1+1e-11*D^12);

HF = ((1-Hd_lvl)^C-1)/((1-0.45)^C-1); %Hematocrit factor

mu_vivo = (1 + ( mu_45 - 1) *HF* (D/(D-1.1))^2) * (D/(D-1.1))^2 ; % 

mu_plasma = 1*1e-3;


%% absolute viscosity
mu_vivo = mu_vivo*mu_plasma;