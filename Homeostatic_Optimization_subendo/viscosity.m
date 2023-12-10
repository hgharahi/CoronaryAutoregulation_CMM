function mu_vivo = viscosity(D)

D = D*1e6;

%% relative viscosities, based on Pries et al. 1994
mu_45 = ( 6*exp(-0.085*D) + 3.2 - 2.44*exp(-0.06*D^0.645)); 

mu_vivo = (1 + ( mu_45 - 1) * (D/(D-1.1))^2) * (D/(D-1.1))^2 ; % 

mu_plasma = 1*1e-3;


%% absolute viscosity
mu_vivo = mu_vivo*mu_plasma;