% modified by Vasilina on August 7, 2018
function [CostFun, Mtotal] = mass_optimiz_tree_2(p,Q,k, lhat)

global rho_w beta gamma R0 alpha_t Ghe1 Ghe2 Ghc Ghm fiber_ang Pim



LamZ = 1;


%%
gamma_v = gamma*viscosity(2*lhat*R0)/0.0035;
kc = mechanical_properties_CA(p-Pim,lhat*R0);

alpha_act = 0.00872;

VSM_act = dWmdx_act(1)*0.3*rho_w;

alpha_total = alpha_t + alpha_act*VSM_act;

beta_e = 1/LamZ*dWedx(1, Ghe1, Ghe2,kc);

[phi_e, phi_c, phi_m] = mass_fracs_2(2*lhat*R0);
SMCtoCOL = phi_m/phi_c;
beta_t = 1/LamZ*(1/(1+SMCtoCOL)*(0.1+ 0.4*sin(fiber_ang)^2 + 0.4*sin(-fiber_ang)^2)*Ghc*dWcdx(Ghc,kc) + SMCtoCOL/(1+SMCtoCOL)*(Ghm*dWmdx(Ghm,kc)+dWmdx_act(1)));
phi_t = 1 - phi_e;

Mtotal = ((p-Pim)*R0*lhat)/(phi_e*beta_e+phi_t*beta_t);

CostFun = 2*pi/(0.3*rho_w)*alpha_total*SMCtoCOL/(SMCtoCOL+1)*phi_t/(phi_e*beta_e+phi_t*beta_t)*((p-Pim)*lhat^2*R0^2) + ...
    2*pi/(0.3*rho_w)*alpha_t*1/(SMCtoCOL+1)*phi_t/(phi_e*beta_e+phi_t*beta_t)*((p-Pim)*lhat^2*R0^2)+...
    beta*R0^2*lhat^2 + gamma_v*Q^2*R0^(-4)*lhat^(-4);


% disp(['     Mass-Radius Optimization at gen. ',num2str(k),': Newton-Raphson iterations # ',num2str(c),', Optimized cost ',num2str(CostFun),]);










