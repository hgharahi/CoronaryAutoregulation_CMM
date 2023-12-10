function [sigma_k, sigma_m] = stress_compute(dwdLtc, dwdLzc, dwdLtm, T_act, Lt, Lz,j, Mm, Mc)

global phi_f rho_w fiber_ang

hm = Mm/((1-phi_f)*rho_w*Lt*Lz);
hc = Mc/((1-phi_f)*rho_w*Lt*Lz);

Tt_c(1) = 2*Lt*dwdLtc/(Lz*hc);

Tz_c(1) = 2*Lz*dwdLzc/(Lt*hc);

Tt_m(1) = (2*Lt*dwdLtm/Lz + T_act)/hm;

alpha = [0, pi/2, atan(Lt/Lz*tan(fiber_ang)), -atan(Lt/Lz*tan(fiber_ang))];

sigma_k = [0,0,0,0];
for k=1:4
    sigma_k(k) = sqrt(Tz_c(1)^2*cos(alpha(k))^2 + Tt_c(1)^2*sin(alpha(k))^2);
end

sigma_m = Tt_m(1);