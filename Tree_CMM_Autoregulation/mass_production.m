Lt = r/R0;
Lz = 1.0;

hm = Mm/((1-phi_f)*rho_w*Lt*Lz);
hc = Mc/((1-phi_f)*rho_w*Lt*Lz);

Tt_c(1) = 2*Lz*dwdLzc/(Lt*hc);

Tz_c(1) = 2*Lt*dwdLtc/(Lz*hc);

Tt_m(1) = (2*Lt*dwdLtc/Lz + T_act)/hm;

alpha = [0, pi/2, atan(Lt/Lz*tan(fiber_ang)), -atan(Lt/Lz*tan(fiber_ang))];

for k=1:4
    sigma_k(k) = sqrt(T1_c(1)^2*cos(alpha(k))^2 + T2_c(1)^2*sin(alpha(k))^2);
end

sigma_m = Tt_m(1);