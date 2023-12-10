function mP = f_mP(sigma, sigma_h, kg, tau, tau_h, kg_sh, mP_b)


mP = mP_b*(1.0+ kg*(sigma/sigma_h -1.0) - kg_sh*(tau/tau_h - 1));
if ( mP<0 )
    mP = 0.0;
end