while err > 10^-12
    Lt = r/R0;
    % h = Mtotal/(rho_w*Lt*Lz);
    
    % elastin: isotropic matrix
    Lt_e = Ghe1*Lt;
    Lz_e = Ghe2*Lz;
    
    dwdLt = Me * Ghe1 * dWedLn(Lt_e,Lz_e, 1);
    d2wddLt = Me * Ghe1^2 * ddWeddLn(Lt_e,Lz_e, 1);
    
    % collagen: four fiber families
    alpha = [0, pi/2, atan(Lt/Lz*tan(fiber_ang)), -atan(Lt/Lz*tan(fiber_ang))];
    dwdLtc = 0;
    d2wddLtc = 0;
    for k = 1:4
        
        Lk = sqrt(Ghc^2*(Lt^2*sin(alpha(k))^2+Lz^2*cos(alpha(k))^2));
        
        dwdLtc = Mk(k) * ( Ghc^2 * Lt * sin(alpha(k))^2/Lk ) * dWkdLn(Lk) + dwdLtc;
        
        d2wddLtc = Mk(k) * ( (2* Ghc^2 * Lt * sin(alpha(k))^2/Lk )^2 * ddWkddLn(Lk) + ...
            Ghc^4 * ( 2*Lz^2*sin(alpha(k))^2*cos(alpha(k))^2/Lk^3 ) *  dWkdLn(Lk) ) + d2wddLtc;
        
    end
    
    dwdLt = dwdLt + dwdLtc;
    d2wddLt = d2wddLt + d2wddLtc;
    
    % SMC: oriented circumferentially, passive and active response
    % Passive response
    Lm = Ghm*Lt;
    
    dwdLtm = (Mm) * ( Ghm) * dWmdLn(Lm);
    
    d2wddLtm = (Mm) * (2 * Ghm^2 * ddWmddLn(Lm));
    
    dwdLt = dwdLt + dwdLtm;
    d2wddLt = d2wddLt + d2wddLtm;
    
    % Active muscle tone (From Biochemomechanical model)
    L_act = 1; % no time evolution involved (see Eq 15, Getachew's draft)
    S = 2*S_basal*(1+tanh( -kton * (tau/tau_h-1) - log(3/2)));
    f = 1 - (Lmax-L_act)^2/(Lmax-Lmin)^2;
    T_act = t_act * hm;
    
    % Function and derivaite
    F = 1/Lz * dwdLt + T_act - P * r;
    dFdr = 1/R0 * 1/Lz * d2wddLt - P;
    
    r_new = r - F/dFdr;
    err = abs((r_new-r)/r);
    r = r_new;
    
end