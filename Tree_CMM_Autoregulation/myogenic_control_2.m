function F = myogenic_control_2(P, Ptm_h, ID, tree)

Fmax = 1;
% C's are calculated based on the napkin calculations! Hamid has a picture
% coeff: a correction coefficient for the corresopnding homeostatic
% pressure from the optimization!
if P<0
    F = 0;
else
    
    Kt = 70*133.32;
    s = 2.1;
    C = 6.3/2;
    coeff = 0.55*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
    %         Kt = 70*133.32;
    %         s = 2.1;
    %         C = 8.58/2;
    %         coeff = .7*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
    
    
    
    
    F = 2*C*(Fmax*(P^s/((Kt)^s+P^s)))/coeff;
    
end
end
