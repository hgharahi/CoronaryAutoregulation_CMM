function F = myogenic_control_subepi(P, Ptm_h, ID)

Fmax = 1;
% C's are calculated based on the napkin calculations! Hamid has a picture
% coeff: a correction coefficient for the corresopnding homeostatic
% pressure from the optimization!

if P<0
    F = 0;
else
    if ID == 'B'
        Kt = 50*133.32;
        s = 1.45;
        C = 6.3/2;
        coeff = 1.85*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
    elseif ID == 'C'
        Kt = 50*133.32;
        s = 1.45;
        C = 8.58/2;
        coeff = 1.8*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
    elseif ID == 'D'
        Kt = 55*133.32;
        s = 1.50;
        C = 7.1716/2;
        coeff = 1.40*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
    else
        Kt = 80*133.32;
        s = 0.8;
        C = 3.5675*1.0;
        coeff = 1.35*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
        
    end
    
    F = 2*C*(Fmax*(P^s/((Kt)^s+P^s)))/coeff;
end

