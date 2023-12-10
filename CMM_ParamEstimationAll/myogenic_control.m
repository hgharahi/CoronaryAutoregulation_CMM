function F = myogenic_control(P, Ptm_h, ID, tree)

Fmax = 1;
% C's are calculated based on the napkin calculations! Hamid has a picture
% coeff: a correction coefficient for the corresopnding homeostatic
% pressure from the optimization!
if strcmp(tree.name,'subendo'==1)
    if P<0
        F = 0;
    else
        if ID == 'B'
            Kt = 60*133.32;
            s = 2.1;
            C = 6.3/2;
            coeff = 0.55*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
        elseif ID == 'C'
            Kt = 70*133.32;
            s = 2.1;
            C = 6.3/2;
            coeff = 0.55*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
            %         Kt = 70*133.32;
            %         s = 2.1;
            %         C = 8.58/2;
            %         coeff = .7*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
        elseif ID == 'D'
            Kt = 70*133.32;
            s = 2.1;
            C = 6.3/2;
            coeff = 0.55*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
            %         Kt = 60*133.32;
            %         s = 2.50;
            %         C = 7.1716/2;
            %         coeff = .70*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
        else
            Kt = 50*133.32;
            s = 2.1;
            C = 6.3/2;
            coeff = 0.55*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
            %         Kt = 30*133.32;
            %         s = 2.4;
            %         C = 3.5675*1.0;
            %         coeff = .65*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
            
        end
        
        F = 2*C*(Fmax*(P^s/((Kt)^s+P^s)))/coeff;
    end
else
    if P<0
        F = 0;
    else
        if ID == 'B'
            Kt = 60*133.32;
            s = 2.1;
            C = 6.3/2;
            coeff = 0.55*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
        elseif ID == 'C'
            Kt = 70*133.32;
            s = 2.1;
            C = 6.3/2;
            coeff = 0.55*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
            %         Kt = 70*133.32;
            %         s = 2.1;
            %         C = 8.58/2;
            %         coeff = .7*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
        elseif ID == 'D'
            Kt = 70*133.32;
            s = 2.1;
            C = 6.3/2;
            coeff = 0.55*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
            %         Kt = 60*133.32;
            %         s = 2.50;
            %         C = 7.1716/2;
            %         coeff = .70*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
        else
            Kt = 50*133.32;
            s = 2.1;
            C = 6.3/2;
            coeff = 0.55*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
            %         Kt = 30*133.32;
            %         s = 2.4;
            %         C = 3.5675*1.0;
            %         coeff = .65*(Fmax*C*(1/((Kt/Ptm_h)^s+1)));
            
        end
        
        F = 2*C*(Fmax*(P^s/((Kt)^s+P^s)))/coeff;
    end
end
