function Smax = Pressure_Dependent_Tension(P, Sm, beta, p0)

%     Sm = 2.45*1e6;
    p0 = p0*133.32;
%     beta = 2.1;
    
    if P > 0
    
        Smax = Sm * P^beta/(P^beta + p0^beta);
        
    else
        
        Smax = 0;
    
    end

end
