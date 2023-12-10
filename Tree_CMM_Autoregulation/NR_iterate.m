function [r0,sigma_k, sigma_m, T_act] = NR_iterate(Me,Mm,Mc,pk1_a,pk2_a,pk3_a,pm_a,r,P,dt,r_p,r_act,j)

global kc Ghe1 Ghe2 Ghc Ghm fiber_ang k_act
 
r_h = r;
Lz = 1.0;
% P = 0;
num_pa = length(pk1_a);
%%
r0 = r;
err = 1000;
max_it = 1000;
num_it = 1;
%% NR iterations
while err > 10^-7 && num_it < max_it
    
    Lt = r0/r_h;
    % h = Mtotal/(rho_w*Lt*Lz);
    %     alpha = [0, pi/2, pi/4, -pi/4];
    
    % elastin: isotropic matrix
    Lz_e2 = Ghe1^2*Lz^2;
    Lt_e2 = Ghe2^2*Lt^2;
    
    dwdLt2 = Me * Ghe2^2 * dWedLn2(Lt_e2,Lz_e2, 2,kc);
    ddwdLt2 = Me * Ghe2^4 * ddWedLn2(Lt_e2,Lz_e2, 2,kc);
    
    % collagen: four fiber families produced at differnet times
    alpha = [0, pi/2, atan(Lt/Lz*tan(fiber_ang)), -atan(Lt/Lz*tan(fiber_ang))];
    dwdLtc = 0;
    dwdLzc = 0;
    ddwd2Ltc = 0;
    dwdLtm = 0;
    ddwd2Ltm = 0;
    for a = 1:num_pa
        
        if (a==1)
            L2_a = r0/r_h;
            wt = 0.5*dt;
        else
            L2_a =  r_p(a)/r_h;
            if a==num_pa
                wt = 0.5*dt;
            else
                wt = dt;
            end
        end
        
        
        
        for k = 1:4
            
            if(k==1)
                pk_a = pk1_a(a);
            elseif(k==2)
                pk_a = pk2_a(a);
            else
                pk_a = pk3_a(a);
            end
            
            Lk2_a =  L2_a^2*sin(alpha(k))^2+Lz^2*cos(alpha(k))^2;
            Lk2 = Lt^2*sin(alpha(k))^2+Lz^2*cos(alpha(k))^2;
            Lk2_n = Ghc^2*Lk2/Lk2_a;
            
            dwdLtc = pk_a * Ghc*Ghc/Lk2_a*sin(alpha(k))^2 * dWkdLn2(Lk2_n,kc) * wt + dwdLtc;
            dwdLzc = pk_a * Ghc*Ghc/Lk2_a*cos(alpha(k))^2 * dWkdLn2(Lk2_n,kc) * wt + dwdLzc;
            
            ddwd2Ltc = pk_a * Ghc*Ghc*Ghc*Ghc/(Lk2_a*Lk2_a)*sin(alpha(k))^4 * ddWkdLn2(Lk2_n,kc) * wt + ddwd2Ltc;
            
        end
        Lm2_a =  L2_a^2;
        Lm2 = Lt^2;
        Lm2_n = Ghm^2*Lm2/Lm2_a;
        
        dwdLtm = pm_a(a) * Ghm*Ghm/Lm2_a * dWmdLn2(Lm2_n,kc) * wt + dwdLtm;
        
        ddwd2Ltm = pm_a(a) * Ghm*Ghm*Ghm*Ghm/(Lm2_a*Lm2_a) * ddWmdLn2(Lm2_n,kc) * wt + ddwd2Ltm;
        
    end
    
    dwdLt2 = dwdLt2 + dwdLtc + dwdLtm;
    ddwdLt2 = ddwdLt2 + ddwd2Ltc + ddwd2Ltm;
    %     Lt_e2
    %     Lk2_n
    %     Lm2_n
    % Active muscle tone (Constant)
    [T_act, dT_actdr] = active_tone(r0, r_h, r_act, Mm, j);
    %     T_act= 0;
    %     dT_actdr = 0;
    % Function and derivaite
    F = 2.0 * Lt * dwdLt2/Lz + T_act - P * r0;
    dFdr = (1.0/r_h) * ( 2*dwdLt2/Lz + 4*Lt*Lt*ddwdLt2/Lz) +dT_actdr - P;
    
    r_new = r0 - F/dFdr;
    err = abs((r_new-r0)/r0);
    
    r_act = r_act + k_act * (r_new - r0);
    r0 = r_new;
 
    num_it = num_it + 1;
end

[sigma_k, sigma_m] = stress_compute(dwdLtc, dwdLzc, dwdLtm, T_act, Lt, Lz,j, Mm, Mc);

%% Derivatives of strain energy function
    function  y=ddWedLn2(Ln_t2, Ln_z2, f, kc)
        
        if f==1 % axial
            y = kc(1) * (1.0/(Ln_t2*Ln_z2^3));
        elseif f==2 % circumferential
            y = kc(1) * (1.0/(Ln_t2^3*Ln_z2));
        elseif  f==3
            y = kc(1) * 2.0 /(Ln_t2*Ln_z2)^3;
        else
            exit('Wrong parameter for ddWddLn');
        end
        
    end

    function  y=ddWkdLn2(Lk2_n,kc)
        
        exp_Q = exp(kc(3)*(Lk2_n-1.0)^2);
        y = kc(2)/2*(1+ 2*kc(3)*(Lk2_n-1)^2)*exp_Q;
        
    end

    function  y=ddWmdLn2(Lm2_n,kc)
        
        exp_Q = exp(kc(5)*(Lm2_n-1.0)^2);
        y = kc(4)/2*(1+ 2*kc(5)*(Lm2_n-1)^2)*exp_Q;
    end

    function  y=dWedLn2(Ln_t2, Ln_z2, axis,kc)
        
        if axis == 1  %axial
            y = kc(1)/2*(1 - 1.0/(Ln_t2*Ln_z2^2));
        elseif axis ==2  %circumferential
            y = kc(1)/2*(1 - 1.0/(Ln_t2^2*Ln_z2));
        else
            exit('wrong parameter in dWedLn');
        end
    end

    function  y=dWkdLn2(Lk2_n,kc)
        
        y=kc(2)/2*(Lk2_n-1)*exp(kc(3)*(Lk2_n-1)^2);
    end

    function  y=dWmdLn2(Lm2_n,kc)
        
        y= kc(4)/2*(Lm2_n-1)*exp(kc(5)*(Lm2_n-1)^2);
    end


end