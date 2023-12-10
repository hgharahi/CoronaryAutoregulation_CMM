function [r0, T_act] = NR_iterate_fit(Me,Mk,Mm,r,hm,P,x)

global kc Ghe1 Ghe2 Ghc Ghm fiber_ang S_basal Lmax Lmin k_act

kc=[0.0 0.0 0.0 0.0 0.0];

kc(1)=x(1);
kc(2) = x(2);
kc(3) = x(3);
kc(4) = x(4);
kc(5) = x(5);
Ghe1 = x(6);
Ghe2 = x(7); 
Ghm = x(8);
Ghc = x(9);
fiber_ang = x(10);

S_basal = x(11)*hm;
Lmax = x(13);
Lmin = x(12);

r_h = r;
rh = r;
Lz = 1.0;

r0 = 1.2*r;
err = 1000;
    alpha = [0, pi/2, atan(tan(fiber_ang)), -atan(tan(fiber_ang))];


%% NR iterations
while err > 10^-6
    
    Lt = r0/r_h;
    
    % h = Mtotal/(rho_w*Lt*Lz);
    %     alpha = [0, pi/2, pi/4, -pi/4];
    
    % elastin: isotropic matrix
    Lz_e2 = Ghe2^2*Lz^2;
    Lt_e2 = Ghe1^2*Lt^2;
    
    dwdLt2 = Me * Ghe1^2 * dWedLn2(Lt_e2,Lz_e2, 2,kc);
    ddwdLt2 = Me * Ghe1^4 * ddWedLn2(Lt_e2,Lz_e2, 2,kc);
    
    % collagen: four fiber families produced at differnet times
    dwdLtc = 0;
    dwdLzc = 0;
    ddwd2Ltc = 0;
    dwdLtm = 0;
    ddwd2Ltm = 0;

        
        
        for k = 1:4
            
            
            Lk2 = Lt^2*sin(alpha(k))^2+Lz^2*cos(alpha(k))^2;
            Lk2_n = Ghc^2*Lk2;
            
            dwdLtc = Mk(k) * Ghc*Ghc*sin(alpha(k))^2 * dWkdLn2(Lk2_n,kc) + dwdLtc;
            dwdLzc = Mk(k) * Ghc*Ghc*cos(alpha(k))^2 * dWkdLn2(Lk2_n,kc) + dwdLzc;
            
            ddwd2Ltc = Mk(k) * Ghc*Ghc*Ghc*Ghc*sin(alpha(k))^4 * ddWkdLn2(Lk2_n,kc) + ddwd2Ltc;
            
        end

        Lm2 = Lt^2;
        Lm2_n = Ghm^2*Lm2;
        
        dwdLtm = Mm * Ghm*Ghm * dWmdLn2(Lm2_n,kc) + dwdLtm;
        
        ddwd2Ltm = Mm * Ghm*Ghm*Ghm*Ghm * ddWmdLn2(Lm2_n,kc) + ddwd2Ltm;

    
    dwdLt2 = dwdLt2 + dwdLtc + dwdLtm;
    ddwdLt2 = ddwdLt2 + ddwd2Ltc + ddwd2Ltm;
    %     Lt_e2
    %     Lk2_n
    %     Lm2_n
    % Active muscle tone (Constant)
    [T_act, dT_actdr] = active_tone_2(r0, 1, rh);
%         T_act= 0;
%         dT_actdr = 0;
    % Function and derivaite
    F = 2.0 * Lt * dwdLt2/Lz + Mm*T_act - P * r0;
%      (1.0/r_h)*(Mm*dT_actdr)
    dFdr = (1.0/r_h) * ( 2*dwdLt2/Lz + 4*Lt*Lt*ddwdLt2/Lz) + (1.0/r_h)*(Mm*dT_actdr) - P;
    
    r_new = r0 - 0.5*F/dFdr;
    err = abs((r_new-r0)/r0);
    
%     r_act = r_act + k_act * (r_new - r0)*0.002;
    r0 = r_new;
end
Lt = r0/r_h;


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