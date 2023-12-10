function P = Find_Pressure(Me,Mk,Mm,r,hm,r_h,x)

global kc Ghe1 Ghe2 Ghc Ghm fiber_ang S_basal Lmax Lmin

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


Lz = 1.0;

r0 = r;
err = 1000;
alpha = [0, pi/2, atan(tan(fiber_ang)), -atan(tan(fiber_ang))];


%% Constructing the tension

Lt = r0/r_h;

% elastin: isotropic matrix
Lz_e2 = Ghe2^2*Lz^2;
Lt_e2 = Ghe1^2*Lt^2;

dwdLt2 = Me * Ghe1^2 * dWedLn2(Lt_e2,Lz_e2, 2,kc);

% collagen: four fiber families produced at differnet times
dwdLtc = 0;
dwdLzc = 0;

for k = 1:4
    
    
    Lk2 = Lt^2*sin(alpha(k))^2+Lz^2*cos(alpha(k))^2;
    Lk2_n = Ghc^2*Lk2;
    
    dwdLtc = Mk(k) * Ghc*Ghc*sin(alpha(k))^2 * dWkdLn2(Lk2_n,kc) + dwdLtc;
    dwdLzc = Mk(k) * Ghc*Ghc*cos(alpha(k))^2 * dWkdLn2(Lk2_n,kc) + dwdLzc;
    
    
end

% Passive SMC

Lm2 = Lt^2;
Lm2_n = Ghm^2*Lm2;

dwdLtm = 0;

dwdLtm = Mm * Ghm*Ghm * dWmdLn2(Lm2_n,kc) + dwdLtm;



dwdLt2 = dwdLt2 + dwdLtc + dwdLtm;

% Active muscle tone (Constant)
[T_act, ~] = active_tone_2(r0, 1, r_h);

% Function and derivaite

P  = (Mm*T_act + 2.0 * Lt * dwdLt2/Lz) / r0 ;




%% Derivatives of strain energy function

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