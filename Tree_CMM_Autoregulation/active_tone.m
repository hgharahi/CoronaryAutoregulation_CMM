function [t_act, dt_actdr] = active_tone(r, r_h, r_act, Mm, j)

global S_basal Lmax Lmin kton Mm_h


S = S_basal * (Mm/Mm_h(j)) * (1+tanh(-kton*(r_h^3/r^3-1)-log(3/2)));

dSdr = (3 * kton *Mm* r_h^3* S_basal* sech(kton* (-1 + r_h^3/r^3) + log(3/2))^2)/(r^4 * Mm_h(j));


L2m_act=r/r_act;

if ( ((Lmax-L2m_act)/(Lmax-Lmin))^2 <= 1 )
    
    %% compute membrane normal stress due to vascular smooth muscle tone - active stress of SMC
    t_act = S*L2m_act * ( 1- ((Lmax-L2m_act)/(Lmax-Lmin))^2 );
    if (r~=0)
        
        dt_actdr = dSdr * L2m_act * (1-((Lmax-L2m_act)/(Lmax-Lmin))^2) +...
            t_act/r + ...
            2 * L2m_act * S * (Lmax-L2m_act)*L2m_act/(r*(Lmax-Lmin)^2);
       
    else
        
        t_act = 0;
        dt_actdr = 0;
        
    end
else
    t_act = 0;
    dt_actdr = 0;
end