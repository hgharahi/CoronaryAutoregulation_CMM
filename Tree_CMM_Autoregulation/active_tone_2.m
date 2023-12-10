function [t_act, dt_actdr] = active_tone_2(r, r_act, r_h)

global S_basal Lmax Lmin kton


S = S_basal;% * sech(kton* (-1 + r_h^3/r^3))^2;

dSdr = 0;%(3 * kton * r_h^3* S_basal* sech(kton* (-1 + r_h^3/r^3) + log(3/2))^2)/(r^4);

L2m_act=r/r_h;

if ( ((Lmax-L2m_act)/(Lmax-Lmin))^2 <= 1 )
    
    %% compute membrane normal stress due to vascular smooth muscle tone - active stress of SMC
%     t_act = S*L1 * ( 1- ((Lmax-L1)/(Lmax-Lmin))^2 );
    if (r~=0)        
        %     dt_actdr = S* ( 1- ((Lmax-L1)/(Lmax-Lmin))^2 ) + ...
        %         2 *  S * (Lmax-L1)*L1/((Lmax-Lmin)^2);
        t_act = S* ( 1- ((Lmax-L2m_act)/(Lmax-Lmin))^2 );
        
        dt_actdr =  dSdr * (1-((Lmax-L2m_act)/(Lmax-Lmin))^2) + 2 *  S * (Lmax-L2m_act)/((Lmax-Lmin)^2);
        
    else
        
        t_act = 0;
        dt_actdr = 0;
    end
else
    t_act = 0;
    dt_actdr = 0;
end