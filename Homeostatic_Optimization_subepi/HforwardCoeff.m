% updated 10-17-17
function [Hfn] = HforwardCoeff(NumModes,omegan,L,P0,cn,Gamma)

    %get coefficient for pressure in forward reflective wave
    Hfn = zeros(NumModes,1);
    for k=2:NumModes;
        Hfn(k)= P0(k)/(1.0 + Gamma(k)*exp(-1i*omegan(k)*2.0*L/cn(k)));
    end
       
end