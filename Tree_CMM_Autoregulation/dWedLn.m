function  y=dWedLn(Ln_t, Ln_z, axis)

global kc 

if axis == 1  %circumferential
    y=kc(1)*(Ln_t - 1.0/(Ln_t^3*Ln_z^2));
elseif axis ==2  %axial
    y= kc(1)*(Ln_z - 1.0/(Ln_z^3*Ln_t^2));
else
    exit('wrong parameter in dWedLn');
end
