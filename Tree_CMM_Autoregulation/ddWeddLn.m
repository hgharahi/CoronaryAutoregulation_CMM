function  y=ddWeddLn(Ln_t, Ln_z, f)

global kc

if f==1
    y = kc(1) * (1.0+3.0/(Ln_t^4*Ln_z^2));
elseif f==2
    y = kc(1) *(1.0+3.0/(Ln_z^4*Ln_t^2));
elseif  f==3
    y = kc(1) * 2.0 /(Ln_t*Ln_z)^3;
else
    exit('Wrong parameter for ddWddLn');
end