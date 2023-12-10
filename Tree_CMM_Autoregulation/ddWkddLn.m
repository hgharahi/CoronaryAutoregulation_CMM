function  y=ddWkddLn(Ln)

global kc

exp_Q = exp(kc(3)*(Ln^2-1.0)^2);
y = kc(2)* (3*Ln^2-1.0+4*kc(3)*(Ln*(Ln^2-1.0))^2)*exp_Q;