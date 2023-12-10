function  y=ddWmddLn(Ln)

global kc

exp_Q = exp (kc(5)*(Ln^2-1.0)^2);
y=kc(4)*(3*Ln^2-1.0+4*kc(5)*(Ln*(Ln^2-1.0))^2)*exp_Q;

