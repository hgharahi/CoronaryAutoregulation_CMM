function  y=dWkdLn(Ln)

global kc

y=kc(2)*Ln*(Ln^2-1.0)*exp(kc(3)*(Ln^2-1.0)^2);
