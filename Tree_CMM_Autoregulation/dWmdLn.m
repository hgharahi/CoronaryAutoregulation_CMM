function  y=dWmdLn(Ln)

global kc

y= kc(4)*Ln*(Ln^2-1.0)*exp(kc(5)*(Ln^2-1.0)^2);