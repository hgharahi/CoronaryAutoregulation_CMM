function y=ddWcddx(x, kc)

dQdx=4.0*kc(3)*x*(x^2-1);

y=kc(2)*(3*x^2-1+(x^3-x)*dQdx)*exp(kc(3)*(x^2-1)^2);

