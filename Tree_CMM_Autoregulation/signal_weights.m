function [phi_s, phi_p, phi_m] = signal_weights(R)

phi_s = 1/(1+exp(-0.2*1e5*(2*R-250*1e-6)));

phi_p = 1/(exp(((2*R-100*1e-6)/(200*1e-6))^2));

phi_m = 1/(exp(((2*R-50*1e-6)/(75*1e-6))^2));


