function [Ks0, Kp0, Km0, phi_s, phi_p, phi_m, a] =  stimuli_params(tree, Ks, Kp, Km)

for i=1:tree.N_gen
    S00(i) = (-log( 1/tree.Act_lvl(i) - 1 ));
    
    [phi_s(i), phi_p(i), phi_m(i)] = signal_weights(tree.R(i));
    
    Kp0(i) = S00(i) / (-phi_m(i)*Km/Kp - phi_p(i) + phi_s(i)*Ks/Kp);
    Ks0(i) = Ks/Kp*Kp0(i);
    Km0(i) = Km/Kp*Kp0(i);
    a(i) = 1;
end