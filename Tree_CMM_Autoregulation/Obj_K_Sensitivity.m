function Y = Obj_K_Sensitivity(K, subendo,  subepi)


global p_in  p_out  M a f_in

Km1 = K(1);
Kp1 = 0;
Ks1 = K(2);

Km2 = K(3);
Kp2 = 0;
Ks2 = K(4);

[q_network1, A1] = tree_Autoreg(subendo, p_in,  p_out,  M, a, Km1, Kp1, Ks1);
[q_network2, A2] = tree_Autoreg(subepi, p_in,  p_out,  M, a, Km2, Kp2, Ks2);

Y = sqrt(sum((q_network1./subendo.q(end)-f_in').^2)) + ...
    sqrt(sum((q_network2./subepi.q(end)-f_in').^2)) + ...                                                    
                                                        0.01*(Ks1-1)^2 + ....
                                                        10*(max(Km1, 0)==0) + ...
                                                        0.01*(Ks2-1)^2 + ...
                                                        10*(max(Km2, 0)==0) + ...                                                        
                                                        10.0*(A1(end,1)) + ...
                                                        10.0*(A2(end,1));



