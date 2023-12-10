function [Lt0, Lz0] = UnloadedConfig(phi, tree, tree_idx, R_h, H)


phi_e = phi(1);
phi_c = phi(2);
phi_m = 1 - phi_e - phi_c;

phi_f = 0.7;

M_total = (1 - phi_f) * tree.rho_w * H;

Me = phi_e * M_total;

Mc = phi_c * M_total;
Mk = [0.1, 0.1, 0.4, 0.4]*Mc;

Mm = phi_m * M_total;

hm = 1/tree.rho_w;% subendo.Thickness(i)*Mm/(Mm+Mc+Me)*3/10;


x0 = [tree.c1(tree_idx), tree.c2(tree_idx), tree.c3(tree_idx), tree.c4(tree_idx), tree.c5(tree_idx), ...
    tree.Ge1(tree_idx), tree.Ge2(tree_idx), tree.Gm(tree_idx), tree.Gc(tree_idx), tree.angle(tree_idx), ...
    tree.Smax(tree_idx), tree.l_min(tree_idx), tree.l_max(tree_idx)];


%% Material Properties
kc = x0(1:5);
Ghe1 = x0(6);
Ghe2 = x0(7);
Ghm = x0(8);
Ghc = x0(9);
fiber_ang = tree.angle(tree_idx);
% Calculate unloading stress from in vivo to a traction-free
% biaxial specimen

Tol_error = 10^(-8);

PL=[0.8;0.8]; %initial guess for PL

fiber_angle =[0.0, pi/2, fiber_ang, -fiber_ang]; %angle wrt the longitudinal direction
error = 1.0;
while error > Tol_error
    F = [0.0;0.0];
    K = [0.0 0.0; 0.0 0.0];
    L_e1 = PL(1)*Ghe1;
    L_e2 = PL(2)*Ghe2;
    F(1) = F(1)+Me*L_e1*dWedx(1,L_e1, L_e2,kc);
    F(2) = F(2)+Me*L_e2*dWedx(2,L_e1, L_e2,kc);
    K(1,1) = K(1,1)+Me*Ghe1*...
        (dWedx(1,L_e1, L_e2,kc)+L_e1*ddWeddx(1,L_e1,L_e2,kc));
    K(2,2) = K(2,2)+Me*Ghe2*...
        (dWedx(2,L_e1, L_e2,kc)+L_e2*ddWeddx(2,L_e1,L_e2,kc));
    K(1,2) = K(1,2)+Me*Ghe2*L_e1*ddWeddx(3,L_e1,L_e2,kc);
    K(2,1) = K(2,1)+Me*Ghe1*L_e2*ddWeddx(3,L_e1,L_e2,kc);
    
    
    L_m1 = PL(1)*Ghm;
    F(1) = F(1)+ Mm*L_m1*dWmdx(L_m1,kc);
    K(1,1) = K(1,1)+Mm*Ghm*...
        (dWmdx(L_m1, kc)+L_m1*ddWmddx(L_m1,kc));
    
    for k=1:4
        ang = fiber_angle(k);
        L_k = Ghc*sqrt(PL(1)^2*sin(ang)^2+PL(2)^2*cos(ang)^2);

        dL_kdL1 = Ghc^2*PL(1)*sin(ang)^2/L_k;
        dL_kdL2 = Ghc^2*PL(2)*cos(ang)^2/L_k;
        ddL_kddL1 = Ghc^4*PL(2)*sin(ang)^2*cos(ang)^2/L_k^3;
        ddL_kddL2 = Ghc^4*PL(1)*sin(ang)^2*cos(ang)^2/L_k^3;
        ddL_kdL1dL2 = -Ghc^4*PL(1)*PL(2)*sin(ang)^2*cos(ang)^2/L_k^3;
        
        F(1) = F(1)+Mk(k)*dWcdx(L_k,kc)*dL_kdL1*PL(1);
        F(2) = F(2)+Mk(k)*dWcdx(L_k,kc)*dL_kdL2*PL(2);
        K(1,1) = K(1,1)+ Mk(k)*(ddWcddx(L_k,kc)*dL_kdL1^2*PL(1)+...
            dWcdx(L_k, kc)*(ddL_kddL1*PL(1)+dL_kdL1));
        K(2,2) = K(2,2)+ Mk(k)*(ddWcddx(L_k,kc)*dL_kdL2^2*PL(2)+...
            dWcdx(L_k, kc)*(ddL_kddL2*PL(2)+dL_kdL2));
        K(1,2) = K(1,2)+ Mk(k)*(ddWcddx(L_k,kc)*dL_kdL1*dL_kdL2*PL(1)+...
            dWcdx(L_k, kc)*ddL_kdL1dL2*PL(1));
        K(2,1) = K(2,1)+ Mk(k)*(ddWcddx(L_k,kc)*dL_kdL1*dL_kdL2*PL(2)+...
            dWcdx(L_k, kc)*ddL_kdL1dL2*PL(2));
    end
    d_PL=K\F;
    PL = PL - d_PL;
    error = (d_PL(1)^2+d_PL(2)^2)/(PL(1)^2+PL(2)^2);
end
        
Lt0 = PL(1);
Lz0 = PL(2);