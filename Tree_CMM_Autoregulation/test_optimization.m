function tau_h = test_optimization(tree)

% intitalize L,HydRes and Table
for j=1:tree.N_gen
    Res(j) = 8*tree.mu(j)*(tree.Length(j))/(pi*tree.Radius(j)^4);
end
Table = SymmetricArterialTree2(tree.q(1),tree.Pout,tree.N_gen,Res);

q_seg(:) = Table(:,2);
Res(:)= Table(:,4);

p_mid(:) = Table(:,3) + (Res(:)/2).*q_seg(:);
tau_h(:) = 4*tree.mu(:).*q_seg(:)./(pi*tree.Radius(:).^3);

%% Check the baseline vs. the current code! if LL=1, then the parameter conversion is correct!
for i = 1:tree.N_gen

        P = p_mid(i) - tree.Pim;
        
    Smax = Pressure_Dependent_Tension(P, tree.A.Smax, tree.A.beta, tree.A.p0);
    
    x = [tree.c1(i), tree.c2(i), tree.c3(i), tree.c4(i), tree.c5(i), ...
        tree.Ge1(i), tree.Ge2(i), tree.Gm(i), tree.Gc(i), tree.angle(i), ...
        tree.Act_lvl(i)*Smax, tree.l_min(i), tree.l_max(i)];
    
    Me = tree.Me(i);
    Mm = tree.Mt(i)*(tree.SMCtoCOL(i)/(tree.SMCtoCOL(i)+1));
    Mc = tree.Mt(i)*(1/(tree.SMCtoCOL(i)+1));
    Mk = [0.1, 0.1, 0.4, 0.4]*Mc;

    hm = 1/tree.rho_w;% subendo.Thickness(i)*Mm/(Mm+Mc+Me)*3/10;
    flg = 0;
    R(i) = NR_iterate_fit(Me,Mk,Mm,tree.Radius(i),hm,P,x, 1.1*tree.Radius(i), tree.Radius(i), flg, tree, i, tree.ID);
    LL(i) = R(i)/tree.Radius(i);
    
end
LL
if norm(LL-1) < 1e-5
    disp('Parameters are correct!');
else
    disp('Parameter conversion is not correct!!!!!');
end