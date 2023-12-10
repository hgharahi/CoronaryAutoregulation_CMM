
for c=1:4
    
    j
    idx = tree_idx(c);
    
    phi_e = y(2*c-1);
    
    phi_c = y(2*c);
    
    phi_m(c) = 1 - phi_e - phi_c;
    
    M_total = (1 - phi_f) * rho_w * H(c);
    
    Me = phi_e * M_total;
    
    Mc = phi_c * M_total;
    Mk = [0.1, 0.1, 0.4, 0.4]*Mc;
    
    Mm = phi_m(c) * M_total;
    
    hm = 1/rho_w;% subendo.Thickness(i)*Mm/(Mm+Mc+Me)*3/10;
    
    
    x0 = [tree.c1(idx), tree.c2(idx), tree.c3(idx), tree.c4(idx), tree.c5(idx), ...
        tree.Ge1(idx), tree.Ge2(idx), tree.Gm(idx), tree.Gc(idx), tree.angle(idx), ...
        tree.Smax(idx), tree.l_min(idx), tree.l_max(idx)]
    
    pressure_tests = pressure_data;
    
    
    
    for i = 1:length(pressure_tests)
        
        PF = Pressure_Dependent_Tension( pressure_tests(i), y(end-1), y(end));
        
        x0(11) = PF;
        
        R_active(i) = NR_iterate_fit(Me , Mk, Mm, R_h(c), hm, pressure_tests(i), x0 );
        
        x0(11) = 0;
        
        R_passive(i) = NR_iterate_fit(Me , Mk, Mm, R_h(c), hm, pressure_tests(i), x0 );
        
    end
    
    PF = Pressure_Dependent_Tension( P_h(c) - Pim(1) , y(end-1), y(end));
    
    x0(11) = PF;
    
    A(c) = Find_Activation(Me , Mk, Mm, R_h(c), hm, P_h(c) - Pim(1), x0 );
    
    PF = Pressure_Dependent_Tension( 100*133.32, y(end-1), y(end));
    
    x0(11) = 0;
    
    R100 = NR_iterate_fit(Me , Mk, Mm, R_h(c), hm, pressure_tests(i), x0 );
    
    figure; hold on;
    plot( pressure_tests , R_active/R100 );
    plot( pressure_tests , R_passive/R100 );
    
    plot( pressure_data , diameter_data_active(:,c), 'd');
    plot( pressure_data , diameter_data_passive(:,c), 'o');
    
    
    scatter(P_h(c) - Pim(1), R_h(c)/R100, '*');
end


figure;
plot(phi_m);