function [A, x0] = eval_single(phi, c,  tree, tree_idx, pressure_data, diameter_data_active, diameter_data_passive, R_h, H, P_h, Pim, beta, p0)


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

pressure_tests = 0:2*133.32:110*133.32;

for i = 1:length(pressure_tests)
    
    PF = Pressure_Dependent_Tension( pressure_tests(i), beta, p0);
    
    x0(11) = PF;
    
    R_active(i) = NR_iterate_fit(Me , Mk, Mm, R_h, hm, pressure_tests(i), x0 );
    
    x0(11) = 0;
    
    R_passive(i) = NR_iterate_fit(Me , Mk, Mm, R_h, hm, pressure_tests(i), x0 );
    
    
end

PF = Pressure_Dependent_Tension( 100*133.32 , beta, p0);

x0(11) = 0;

R100 = NR_iterate_fit(Me , Mk, Mm, R_h, hm, pressure_tests(i), x0 );


PF = Pressure_Dependent_Tension( P_h - Pim , beta, p0);

x0(11) = PF;

A = Find_Activation(Me , Mk, Mm, R_h, hm, P_h - Pim, x0 );

x0(11) = A;

names = {'A','LA','IA','SA'};


h = figure; hold on;
plot( pressure_tests/133.32 , R_active/R100 , '--r','LineWidth',3);
plot( pressure_tests/133.32 , R_passive/R100, '--r','LineWidth',3);

sz = 50;
scatter( pressure_data/133.32 , diameter_data_active(:,c), sz, 'ok', 'filled');
scatter( pressure_data/133.32 , diameter_data_passive(:,c), sz, 'ok');

xlabel('p_{tm} (mmHg)','FontSize',16);
ylabel('D/D_{100}','FontSize',16);
box on;

print(h, '-dpdf', [names{c},'_',tree.name,'.pdf']);

