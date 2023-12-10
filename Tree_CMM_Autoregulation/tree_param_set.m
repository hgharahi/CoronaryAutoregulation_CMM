function tree = tree_param_set(tree)

tree.R = tree.Radius;
tree.r_act = tree.R;
tree.r_h = tree.R;
tree.dt = 0.5*6.9444e-04;
tree.t = 0;
for i=1:tree.N_gen
    if tree.R(i) < 2.5e-5
        tree.ID(i) = 'E';
    elseif tree.R(i) >= 2.5e-5 && tree.R(i) < 5.0e-5
        tree.ID(i) = 'D';
    elseif tree.R(i) >= 5.0e-5 && tree.R(i) < 9.5e-5
        tree.ID(i) = 'C';
    elseif tree.R(i) >= 9.5e-5
        tree.ID(i) = 'B';
    else
        disp('error in parameter assignment');
    end
end
tree.Cap_res = (tree.Pout - 20*133.32)/tree.q(end);

a = 1;

tree.tau_h = test_optimization(tree);