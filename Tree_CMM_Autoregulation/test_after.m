j = 1;

tau(j) = 4*mu*q_seg(j)/(pi*r(c-1,j)^3);

p_mid(j) = Table(j,3)+ Res(j)/2*q_seg(j);


err = 1000;
num_it = 1;
while err > 10^-6 || num_it < 5
    
    [r(c,j), sigma_k, sigma_m, T_act] = NR_iterate( Me(1,j),pk1_a(:,j),pk2_a(:,j),pk3_a(:,j),pm_a(:,j),r(1,j),p_mid(j),dt,r_p(:,j),r_act(j), j);
    Lt = r(c,j)/R0;
    
    for k=1:4
        mPk(k) = f_mP(sigma_k(k), sigma_h, kg, tau(j), tau_h, kg_sh, m_basal_c(k));
    end
    
    mPm = f_mP(sigma_m, sigma_h, kg, tau(j), tau_h, kg_sh, m_basal_m);
    
    [Mk_new(:,j), Mm_new(:,j), pk1_a(:,j), pk2_a(:,j), pk3_a(:,j), pm_a(:,j)] = ...
        age_distribution(pk1_p(:,j), pk2_p(:,j), pk3_p(:,j), pm_p(:,j), mPk, mPm, num_pa, dt, kq_c);
    
    err = sqrt(sum((Mk_new(:,j)-Mk(:,j)).^2)/sum(Mk(:,j).^2) + (Mm_new(:,j)-Mm(:,j))^2/Mm(:,j)^2);
    Mk(:,j) = Mk_new(:,j);
    Mm(:,j) = Mm_new(:,j);
    Mc(:,j) = sum(Mk(:,j));
    total_M = (Me(j)+Mm(j)+Mc(j));
    h(j) = total_M / ((1-phi_f)*rho_w*Lz*Lt);
    hm(j) = Mm(j)/((1-phi_f)*rho_w*Lt*Lz);
    
    num_it = num_it + 1;
end