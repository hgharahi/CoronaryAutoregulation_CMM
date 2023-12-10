function [Mk, Mm, pk1_a, pk2_a, pk3_a, pm_a] = age_distribution(pk1_p, pk2_p, pk3_p, pm_p, m_basal_c, m_basal_m, num_pa, dt, kq_c)

pk1_a(1) = m_basal_c(1);
pk2_a(1) = m_basal_c(2);
pk3_a(1) = m_basal_c(3);

pm_a(1) = m_basal_m;

Mm = 0;
Mk(:) = [0, 0, 0, 0];

for i=2:num_pa
    
    pk1_a(i) = pk1_p(i-1)*exp(-0.5*dt*(mu_F(i*dt,kq_c)+mu_F((i-1)*dt,kq_c)));
    pk2_a(i) = pk2_p(i-1)*exp(-0.5*dt*(mu_F(i*dt,kq_c)+mu_F((i-1)*dt,kq_c)));
    pk3_a(i) = pk3_p(i-1)*exp(-0.5*dt*(mu_F(i*dt,kq_c)+mu_F((i-1)*dt,kq_c)));
    
    pm_a(i) = pm_p(i-1)*exp(-0.5*dt*(mu_F(i*dt,kq_c)+mu_F((i-1)*dt,kq_c)));

    Mk(1) = Mk(1) + 0.5*dt*(pk1_a(i)+pk1_a(i-1));
    Mk(2) = Mk(2) + 0.5*dt*(pk2_a(i)+pk2_a(i-1));
    Mk(3) = Mk(3) + 0.5*dt*(pk3_a(i)+pk3_a(i-1));
    
    Mm = Mm + 0.5*dt*(pm_a(i)+pm_a(i-1));
end

Mk(4) = Mk(3);