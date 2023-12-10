pk1_a = zeros(num_pa,1);
pk2_a = zeros(num_pa,1);
pk3_a = zeros(num_pa,1);

pk1_p = zeros(num_pa,1);
pk2_p = zeros(num_pa,1);
pk3_p = zeros(num_pa,1);

pm_a = zeros(num_pa,1);
pm_p = zeros(num_pa,1);
r_p = zeros(num_pa,1);

a_mean = 0;

pk1_p(1) = 1.0;
pk2_p(1) = 1.0;
pk3_p(1) = 1.0;
pm_p(1) = 1.0;
r_p(1) = R0;

for i=2:num_pa

    pk1_p(i) = pk1_p(i-1)*exp(-0.5*dt*(mu_F(i*dt,kq_c)+mu_F((i-1)*dt,kq_c)));
    pk2_p(i) = pk2_p(i-1)*exp(-0.5*dt*(mu_F(i*dt,kq_c)+mu_F((i-1)*dt,kq_c)));
    pk3_p(i) = pk3_p(i-1)*exp(-0.5*dt*(mu_F(i*dt,kq_c)+mu_F((i-1)*dt,kq_c)));
    
    pm_p(i) = pm_p(i-1)*exp(-0.5*dt*(mu_F(i*dt,kq_c)+mu_F((i-1)*dt,kq_c)));
    
    a_mean = a_mean + 0.5*dt*(pk2_p(i-1)+pk2_p(i));
    r_p(i) = R0;
end

pk1_p = pk1_p/a_mean;
pk2_p = pk2_p/a_mean;
pk3_p = pk3_p/a_mean;

pm_p = pm_p/a_mean;



