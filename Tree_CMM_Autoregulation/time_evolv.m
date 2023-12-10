function [pk1_p, pk2_p, pk3_p, pm_p, r_p] = time_evolv(pk1_a, pk2_a, pk3_a, pm_a, r_p, r)

num_pa = length(pk1_a);

for i=1:num_pa

    pk1_p(i) = pk1_a(i);
    pk2_p(i) = pk2_a(i);
    pk3_p(i) = pk3_a(i);
    
    pm_p(i) = pm_a(i);
    
    if i == num_pa
        r_p(1) = r;
    else
        r_p( num_pa - i + 1) = r_p( num_pa - i );
    end
end