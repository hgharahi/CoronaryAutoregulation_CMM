function err = objective_all(x, subepi, subendo, idx, pressure_data, diameter_data_active, diameter_data_passive, R_h, H, P_h, Pim)

beta = x(end-1); 
p0 = x(end);

j = 1;
for c=1:4
    phi(1) = x(2*j-1);
    phi(2) = x(2*j);
    
    er (j) = objective_single(phi, c, subepi, idx(c), pressure_data, diameter_data_active, diameter_data_passive, R_h(j), H(j), P_h(j), Pim(1), beta, p0);
    
    j = j+1;
end

for c=1:4
    phi(1) = x(2*j-1);
    phi(2) = x(2*j);
    
    er (j) = objective_single(phi, c, subendo, idx(c), pressure_data, diameter_data_active, diameter_data_passive, R_h(j), H(j), P_h(j), Pim(1), beta, p0);
    
    j = j+1;
end

err_col = sum(abs(x(2:2:16) - 0.5));

err = sum(er) + 0.1*err_col + ...
    0.5*abs(x(4)-0.2) + 0.5*abs(x(12)-0.2) +...
    0.5*abs(x(13)-0.2) + 0.5*abs(x(5)-0.2);