for j=1:subendo.N_gen
    Res(j) = 8*subendo.mu(j)*(subendo.Length(j))/(pi*subendo.Radius(j)^4);
end


subendo.Cap_res = (subendo.Pin - subendo.Pout)/subendo.q(end);

Res(subendo.N_gen+1) = subendo.Cap_res;