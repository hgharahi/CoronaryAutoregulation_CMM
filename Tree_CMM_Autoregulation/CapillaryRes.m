subendo.Cap_res = (subendo.Pout - 30*133.32)/subendo.q(end);


for j=1:subendo.N_gen
    Res(j) = 8*subendo.mu(j)*(subendo.Length(j))/(pi*subendo.Radius(j)^4);
end


% subendo.Cap_res = (subendo.Pin - subendo.Pout)/subendo.q(end);

Res(subendo.N_gen+1) = 2*subendo.Cap_res;

[q, p_term] = SymmetricArterialTree_PinPout_CapRes(subendo.N_gen,subendo.Pin, 30*133.32 ,Res);