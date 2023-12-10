global N_gen 

%% initialization
for k=1:N_gen
    lhat0(k) = (1/(2^(k-1)))^(1/2.75);
end

options = optimset('MaxFunEvals',16000,'MaxIter',16000,'TolFun',1e-10,'TolX',1e-8);

[lhat] = fminsearch(@objective, lhat0, options)

[Cost,q,p_term,Mtotal,Table,Radius,SMCtoCOL,Me,Mt,p_mid,ksi,Length] = evaluate(lhat);

