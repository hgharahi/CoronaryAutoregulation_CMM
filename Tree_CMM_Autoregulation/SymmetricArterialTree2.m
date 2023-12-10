% modified by Vasilina on April 19, 2018
% -Table(k,4) is substituted by Res(k)
function Table = SymmetricArterialTree2(q,p_term,N_gen,Res)

%% Bifurcations
% p(k) = p(k+1) + dp(k+1); dp(k)=q/(2^(k-1))*Res(k)

Table = zeros(N_gen,4); %VF

dp = 0;
p=zeros(N_gen,1);
for k = N_gen:-1:1
    p(k) = p_term + dp; %terminal pressure at k gen, p_term at N_gen
    dp = +q/(2^(k-1))*Res(k) + dp;
end

for k = 1:N_gen
    Table(k,:) = [k, q/(2^(k-1)), p(k) , Res(k)];
end
 