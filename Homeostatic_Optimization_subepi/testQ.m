% test half split
% qInpTime = zeros(Nt,Ngen);
% qTermTime = zeros(Nt,Ngen);
% qSteady(k)-qSteady(k+1)*2==0

for k=1:N_gen-1
    test1=qTermTime(:,k)-qInpTime(:,k+1)*2
end