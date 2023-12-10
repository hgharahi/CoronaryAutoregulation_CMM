function rp = NR_interpolation(Me,Mk,Mm,r,hm,P, r_act, rh, subendo, i, ID)

PF = myogenic_control(subendo.Pmid(i) - subendo.Pim , subendo.Pmid(i) - subendo.Pim, ID);
% if a==0
%     PF = 0;
% end
x1 = mechanical_params(subendo, i, PF, ID);

flg = 1;
rp = NR_iterate_fit(Me,Mk,Mm,r,hm,P,x1, r_act, rh, flg, subendo, i, ID);

