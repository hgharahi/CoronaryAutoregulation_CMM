for j = 1:length(Pin)
    
%     er = 1000;
%     c = 0;
%     while er > 10^-6
        
        
%         c = c+1;
%         t(c+1) = t(c) + dt;
        
                
        PF = 6*myogenic_control(Pin(j) - subendo.Pim , subendo.Pmid(i) - subendo.Pim, ID);
        if a==0
            PF = 0;
        end
%           x = mechanical_params(subendo, i, PF, ID);
        x = [113*0.1*10/3 333*0.1*10/3 4.64 78*0.1*10/3 3.41 1.12 1.12 1.07 0.95 1.063 0000*10/3 0.7 1.7];
%         Me = subendo.Me(i);
%         Mm = subendo.Mt(i)*(subendo.SMCtoCOL(i)/(subendo.SMCtoCOL(i)+1));
%         Mc = subendo.Mt(i)*(1/(subendo.SMCtoCOL(i)+1));
%         Mk = [0.1, 0.1, 0.4, 0.4]*Mc;
        
        P = Pin(j) - subendo.Pim;
        
        hm = 1/subendo.rho_w;% subendo.Thickness(i)*Mm/(Mm+Mc+Me)*3/10;
        
        flg = 0;
        R1(i) = NR_iterate_fit(Me,Mk,Mm,subendo.Radius(i),hm,P,x, r_act,r_h, flg, subendo, i, ID);
        LL(i) = R1(i)/R(i);
        
        
        
        er = sqrt(sum(LL(i)-1).^2);
        R = R1;
        
%     end
    R(i)/subendo.Radius(i)
    scatter(abs(R(i)/R0),(Pin(j)-subendo.Pim)/133.32,'.r');hold on;
end