function Autoregulation(tree, p_in, p_out, af, ap, am, a, marker)

global k_act kton

mu = tree.mu;

for j = 1:length(p_in)
    
    er = 1000;
    c = 0;
    while er > 10^-6
        
        
        c = c+1;
        tree.t(c+1) = tree.t(c) + tree.dt;
        
        for i=1:tree.N_gen
            
            Res(i) =8*mu*(tree.Length(i))/(pi*tree.R(i)^4);
            
        end
        
        [q_seg, p_term] = SymmetricArterialTree_PinPout(tree.N_gen,p_in(j),p_out,Res);
        
        p_mid(:) = p_term(:) + (Res(:)/2).*q_seg(:);
        tau(:) = 4*mu*q_seg(:)./(pi*tree.R(:).^3);
        
        %     M = metabolic_control(i, q_seg(end));
        
        
        for i = 1:tree.N_gen
            
            
            if af == 1
                FF = shear_control(tau(i), tree.tau_h(i), tree.ID(i), tree.Act_lvl(i));
            else
                FF = tree.Act_lvl(i);
            end
            
            if ap ~= 0
                PF = ap*1/tree.Act_lvl(i)*0.5*myogenic_control(p_mid(i) - tree.Pim , tree.Pmid(i) - tree.Pim, tree.ID(i));%% the  multiplier 0.5 here is to get the basal tone!
            else
                PF = 1;
            end

           if am == 1
               
               if i == tree.N_gen-1 || i == tree.N_gen
                   
                   MF = 1/tree.gol;
               else
                   MF = 1;
               end
           else
               MF = 1;
           end
           
            x = mechanical_params(tree, i, a*MF*FF*PF, tree.ID(i));
            
            %             Sb(i) = PF*subendo.Sbasal(i);
            
            Me = tree.Me(i);
            Mm = tree.Mt(i)*(tree.SMCtoCOL(i)/(tree.SMCtoCOL(i)+1));
            Mc = tree.Mt(i)*(1/(tree.SMCtoCOL(i)+1));
            Mk = [0.1, 0.1, 0.4, 0.4]*Mc;
            P = p_mid(i) - tree.Pim;
            hm = 1/tree.rho_w;% subendo.Thickness(i)*Mm/(Mm+Mc+Me)*3/10;
            
            R1(i) = NR_iterate_fit( Me, Mk, Mm, tree.Radius(i), hm, P, x, tree.r_act(i), tree.r_h(i));
            LL(i) = R1(i)/tree.R(i);
            
            tree.r_act(i) = tree.r_act(i) + k_act*(R1(i)-tree.r_act(i));
            
        end
        
        er = sqrt(sum(LL-1).^2);
        tree.R = R1;
%         R1;
    end
%     tree.R./tree.Radius
    %     for i = 1:subendo.N_gen
    q_network(j) = q_seg(1);
%     scatter((p_in(j))/133.32,q_seg(1)/tree.q(1),'ok','MarkerFaceColor',i/tree.N_gen*[0 0 1]);hold on;
    %     end
end
%     semilogx(R,p_mid/133.32,'*');
LL = tree.R./tree.Radius;
plot((p_in)/133.32,q_network./tree.q(1),marker,'LineWidth',2);hold on;
%
% t(end)
% c
% tau./tree.tau_h
