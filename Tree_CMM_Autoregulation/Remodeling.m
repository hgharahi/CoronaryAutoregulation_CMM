function [LL, q_network] = Remodeling(tree, p_in, p_out, af, ap, am)

global k_act

Kp = ap*0.5;
Km = am*0.5;
Ks = af*0.5;

% MF = 1;
t_final = 0.0417;

M = 70;
[tree.q_incr, tree.gol] = metabolic_control(M);
q_target = tree.q(tree.N_gen) * tree.q_incr;

mu = tree.mu;


for i=1:tree.N_gen
    
    Res(i) =8*mu*(tree.Length(i))/(pi*tree.R(i)^4);
    
end

[q_seg, p_term] = SymmetricArterialTree_PinPout(tree.N_gen,tree.Pin,tree.Pout,Res);

p_mid(:) = p_term(:) + (Res(:)/2).*q_seg(:);
tau(:) = 4*mu*q_seg(:)./(pi*tree.R(:).^3);
for i = 1:tree.N_gen
    
    
    
    FF0(i) = shear_control(tau(i), tree.tau_h(i), tree.ID(i), tree.Act_lvl(i));
    
    PF0(i) = 1/tree.Act_lvl(i)*0.5*myogenic_control(p_mid(i) - tree.Pim , tree.Pmid(i) - tree.Pim, tree.ID(i));%% the  multiplier 0.5 here is to get the basal tone!
    
    MF0(i) = 1;
    
    S0(i) = 0;
end


M = 70;
%         [tree.q_incr, tree.gol] = metabolic_control(M);
tree.q_incr = M/70;
q_target = tree.q(end) * tree.q_incr;

j = 1;
S = S0;
er = 1000;
% c = 0;
for j = 1 :length(p_in)
    clf
    tree.t = 0;
    c = 0;
    while tree.t < t_final
        
        c = c+1;
        tree.t(c+1) = tree.t(c) + tree.dt;
        
        %     if tree.t(c) > 0.0417 && tree.t(c) < 0.0417 + tree.dt
        %         M = 70;
        %         [tree.q_incr, tree.gol] = metabolic_control(M);
        %         q_target = tree.q(end) * tree.q_incr;
        %     end
        
        
        for i=1:tree.N_gen
            
            Res(i) =8*mu*(tree.Length(i))/(pi*tree.R(i)^4);
            
        end
        
        [q_seg, p_term] = SymmetricArterialTree_PinPout(tree.N_gen,p_in(j),p_out,Res);
        
        p_mid(:) = p_term(:) + (Res(:)/2).*q_seg(:);
        tau(:) = 4*mu*q_seg(:)./(pi*tree.R(:).^3);
        
        %     M = metabolic_control(i, q_seg(end));
        
        
        for i = 1:tree.N_gen
            
            
            
            %         FF = shear_control(tau(i), tree.tau_h(i), tree.ID(i), tree.Act_lvl(i));
            
            PF = 1/tree.Act_lvl(i)*0.5*myogenic_control(p_mid(i) - tree.Pim , tree.Pmid(i) - tree.Pim, tree.ID(i))/PF0(i);%% the  multiplier 0.5 here is to get the basal tone!
            
            A = 1.5*1/(1+exp(-S(i)));
            x = mechanical_params(tree, i, A*PF, tree.ID(i));
            
            %             Sb(i) = PF*subendo.Sbasal(i);
            
            Me = tree.Me(i);
            Mm = tree.Mt(i)*(tree.SMCtoCOL(i)/(tree.SMCtoCOL(i)+1));
            Mc = tree.Mt(i)*(1/(tree.SMCtoCOL(i)+1));
            Mk = [0.1, 0.1, 0.4, 0.4]*Mc;
            P = p_mid(i) - tree.Pim;
            hm = 1/tree.rho_w;% subendo.Thickness(i)*Mm/(Mm+Mc+Me)*3/10;
            
            R1(i) = NR_iterate_fit( Me, Mk, Mm, tree.Radius(i), hm, P, x, tree.r_act(i), tree.r_h(i));
            LL(i) = R1(i)/tree.R(i);
            
            %             tree.r_act(i) = tree.r_act(i) + k_act*(R1(i)-tree.r_act(i));
            
            %         S(i) = Km*(q_seg(end) - q_target)/q_target*tree.dt + Ks*(tau(i) - tree.tau_h(i))/tree.tau_h(i)*tree.dt + Kp*(p_mid(i) - tree.Pmid(i))/tree.Pmid(i)*tree.dt ;
        end
        
        for i = 1:tree.N_gen
            S(i) = Km*(q_seg(end) - q_target)/q_target + Ks*(tau(i) - tree.tau_h(i))/tree.tau_h(i) -  Kp*(p_mid(i) - tree.Pmid(i))/tree.Pmid(i) ;
        end
        
        tree.R = R1;
%         if mod(c,1) == 0
%             scatter(tree.t(c),q_seg(end)/q_target,'.k');hold on;
%             pause(0.0001);
%         end
        %         R1;
    end
    %     tree.R./tree.Radius
    %     for i = 1:subendo.N_gen
    q_network(j) = q_seg(1);
    % scatter((p_in(j))/133.32,tree.t,'ok');hold on;
end

%     semilogx(R,p_mid/133.32,'*');
LL = tree.R./tree.Radius;
% plot((p_in)/133.32,q_network./tree.q(1),marker,'color',0.1*[1 1 1]);hold on;
%
% t(end)
% c
% tau./tree.tau_h
