function [LL, q_network, PMID, TAU, Active] = Remodeling2t(tree, p_in, p_out, M, af, ap, am, a)

global k_act

Kp = ap*0.5;
Km = am*0.5;
Ks = af*0.5;

% MF = 1;
t_final = 0.52;

M0 = 70;
[tree.q_incr, tree.gol] = metabolic_control(M0);
q_target = tree.q(tree.N_gen) * tree.q_incr;

mu = tree.mu;


for i=1:tree.N_gen
    
    Res(i) =8*mu(i)*(tree.Length(i))/(pi*tree.R(i)^4);
    
end

[q_seg, p_term] = SymmetricArterialTree_PinPout(tree.N_gen,tree.Pin,tree.Pout,Res);

p_mid(:) = p_term(:) + (Res(:)/2).*q_seg(:);
tau(:) = 4*mu(:).*q_seg(:)./(pi*tree.R(:).^3);
for i = 1:tree.N_gen
    
    
    
    FF0(i) = shear_control(tau(i), tree.tau_h(i), tree.ID(i), tree.Act_lvl(i));
    
    Tone_max0(i) = .5*myogenic_control(p_mid(i) - tree.Pim , tree.Pmid(i) - tree.Pim, tree.ID(i));%% the  multiplier 0.5 here is to get the basal tone!
    
    MF0(i) = 1;
    
    S00(i) = (-log( 1/tree.Act_lvl(i) - 1 ));
    %     S0(i) = 0;
end


% M0 = 210;
%         [tree.q_incr, tree.gol] = metabolic_control(M);
tree.q_incr = M/70;
q_target_f = tree.q(end) * tree.q_incr;

q_target = q_seg(end);

S0 = S00;
eps0 = 0;

% er = 1000;
% c = 0;
for j = 1 :length(p_in)
    
    eps_error = 10;
    S0 = S00;
    tree.R = tree.Radius;
    %         clf
    % finding the CEP pressure associated with the given mean
    % pressure at the inlet. The derivations in the notebook.
%     C = tree.Pin/tree.Pim;
%     alpha = p_in(j)/tree.Pin;
    

    
%     beta = C*alpha - 50*133.32/(tree.Pim);
    
    Pim = tree.Pim;%*beta;
    %     Pim/133.32
    
    tree.t = 0;
    c = 0;
    %     while tree.t < t_final
    while eps_error >  1e-6
        %         eps =  (q_seg(end) - q_target_f)/q_target_f
        %
        %                 if eps>0
        %                     eps = 0;
        %                 end
        c = c+1;
        tree.t(c+1) = tree.t(c) + tree.dt;
        
        err = 1000;
        while err > 1e-6
            for i=1:tree.N_gen
                
                mu(i) = viscosity(2*tree.R(i));
                Res(i) =8*mu(i)*(tree.Length(i))/(pi*tree.R(i)^4);
                
            end
            
            [q_seg, p_term] = SymmetricArterialTree_PinPout(tree.N_gen,p_in(j),p_out,Res);
            
            p_mid(:) = p_term(:) + (Res(:)/2).*q_seg(:);
            tau(:) = 4*mu(:).*q_seg(:)./(pi*tree.R(:).^3);
            
            
            for i = 1:tree.N_gen
                
                
                
                Tone_max(i) = 1/tree.Act_lvl(i)*0.5*myogenic_control(p_mid(i) - Pim , tree.Pmid(i) - tree.Pim, tree.ID(i));%% the  multiplier 0.5 here is to get the basal tone!
                
                A = a*1/(1+exp(-S0(i)));
                if a ==2
                    A=1;
                end
                x = mechanical_params(tree, i, A*Tone_max(i), tree.ID(i));
                
                %             Sb(i) = PF*subendo.Sbasal(i);
                
                Me = tree.Me(i);
                Mm = tree.Mt(i)*(tree.SMCtoCOL(i)/(tree.SMCtoCOL(i)+1));
                Mc = tree.Mt(i)*(1/(tree.SMCtoCOL(i)+1));
                Mk = [0.1, 0.1, 0.4, 0.4]*Mc;
                
                P = p_mid(i) - Pim;
                hm = 1/tree.rho_w;% subendo.Thickness(i)*Mm/(Mm+Mc+Me)*3/10;
                
                R1(i) = NR_iterate_fit( Me, Mk, Mm, tree.Radius(i), hm, P, x, tree.r_act(i), tree.r_h(i));
                LL(i) = R1(i)/tree.R(i);
                
                %             tree.r_act(i) = tree.r_act(i) + k_act*(R1(i)-tree.r_act(i));
                
                %         S(i) = Km*(q_seg(end) - q_target)/q_target*tree.dt + Ks*(tau(i) - tree.tau_h(i))/tree.tau_h(i)*tree.dt + Kp*(p_mid(i) - tree.Pmid(i))/tree.Pmid(i)*tree.dt ;
            end
            %             tree.R = R1;
            %             eps =  ( q_target_f - q_seg(end) )/q_target_f;
            % %             q_target_f
            % %             q_seg(end)
            %             for i = 1:tree.N_gen
            %                 [phi_s(i), phi_p(i)] = signal_weights(tree.R(i));
            %                 S0(i) =  -Km*eps - phi_s(i) * Ks *(tau(i) - tree.tau_h(i))/tree.tau_h(i) + phi_p(i) * Kp *(p_mid(i) - tree.Pmid(i))/tree.Pmid(i) + S0(i) ;%
            %             end
            %
            err = sqrt(sum(LL-1).^2);
            tree.R = R1;
            
            %         R1;
        end
        tree.R = R1;
        eps =  ( q_target_f - q_seg(end) )/q_target_f;
        %             q_target_f
        %             q_seg(end)
        for i = 1:tree.N_gen
            [phi_s(i), phi_p(i), phi_m(i)] = signal_weights(tree.R(i));
            S0(i) =  -phi_m(i)*Km*eps - phi_s(i) * Ks *(tau(i) - tree.tau_h(i))/tree.tau_h(i) + phi_p(i) * Kp *(p_mid(i) - tree.Pmid(i))/tree.Pmid(i) + S0(i) ;%
        end
        if mod(c,10) == 0
            scatter(tree.t(c),q_seg(end)/q_target_f,'.k');hold on;
            %                 scatter(tree.t(c),mean((tau(end))./tree.tau_h(end)),'.k');hold on;
            pause(0.0001);
        end
        %         for i = 1:tree.N_gen
        eps_error = abs(eps - eps0);
        eps0 = eps;
    end
    
    
    %         q_network(j,c) = q_seg(1);
    %     end
    q_network(j) = q_seg(1);
    %     tree.R./tree.Radius
    %     for i = 1:subendo.N_gen
    
    %     scatter((p_in(j))/133.32,tree.t,'ok');hold on;
    for i=1:tree.N_gen
        Active(i,j) = a*1/(1+exp(-S0(i)));
        TAU(i,j) = tau(i);
        PMID(i,j) = p_mid(i);
    end
    A
    %     p_in(j)/133.32
        
end

%     semilogx(R,p_mid/133.32,'*');
% S
LL = tree.R./tree.Radius;
% plot((p_in)/133.32,q_network./tree.q(1),marker,'color',0.1*[1 1 1]);hold on;
%
% t(end)
% c
% tau./tree.tau_h
