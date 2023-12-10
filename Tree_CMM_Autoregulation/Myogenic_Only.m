function [q_network, A, R, SMAX, TAU, PMID, S, SIGMA_H] = Myogenic_Only(K, tree)


global p_in  p_out  M a f_in

Km = K(1);
Kp = K(2);
Ks = K(3);

mu = tree.mu;
Pim = tree.Pim;

% This loop : 1) finds the equivalent basal stimulus given the activation.
% 2) finds the weights corresponding to each of the control mechanism. 3)
% computes the resistance in each vessel given the homeostatic sizes of
% vessels.

for i=1:tree.N_gen
    
    s0(i) = ( - log( 1/tree.Act_lvl(i) - 1 ));
    
    [phi_s(i), phi_p(i), phi_m(i)] = signal_weights(tree.R(i));
        
    Res(i) =8*mu(i)*(tree.Length(i))/(pi*tree.R(i)^4);
    
end

Res(tree.N_gen+1) = 2*tree.Cap_res; %% adds the capillary resistance (R  boundary condition)

[q_seg, p_term] = SymmetricArterialTree_PinPout_CapRes(tree.N_gen,tree.Pin,p_out,Res);

Res(end) = []; %% removes the capillary resistace for fast computing the pressure and shear stress in each segment of the tree
p_mid(:) = p_term(:) + (Res(:)/2).*q_seg(:);
tau(:) = 4*mu(:).*q_seg(:)./(pi*tree.R(:).^3);

tree.q_incr = M/70;
q_target_f = tree.q(end) * tree.q_incr;

q_target = q_seg(end);


for j = 1 :length(p_in)
    
    eps_error = 10;
       
    err = 1000;
    while err > 1e-5
        
        
        % Computes the viscosity and
        % resistance iven the current geometry of the model.
        for i = 1:tree.N_gen
            
           
            mu(i) = viscosity(2*tree.R(i));
            Res(i) =8*mu(i)*(tree.Length(i))/(pi*tree.R(i)^4);
            
        end
        
        Res(tree.N_gen+1) = 2*tree.Cap_res; %% adds the capillary resistance (R  boundary condition)
        
        [q_seg, p_term] = SymmetricArterialTree_PinPout_CapRes(tree.N_gen,p_in(j), p_out ,Res);
        
        Res(end) = []; %% removes the capillary resistace for fast computing the pressure and shear stress in each segment of the tree
        
        p_mid(:) = p_term(:) + (Res(:)/2).*q_seg(:);
        tau(:) = 4*mu(:).*q_seg(:)./(pi*tree.R(:).^3);
        

        eps = (q_seg(end)-q_target_f)/q_target_f; % finds the deviation of terminal flow from its homeostatic baseline
%                 if eps > 0
%                     eps = 0;
%                 end
        
        % Computes the stimuli coming from each mechanism and adds them all
        % up (plus the baseline stimulus) to find the stimulus given the
        % current hemodynamic state. 
        
        for i = 1:tree.N_gen
            
            s_m(i) =  phi_m(i)*( Km * eps );
            s_tau(i) = - phi_s(i) * (Ks *(tau(i) - tree.tau_h(i))/tree.tau_h(i));
            s_p(i) = + phi_p(i) * (Kp *(p_mid(i) - tree.Pmid(i))/tree.Pmid(i) );
            
            s_total(i) = s_m(i) + s_tau(i) + s_p(i) + s0(i) ;
            
%             A(i,j) = a*1/(1+exp(-s_total(i)));
            A(i,j) = 1;
            
            P = p_mid(i) - Pim;
            
            Smax = Pressure_Dependent_Tension(P, tree.A.Smax, tree.A.beta, tree.A.p0);
            
            
            if a == 2
                A(i,j) = 1;
            end
            
            x = mechanical_params(tree, i, A(i,j)*Smax, tree.ID(i));
                       
            Me = tree.Me(i);
            Mm = tree.Mt(i)*(tree.SMCtoCOL(i)/(tree.SMCtoCOL(i)+1));
            Mc = tree.Mt(i)*(1/(tree.SMCtoCOL(i)+1));
            Mk = [0.1, 0.1, 0.4, 0.4]*Mc;
            
            hm = 1/tree.rho_w;% subendo.Thickness(i)*Mm/(Mm+Mc+Me)*3/10;
            
            flg = 0;
            R1(i) = NR_iterate_fit( Me, Mk, Mm, tree.Radius(i), hm, P, x, tree.R(i), tree.r_h(i), flg, tree, i, tree.ID);
            LL(i) = R1(i)/tree.R(i);
            
            
        end
        
        err = sqrt(sum(LL-1).^2);
        tree.R = 0.1*R1 + 0.9*tree.R;
        
    end
    tree.R = R1;
    
    for i=1:tree.N_gen
        
        
        R(i,j) = R1(i);
        P = p_mid(i) - Pim;
        
        SMAX(i,j) = Pressure_Dependent_Tension(P, tree.A.Smax, tree.A.beta, tree.A.p0);
        
        S(i,j) = A(i,j)*SMAX(i,j);
        TAU(i,j) = tau(i);
        PMID(i,j) = p_mid(i);
        
        SIGMA_H(i,j) = p_mid(i)*tree.R(i)/(tree.Thickness(i))*(R(i,j)/tree.R(i))^2;
        
    end
    q_network(j) = q_seg(end);
    
end





