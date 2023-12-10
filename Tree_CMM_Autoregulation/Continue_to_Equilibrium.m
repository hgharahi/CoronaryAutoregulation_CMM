%% G&R with Hemodynamics update
t_final =  350 + t(c);
while t < t_final
       
    t(c+1) = t(c) + dt;
    c = c+1;
    
    for j=1:N_gen
        R(j) =r(c-1,j);
        L(j) = LengthSegmentk(R(j));
        Res(j) = 8*mu*(L(j))/(pi*R(j)^4);
    end
    Table = SymmetricArterialTree2(Q_parent,p_terminal,N_gen,Res);
    
    q_seg(:) = Table(:,2);
    Res(:)= Table(:,4);

%     disp([num2str(t(c-1)),'    ',num2str(r(c-1,1)),'    ',num2str(r(c-1,2))]);

    for j=1:N_gen
        
%         j
%         
        tau(j) = 4*mu*q_seg(j)/(pi*r(c-1,j)^3);
        
        p_mid(j) = Table(j,3)+ Res(j)/2*q_seg(j);
        
       
        err = 1000;
        num_it = 1;
        while err > 10^-6 && num_it < 5
            
            [r(c,j), sigma_k, sigma_m, T_act] = NR_iterate( Me(1,j),Mm(1,j),Mc(1,j),pk1_a(:,j),pk2_a(:,j),pk3_a(:,j),pm_a(:,j),r(1,j),p_mid(j),dt,r_p(:,j),r_act(j), j);
            Lt = r(c,j)/R(j);
            
            for k=1:4
                mPk(k) = f_mP(sigma_k(k), sigma_h_c, kg, tau(j), tau_h, kg_sh, m_basal_c(k,j));
            end
            
            mPm = f_mP(sigma_m, sigma_h_m, kg, tau(j), tau_h, kg_sh, m_basal_m(j));
            
            [Mk_new(:,j), Mm_new(:,j), pk1_a(:,j), pk2_a(:,j), pk3_a(:,j), pm_a(:,j)] = ...
                age_distribution(pk1_p(:,j), pk2_p(:,j), pk3_p(:,j), pm_p(:,j), mPk, mPm, num_pa, dt, kq_c);

            err = sqrt(sum((Mk_new(:,j)-Mk(:,j)).^2)/sum(Mk(:,j).^2) + (Mm_new(:,j)-Mm(:,j))^2/Mm(:,j)^2);
            Mk(:,j) = Mk_new(:,j);
            Mm(:,j) = Mm_new(:,j);
            Mc(:,j) = sum(Mk(:,j));
            total_M = (Me(j)+Mm(j)+Mc(j));
            h(j) = total_M / ((1-phi_f)*rho_w*Lz*Lt);
            hm(j) = Mm(j)/((1-phi_f)*rho_w*Lt*Lz);
            
            num_it = num_it + 1;
        end
        
        for k=1:4
            mPk_p(k) = mPk(k);
            sigma_k_p(k) = sigma_k(k);
        end
        mPm_p = mPm;
        sigma_m_p = sigma_m;
        
        [pk1_p(:,j), pk2_p(:,j), pk3_p(:,j), pm_p(:,j), r_p(:,j)] = time_evolv(pk1_a(:,j), pk2_a(:,j), pk3_a(:,j), pm_a(:,j), r_p(:,j), r(c,j));
        
        
        
        r_act(j) = r_act(j) + k_act*(r(c,j)-r_act(j))*dt;
    end
    
    if mod(c-1,5) == 0
        clf;
        for j=1:N_gen
            
            scatter(t,r(:,j)/(R0^3/(2^(j-1)))^(1/3),'.'); hold on;
            pause(0.001);
            %         disp(num2str(er))
        end
        saveas(gcf,'radius.png');
        %                     figure(3);
        %                     scatter(t(c),T_act/hm,'.r'); hold on;
        %                     pause(0.1)
        %                     figure(4);
        %                     scatter(t(c),P*r(c)/h,'.r'); hold on;
        %                     pause(0.1)
    end
end


