% for j = 1:length(Pin)
%     j
%     er = 1000;
%     c = 0;
%     while er > 10^-6
%         
%         
%         c = c+1;
%         t(c+1) = t(c) + dt;
%         
%         for i=1:subendo.N_gen
%             
%             Res(i) =8*mu*(subendo.Length(i))/(pi*R(i)^4);
%             
%         end
%         
%         [q_seg, p_term] = SymmetricArterialTree_PinPout(subendo.N_gen,Pin(j),Pout,Res);
%         
%         p_mid(:) = p_term(:) + (Res(:)/2).*q_seg(:);
%         tau(:) = 4*mu*q_seg(:)./(pi*R(:).^3);
%         
%         %     M = metabolic_control(i, q_seg(end));
%         
%         
%         for i = 1:subendo.N_gen
%             
%            
%             
%             FF = shear_control(tau, tau_h);
%             PF = 0.5*myogenic_control(p_mid(i) - subendo.Pim , subendo.Pmid(i) - subendo.Pim, ID(i)); %% the  multiplier 0.5 here is to get the basal tone!
%             %%% Danger! Danger! PF must be 1 for the basal condition, right now, it is not!
%             %%% You have got to tune it again! 
%             x = mechanical_params(subendo, i, PF, 1, ID(i));
%             
% %             Sb(i) = PF*subendo.Sbasal(i);
%             
%             Me = subendo.Me(i);
%             Mm = subendo.Mt(i)*(subendo.SMCtoCOL(i)/(subendo.SMCtoCOL(i)+1));
%             Mc = subendo.Mt(i)*(1/(subendo.SMCtoCOL(i)+1));
%             Mk = [0.1, 0.1, 0.4, 0.4]*Mc;
%             P = p_mid(i) - subendo.Pim;
%             hm = 1/subendo.rho_w;% subendo.Thickness(i)*Mm/(Mm+Mc+Me)*3/10;
%             
%             R1(i) = NR_iterate_fit( Me, Mk, Mm, subendo.Radius(i), hm, P, x, r_act(i), r_h(i));
%             LL(i) = R1(i)/R(i);
%             
%             r_act(i) = r_act(i) + k_act*(R1(i)-r_act(i));
%             
%         end
%         
%         er = sqrt(sum(LL-1).^2);
%         R = R1
%         
%     end
%     %     for i = 1:subendo.N_gen
%     scatter((Pin(j))/133.32,q_seg(1)/subendo.q(1),'ok','MarkerFaceColor',i/subendo.N_gen*[0 0 1]);hold on;
%     %     end
% end
% %     semilogx(R,p_mid/133.32,'*');
% LL = R./subendo.Radius
% %
% % t(end)
% % c
% % tau./tau_h











