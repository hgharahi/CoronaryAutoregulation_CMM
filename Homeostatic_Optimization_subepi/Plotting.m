%close all;
%% ========================== FIGURES =====================================
% Output for the presentatoins and paper
% =========================================================================
% modified on October 10 2018
%
%% ************************** homeostatic optimization results ************

% for generations 
k = 1:1:N_gen;  
%------------------------h/2R ratio --------------------------------------
h2=figure;
plot(2*Radius(k)*100,ratio(k),'o-','LineWidth',2);hold on;
for i=1:N_gen
plot(2*Radius(i)*100,HtoR(Radius(i))/2,'+b','LineWidth',1);
end
D_choy_subendo = [(24.4+48.2)/2, (48.3+101.5)/2];
HtoR_choy_subendo = [0.07, 0.05];
scatter(D_choy_subendo*1e-6*100, 2*HtoR_choy_subendo,'or','LineWidth',1);
set(gca,'XDir','reverse');xlabel({'D (cm)'},'FontSize',12);
ylabel({'h/D'},'FontSize',12);hold on;
legend({'Thickness-to-diameter ratio','Gou & Kassab 2003','Choy & Kassab 2009'},'FontSize',12,'Location','best');
grid on;
% axis([0 2*Radius(1)*100 0 0.2]);
FigModify(h2,'ratioHto2R')

% %------------------------Lenght vs gen-------------------------------------
% h3=figure;
% set(gca,'xscale'); xlim([0 N_gen]); hold on;
% plot(k, Length(k).*100,'o-','LineWidth',2);   % factor 0.5 from Olufsen
% plot(1, 34.7754*0.1,'*','LineWidth',2); 
% legend({'Linear','LAD from Mittal et al 2005'},'FontSize',12);
% xlabel({'gen.'},'FontSize',12);ylabel({'Length (cm)'},'FontSize',12);
% grid on;
% FigModify(h3,'LvsGen')
% 
% %----------------------- area ratio----------------------------------------
% %(A1+A2)/A0=2A1/A0=2R1^2/R0^2
% h4=figure; 
% set(gca,'xscale'); xlim([0 N_gen]); hold on;
% eta01=1.2; eta02=1.3;
% eta=ones(1,N_gen);
% for m=2:N_gen
%     eta(m)=2*Radius(m)^2/Radius(m-1)^2;
% end
% etaR1(1:N_gen)=eta01;etaR2(1:N_gen)=eta02;
% set(gca,'xscale'); xlim([0 N_gen]); hold on;
% plot(2:N_gen, eta(2:N_gen),'o--','LineWidth',2); 
% plot(2:N_gen,etaR1(2:N_gen),'r--','LineWidth',2);
% plot(2:N_gen,etaR2(2:N_gen),'r--','LineWidth',2);
% legend({'Optimization Results','open-end: 1.2< area ratio <1.3 [Hollander-2001]'},...
%     'FontSize',12);
% xlabel({'gen.'},'FontSize',12);ylabel({'Daughter-to-parent area ratio'},'FontSize',12);
% axis([1 N_gen 1.1 1.4]);
% grid on;
% FigModify(h4,'AreaRatio')
% 
% %-----------------------Eh/R0----------------------------------------------
% h5=figure;
% set(gca,'xscale'); xlim([0 N_gen]); hold on;
% % EhrKall(1:N_gen)=EhrK; EhrYall(1:N_gen)=EhrY; EhrQall(1:N_gen)=EhrQ;
% plot(k, Thickness(k).*YoungMod_tt(k)./RzeroP(k)./1000,'o--','LineWidth',2); 
% % plot(k, EhrKall./1000,'--','LineWidth',2); 
% % plot(k, EhrYall./1000,'-.','LineWidth',2);
% % errorbar(1,EhrEx_mean/1000,abs(min(EhrEx_err/1000)),abs(max(EhrEx_err/1000)),...
% %     'o','LineWidth',2);
% % legend({'Optimization Results','Eh/R0 [Krenz-2003]','Eh/R0 [Yen-1990]',...
% %     'E_{\theta\theta}h/R0 (MSU Experiments)'},'Location','NorthEast',...
% %     'FontSize',12);
% xlabel('gen.','FontSize',12);ylabel('E_{\theta\theta}h/R0 (kPa)',...
%     'FontSize',12);
% grid on;
% %axis([0 N_gen 2 14]);
% FigModify(h5,'StructuralStiffness')
% 
% %---------------homeostatic values vs radius (and generations)-------------
% h6=figure;
% 
% set(gca,'XDir','reverse');
% xlabel('d (cm)');ylabel('E_{\theta\theta} (kPa)');hold on;
% plot(2*Radius(k)*100,YoungMod_tt(k)/1000,'o-'); %,'LineWidth',1.5);
% FigModify(h6,'Ett')

h61=figure;
subplot(3,1,1); set(gca,'XDir','reverse');xlabel('d (cm)');ylabel('\tau (Pa)');hold on;
plot(2*Radius(k)*100,shear(k),'o-'); %,'LineWidth',1.5);

subplot(3,1,2); set(gca,'XDir','reverse');xlabel('d (cm)');ylabel('Pmid (mmHg)'); hold on;
plot(2*Radius(k)*100,Pmid(k)/133.32,'o-'); %,'LineWidth',1.5);
subplot(3,1,3); set(gca,'XDir','reverse');xlabel('d (cm)');ylabel('\sigma_h (kPa)');hold on;
plot(2*Radius(k)*100,sigma_h(k)/1000,'o-'); %,'LineWidth',1.5);
% saveas(gcf,'hemodynamics.tiff')

h62=figure;
set(gca,'xscale');xlabel('gen.');ylabel(' Radius exponent \xi');hold on;
plot(1:N_gen-1,ksi(1:N_gen-1),'o-','LineWidth',2); %,'LineWidth',1.5);
xlim([0 N_gen]);
grid on;
% FigModify(h62,'ksi')
hold off;
% 

h70=figure;
set(gca,'XDir','reverse');
set(gca,'xscale');xlabel('R');ylabel(' Activation');hold on;
plot(2*100*Radius(k),Act_lvl(k),'o-','LineWidth',2); %,'LineWidth',1.5);
grid on;
% FigModify(h62,'ksi')
hold off;