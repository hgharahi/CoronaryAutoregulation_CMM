function ParamEst_plots(q_network, A, R, SMAX, TAU, PMID, S, p_in, tree)

markers = {'-','--','-.','.-','.','+','x'};
TYPEs = {'Small Artery';'Large Arteriole';'Intermediate Arteriole';'Small Arteriole'};


h1 = figure;
c = 0;
for i=1:3:tree.N_gen
    c = c + 1;
    plot(p_in/133.32,TAU(i,:)/tree.tau_h(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = [TYPEs{c},' (',num2str(round(2*tree.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('$\tau/\tau_{h}$','Interpreter','latex','FontSize',16); 
title(['S',tree.name(2:end),'cardial Tree'],'FontSize',12); grid on;
box on;
print(h1, '-dmeta', ['tau_',tree.name,'.emf']);
movefile(['tau_',tree.name,'.emf'],['Figs\tau_',tree.name,'.emf']);

h2 = figure;
c = 0;
for i=1:3:tree.N_gen
    c = c + 1;
    plot(p_in/133.32,A(i,:),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} =  [TYPEs{c},' (',num2str(round(2*tree.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('A','FontSize',12); 
title(['S',tree.name(2:end),'cardial Tree'],'FontSize',12); grid on;
box on;
axis([20 160 0 1]);
print(h2, '-dmeta', ['A_',tree.name,'.emf']);
movefile(['A_',tree.name,'.emf'],['Figs\A_',tree.name,'.emf']);

h3 = figure;
tau_leg = {};
c = 0;
for i=1:3:tree.N_gen
    c = c + 1;
    plot(p_in/133.32,R(i,:)/tree.Radius(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = [TYPEs{c},' (',num2str(round(2*tree.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('$D/D_{h}$','Interpreter','latex','FontSize',12); 
title(['S',tree.name(2:end),'cardial Tree'],'FontSize',12); grid on;
box on;
axis([20 160 0.6 1.2]);
print(h3, '-dmeta', ['D_',tree.name,'.emf']);
movefile(['D_',tree.name,'.emf'],['Figs\D_',tree.name,'.emf']);

h4 = figure;
tau_leg = {};
c = 0;
for i=1:3:tree.N_gen
    c = c + 1;
    plot(p_in/133.32,PMID(i,:)/tree.p_mid(i),markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} = [TYPEs{c},' (',num2str(round(2*tree.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('p/p_{h}','FontSize',12); 
title(['S',tree.name(2:end),'cardial Tree'],'FontSize',12); 
grid on;
box on;
print(h4, '-dmeta', ['Pmid_',tree.name,'.emf']);
movefile(['Pmid_',tree.name,'.emf'],['Figs\Pmid_',tree.name,'.emf']);

h5 = figure;
tau_leg = {};

c = 0;
for i=1:3:tree.N_gen
    c = c + 1;
    plot(p_in/133.32,SMAX(i,:)/1e6,markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} =  [TYPEs{c},' (',num2str(round(2*tree.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('S_p (MPa)','FontSize',12); 
title(['S',tree.name(2:end),'cardial Tree'],'FontSize',12); grid on;
box on;
print(h5, '-dmeta', ['SMAX_',tree.name,'.emf']);
movefile(['SMAX_',tree.name,'.emf'],['Figs\SMAX_',tree.name,'.emf']);

h6 = figure;
tau_leg = {};

c = 0;
for i=1:3:tree.N_gen
    c = c + 1;
    plot(p_in/133.32,S(i,:)/1e6,markers{c},'color',[0 i/12*114/255 (12-i)/12*189/255],'LineWidth',1.5);hold on;
    tau_leg{c} =  [TYPEs{c},' (',num2str(round(2*tree.Radius(i)*1e6)),' \mum)'];
end
legend(tau_leg,'Location','Best','FontSize',12); 
xlabel('p_{in} (mmHg)','FontSize',12);  ylabel('S (MPa)','FontSize',12); 
title(['S',tree.name(2:end),'cardial Tree'],'FontSize',12); grid on;
box on;
print(h6, '-dmeta', ['Sp_',tree.name,'.emf']);
movefile(['Sp_',tree.name,'.emf'],['Figs\Sp_',tree.name,'.emf']);



