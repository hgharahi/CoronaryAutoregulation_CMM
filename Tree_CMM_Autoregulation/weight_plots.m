close all;clear;clc;

figure;
R = logspace(-3.6,-5,10000);

for i=1:length(R)
[phi_s(i), phi_p(i), phi_m(i)] = signal_weights(R(i));
end

h1 = figure(1);
marker = '-';
plot(2*R*1e6,phi_s,marker,'color',0.0*[1 1 1],'LineWidth',1.5);
hold on;
marker = '--';
plot(2*R*1e6,phi_m,marker,'color',0.0*[1 1 1],'LineWidth',1.5);



xlabel('D (\mum)','FontSize',12);
ylabel('\phi','FontSize',12);
legend('\phi_\tau','\phi_m','Location','best','FontSize',12);
set(gca,'FontSize',12);
axis([0, 2*max(R)*1e6 0 1.1])
set(gca,'XDir','reverse');
% Inset = imread('DiameterDependentResponse.PNG');
% ax2 = axes('Position',[.1 .25 .45 .45]);
% imshow(Inset)
% grid on;
box on;
print(h1, '-dmeta', ['weights','.emf']);
movefile('weights.emf','Figs\weights.emf');