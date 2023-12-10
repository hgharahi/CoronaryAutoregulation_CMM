clear; clear -global; close all; clc;

rho_w = 1060;
subepi = load('subepi_Umich.mat');
subendo = load('subendo_Umich.mat');
%% Data from Liao and Kuo

D100 = [275 175 108 69.5]*1e-6; % Diameters @ p = 100 mmHg
% D_h = [255 155 75 50]*1e-6;
% D_h = [D_h D_h];

tree_idx = [1 4 7 11];

% subepi.Gc(tree_idx(3)) = 0.97;
% subendo.Gc(tree_idx(3)) = 0.97;

subepi.Gc(tree_idx(4)) = 0.98;
subepi.Gc(tree_idx(3)) = 1.00;

subendo.Gc(tree_idx(4)) = 0.97;
subendo.Ge1(tree_idx(4)) = 1.02;
subendo.Ge2(tree_idx(4)) = 1.02;

P_h(1:4) = subepi.Pmid(tree_idx);
D_h(1:4) = 2*subepi.Radius(tree_idx);
P_h(5:8) = subendo.Pmid(tree_idx);
D_h(5:8) = 2*subendo.Radius(tree_idx);

Pim(1) = subepi.Pim; % 13 mmHg for subepi, 48 mmHg for subendo
Pim(2) = subendo.Pim; % 

R_h = D_h/2;
R_act = R_h;

ID = {'B','C','D','E'}; % Artery, Large Arteriole, Intermediate Arteriole, Small Arteriole respectively

pressure_data(:,1) = xlsread('myogenic_response',3,'A2:A8')*133.32; % Pa

for i = 1:4
    
    diameter_data_active(:,i) = xlsread('myogenic_response',1,[ID{i},'2:',ID{i},'8']);
    
    diameter_data_passive(:,i) = xlsread('myogenic_response',3,[ID{i},'2:',ID{i},'8']);
    
    plot( pressure_data , diameter_data_active(:,i)); hold on;
    plot( pressure_data , diameter_data_passive(:,i));
end

%% One vessel estimation


Dmm = D_h*10^3; % meters to mm to use Gou and Kassab's function
H = (41.1*Dmm + 3.2)/1e6; %% From Gou and Kassab 2003

phi_f = 0.7; % percent occupied by interestetial fluid in arterial wall

fun = @(x)objective_all(x, subepi, subendo, tree_idx, pressure_data, diameter_data_active, diameter_data_passive, R_h, H, P_h, Pim);

x0 = [0.2 0.5 0.2 0.5 0.2 0.5 0.2 0.5 ...
      0.2 0.5 0.2 0.5 0.2 0.5 0.2 0.5 ...
      2.0 40*133.32];
options = optimset('MaxFunEvals',16000,'MaxIter',16000,'TolFun',1e-5,'TolX',1e-5,'Display','iter');
y = fminsearch(fun, x0, options)

close all;
beta = y(end-1); 
p0 = y(end);

j = 1;

for c=1:4
    
    phi(1) = y(2*j-1);
    phi_e(j) = phi(1);
    
    phi(2) = y(2*j);
    phi_c(j) = phi(2);
    
    phi_m(j) = 1 - phi(1) - phi(2);
    
    [A(j), X(:,j)] = eval_single(phi, c, subepi, tree_idx(c), pressure_data, diameter_data_active, diameter_data_passive, R_h(j), H(j), P_h(j), Pim(1), beta, p0);
    
    j = j+1;
    
end


names = {'A','LA','IA','SA'};

h7 = figure;

epi_e = [phi_e(1);phi_e(2);phi_e(3);phi_e(4)];
epi_c = [phi_c(1);phi_c(2);phi_c(3);phi_c(4)];
epi_m = [phi_m(1);phi_m(2);phi_m(3);phi_m(4)];
plot(100*epi_e,'-^','LineWidth',2.0);hold on
plot(100*epi_c,'-o','LineWidth',2.0);
plot(100*epi_m,'-s','LineWidth',2.0,'color',[0 12/12*114/255 (0)/12*189/255]);


set(gca,'xtick',[1:4],'xticklabel',names,'Fontsize',12);
ylabel('% mass fractions','Fontsize',12);
axis([0.5 4.5 0 80]);
legend('elastin','collagen','SMC','location','northwest','Fontsize',12);
title('subepicardial vessels','Fontsize',12);
print(h7, '-dmeta', ['subenpi_fractions','.emf']);

for c=1:4
    
    phi(1) = y(2*j-1);
    phi_e(j) = phi(1);
    
    phi(2) = y(2*j);
    phi_c(j) = phi(2);
    
    phi_m(j) = 1 - phi(1) - phi(2);
    
    [A(j), X(:,j)] = eval_single(phi, c, subendo, tree_idx(c), pressure_data, diameter_data_active, diameter_data_passive, R_h(j), H(j), P_h(j), Pim(2), beta, p0);
    
    j = j+1;
    
end

h8 = figure;

endo_e = [phi_e(5);phi_e(6);phi_e(7);phi_e(8)];
endo_c = [phi_c(5);phi_c(6);phi_c(7);phi_c(8)];
endo_m = [phi_m(5);phi_m(6);phi_m(7);phi_m(8)];
plot(100*endo_e,'-^','LineWidth',2.0);hold on
plot(100*endo_c,'-o','LineWidth',2.0);
plot(100*endo_m,'-s','LineWidth',2.0,'color',[0 12/12*114/255 (0)/12*189/255]);


set(gca,'xtick',[1:4],'xticklabel',names,'Fontsize',12);
ylabel('% mass fractions','Fontsize',12);
axis([0.5 4.5 0 80]);
legend('elastin','collagen','SMC','location','northwest','Fontsize',12);
title('subendocardial vessels','Fontsize',12);
print(h8, '-dmeta', ['subendo_fractions','.emf']);








