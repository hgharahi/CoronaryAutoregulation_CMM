% modified by Vasilina. April 23, 2018
function [Mt,Me,R,L,Table,p_mid] = TreeCMM(Q_parent,p_terminal,N_gen) 
% This function solves the optimization problem of the arterial tree
% Given Q_parent,p_terminal,phi_e,N_gen

%% initialization
global R0 mu SMCtoCOL phi_e

% preallocate arrays
Mt = zeros(1,N_gen);
Me = zeros(1,N_gen);
Mtotal = zeros(1,N_gen);
R = zeros(1,N_gen);
L = zeros(1,N_gen);
p_mid = zeros(1,N_gen); 
Res = zeros(1,N_gen);   

% initialization
lhat = ones(1,N_gen);

% intitalize R,L,HydRes and Table
for k=1:N_gen
    R(k) = R0*lhat(k);
    L(k) = LengthSegmentk(R(k));
    Res(k) = 8*mu*(L(k))/(pi*R(k)^4);
end

% %initialization of hemodynamics
Table = SymmetricArterialTree2(Q_parent,p_terminal,N_gen,Res); 
    
Y0 = Table(:,3); % terminal pressures

%% Solve for terminal pressure to fit mass and geometry relations
% for lhat,R,L,Table

err = 1000;
c = 1;
while err > 1e-6
    for j=1:N_gen
        % initialization for each segment (generation)
        q_seg = Table(j,2);
        Res(j)= Table(j,4);
        %average pressure = (pTermSteady(k)+ pInpSteady(k))/2
        p_mid(j) = Table(j,3)+ Res(j)/2*q_seg; 

        % optimization of the mass
        [Mtotal(j), lhat(j), phi_e(j), SMCtoCOL(j)] = single_seg_CMM(p_mid(j),q_seg,j);   

        % update geometry R,L
        R(j) = R0*lhat(j);
        L(j) = LengthSegmentk(R(j));

        % update resistance
        Res(j) = 8*mu*(L(j))/(pi* R(j)^4); 
        Table(j,4) = Res(j);
    end
    % get steady hemodynamics
    Table = SymmetricArterialTree2(Q_parent,p_terminal,N_gen,Res); 

    % update terminal pressure
    Y = Table(:,3);

    % compute residual error
    err = norm(abs(Y-Y0));

    % update reference solution
    Y0 = Y;
    c = c+1;
end

disp(c);
