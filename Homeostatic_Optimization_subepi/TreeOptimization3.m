% modified by Vasilina, August 7, 2018
function [Mt,Me,SMCtoCOL,R,L,Table,p_mid] = TreeOptimization4(Pin,Pout,N_gen) 
% This function solves the optimization problem of the arterial tree
% Given: q_parent,p_terminal,N_gen
% Output:

%% initialization
global R0 mu 

% preallocate arrays
Mt = zeros(1,N_gen);
Me = zeros(1,N_gen);
Mtotal = zeros(1,N_gen);
SMCtoCOL = zeros(1,N_gen); 
phi_e = zeros(1,N_gen); 

R = zeros(1,N_gen);
L = zeros(1,N_gen);
p_mid = zeros(1,N_gen); 
Res = zeros(1,N_gen);   
Table = zeros(N_gen,4);

% initialization
lhat = ones(1,N_gen);
lhat0 = lhat;

% intitalize R,L,HydRes and Table
for k=1:N_gen
    R(k) = (R0^3/(2^(k-1)))^(1/3);
    L(k) = LengthSegmentk(R(k));
    Res(k) = 8*mu*(L(k))/(pi*R(k)^4);
end

% %initialization of hemodynamics
[q, p_term]= SymmetricArterialTree_PinPout(N_gen,Pin,Pout,Res); 

%VF: initialize 
for k=1:N_gen
    p_mid(k) = p_term(k)+ Res(k)*q(k)/2;
end

% p_term = Table(:,3);
Y0 = q;

%% Solve for terminal pressure to fit mass and geometry relations
% for lhat,R,L,Table

err = 1000;
c = 0; tol=1e-8;
while err > tol
    for j=1:N_gen

        % optimization of the mass at each segment
        % use previous generation result as initialization
        if j>1 
            lhat0(j) = lhat(j); 
        end
        [Mtotal(j), lhat(j), phi_e(j), SMCtoCOL(j)] = mass_optimiz_tree(p_mid(j),q(j),j,lhat0(j));   
        Mt(j) = (1 - phi_e(j))*Mtotal(j);
        Me(j) = phi_e(j)*Mtotal(j);

        % update geometry R,L
        R(j) = R0*lhat(j);
        L(j) = LengthSegmentk(R(j));

        % update resistance
        Res(j) = 8*mu*(L(j))/(pi* R(j)^4); 
%         Table(j,4) = Res(j);

    end
    % update steady hemodynamics for entire tree
    % with new resitance
    [q, p_term] = SymmetricArterialTree_PinPout(N_gen,Pin,Pout,Res); 

    % update terminal pressure
    Y = q;

    % compute residual error
    err = norm(abs(Y-Y0));

    % update reference solution
    Y0 = Y;
    c = c+1;
    
    % for output
    for j=1:N_gen
        Table(j,:) = [j, q(j), p_term(j), Res(j)];
        %VF: update October 10, 2018, moved from the beginning to the end
        %of iterations
        p_mid(j) = p_term(j)+ Res(j)*q(j)/2; 
    end
    plot(R,p_term/133.32);
    pause(0.1)
end

disp(['Tree Pressure Optimization iterations # ',num2str(c),', tolerance ',num2str(tol),]);
