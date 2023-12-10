function [Cost,q,p_term,Mtotal,Table,R,SMCtoCOL,Me,Mt,p_mid,ksi,L] = evaluate(lhat)

global R0 Pin Pout N_gen
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


% intitalize R,L,HydRes and Table
for k=1:N_gen
    R(k) = R0*lhat(k);
    L(k) = LengthSegmentk(R(k));
    mu(k) = viscosity(2*R(k));
    Res(k) = 8*mu(k)*(L(k))/(pi*R(k)^4);
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
Cost = 0;
while err > tol
    for j=1:N_gen

        [C, Mtotal(j)] = mass_optimiz_tree_2(p_mid(j),q(j),j,lhat(j));
        [phi_e(j), phi_c, phi_m] = mass_fracs_2(2*lhat(j)*R0);
        SMCtoCOL(j) = phi_m/phi_c;
        Mt(j) = (1 - phi_e(j))*Mtotal(j);
        Me(j) = phi_e(j)*Mtotal(j);

        % update geometry R,L
        R(j) = R0*lhat(j);
        L(j) = LengthSegmentk(R(j));
        mu(j) = viscosity(2*R(j));
        % update resistance
        Res(j) = 8*mu(j)*(L(j))/(pi* R(j)^4); 
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
    Cost = Cost + 2^(j-1)*C*L(j);
end

for k=1:N_gen-1
    ksi(k) = 1/(log2(R(k))-log2(R(k+1))); 
end

