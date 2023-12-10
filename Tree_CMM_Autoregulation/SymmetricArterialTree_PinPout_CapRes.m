% modified by Vasilina on August 2, 2018
function [q, p_term] = SymmetricArterialTree_PinPout_CapRes(N_gen,Pin,Pout,Res)
% Find steady flow and terminal pressure at each generation
% symmetric bifurcation: unknowns: terminal pressure after each generation,
% and the first flow, since the rest can be found using the symmetry qk=q1/2^(k-1). The
% equations are written in terms of inlet pressure at each generation which
% is equivalent to the terminal pressure of the previous generation.

p_term = zeros(N_gen+1,1);
q = zeros(N_gen+1,1);

A = zeros(N_gen+3,N_gen+3);

for k =  1:N_gen+1
    A(k+1,1) = Res(k)/2^(k-1);
    A(k+1,k+1) = -1;
    A(k+1,k+2) = 1;
end

A(1,2) = 1;
A(k+2,k+2) = 1;

c = [Pin;zeros(N_gen+1,1);Pout];

b = A\c;

q(1) = b(1);
for k = 1:N_gen+1
   
   q(k) = b(1)/(2^(k-1));
   p_term(k) = b(k+2);
    
end

p_term(N_gen+1) = [];
q(N_gen+1) = [];
    
    
    
    
