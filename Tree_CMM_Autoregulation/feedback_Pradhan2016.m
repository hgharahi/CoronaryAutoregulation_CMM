function [F1, gfb1, gol] = feedback_Pradhan2016(M)
%% this sscript numerically finds the feedback resposne to the oxygen metabolism in the myocites
syms x1 x2 x3 x4
Vc = 0.04;  %capillary volume density in the myocardium ml/g
H = 0.45;    %Hematocrit level
Sa = 0.955; %Arterial oxygen saturation
Co = 476;   %Oxygen carrying capacity of red blood cells micro-litre O2/ml

J0 = 283.4; %Basal rate of ATP release into plasma, microM/min
S0 = 0.1; %Basal Oxygen saturation
Ta = 28.15; %Arterial plasma ATP saturation nM
Km = 22.77; %parameter relating sympathetic pathway response and myocardial oxygen metabolism mico-litre O2/min/micro-g/nM
M0 = 33.54; %myocardial oxygen metabolism at zero N.
Kmyo = 0;
Kol = 1.8*1e-3;

Kfb = 2.3507*1e-3;
g0 = 0.63*1e-3;

% N = 10;
% M = 170; 
P = 100;

% M = Km*N + M0
N = (M - M0)/Km;
gol = Kol*N;

F = x1;
Sv = x2;
Tv = x3;
gfb = x4;

assume(Sv, 'positive')
[F1,Sv1,TV1,gfb1] = vpasolve([0 == Tv - ( Ta + ( Vc*H*J0*S0 )/( F*(Sa-Sv) )* exp(-Sa/S0) *( exp((Sa-Sv)/S0)-1 )), ...
                        0 == gfb - Kfb * ( Tv-Ta ),...
                        0 == Sv - (Sa - M/( F*Co*H )),...
                        0 == F - (gfb + gol + g0)*P] , [F,Sv,Tv,gfb]);

