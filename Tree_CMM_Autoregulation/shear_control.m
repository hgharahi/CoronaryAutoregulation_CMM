function F = shear_control(tau, tau_h, ID, Act_lvl)

% F0 = 2*(1-Act_lvl);
% Kt = tau_h * Act_lvl/(1 - Act_lvl);
% F = (1 - Fmax*(tau/(Kt+tau)));

if ID == 'B'
    s = 1.00;
    Kt = tau_h * (Act_lvl/(1 - Act_lvl))^(1/s);
    F = (1 - (1/((Kt/tau)^s+1)));
elseif ID == 'C'
    s = 1.00;
    Kt = tau_h * (Act_lvl/(1 - Act_lvl))^(1/s);
    F = (1 - (1/((Kt/tau)^s+1)));
elseif ID == 'D'
    s = .50;
    Kt = tau_h * (Act_lvl/(1 - Act_lvl))^(1/s);
    F = (1 - (1/((Kt/tau)^s+1)));
else
    s = 0.0;
%     Kt = tau_h * (Act_lvl/(1 - Act_lvl))^(1/s);
    F = Act_lvl;%(1 - (1/((Kt/tau)^s+1)));
end

