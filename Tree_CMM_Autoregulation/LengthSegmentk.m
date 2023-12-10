function [L]=LengthSegmentk(R)
% Estimation of length of each segment from Qureshi et al. (2014)
% (Empirical relations given in cm)
% symmetric tree (one representative segment per generation)
% R(k) radius of the k-generation's segment
    % VF: correction for lhat(i,k)*100

    % from cm to m
%     if (R<0.00005)
%         L = 1.79*(R*100)^(0.47)/100;
%     else
%         L = 15.75*(R*100)^1.10/100;
%     end

    % from [Olufsen 2012] for 4-12 order in Huang paper
    % empirical relations given in mm = 0.001m
    L = 42*(R*1000)/1000;
%     L = LOluf; %this is a test

%     L = (0.008069 *(R*1e6)^1.178 - 0.1508)/1000;  % From pig morphometry data for LAD, added on Oct 31st
end
% 
% %test:
% figure()
% plot(1:25,12.4.*(Radius.*1000).^1.1./1000,'r-'); hold on
% plot(1:25,12.4.*(Radius.*1000)./1.5/1000,'b--');