function [L]=LengthSegmentk(R)
% Estimation of length of each segment from Qureshi et al. (2014)
% (Empirical relations given in cm)
% symmetric tree (one representative segment per generation)
% R(k) radius of the k-generation's segment

% from cm to m
%     if (R<0.00005)
%         L = 1.79*(R*100)^(0.47)/100;
%     else
%         L = 15.75*(R*100)^1.10/100;
%     end

% from [Olufsen 2012] for 4-12 order in Huang paper
% empirical relations given in mm = 0.001m
%     L = 12.4*(R*1000)^1.1/1000;

% scale by factor half to get more physiological length (shorter)
% as seen in larger vessels
%     L = 10*LOluf;

%     L = (0.007474 *(2*R*1e6).^1.172 + 0.2235)/1000;  % From pig morphometry data for LAD, added on Nov14th
    L = 14.5*2*R + 0.2235/100;
%
% %test:
% figure()
% plot(1:25,12.4.*(Radius.*1000).^1.1./1000,'r-'); hold on
% plot(1:25,12.4.*(Radius.*1000)./1.5/1000,'b--');