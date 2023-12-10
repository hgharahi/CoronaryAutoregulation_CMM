% Flow and pressure in freq. domain are Hermitian functions, i.e should be
% conjugate symmetric. Thus only the right half of spectrum is considered
% with first component related to zero frequency.
%  - Here N = Nf-1 to exclude repeating boundary point
%  - Scale by N to get steady flow at Q(1)
%  - no factor 2 to use with ifft(Q,'symmetric') - self-conjugate Q
%
% (by Vasilina, 6-19-2018)

function [Qn] = FilterAndFrequencyDomain(f,NumModes,N)
   
    % Obtain the discrete Fourier Transform 
    Fn=fft(f(1:N)); 

    % flow coeffs
    Qn = zeros(NumModes,1);
    
    % scale by the number of points so that the magnitude does not depend 
    % on the length of the signal
    % zero-frequency, steady flow
    Qn(1) = Fn(1)/N; 
    
    %(2:NumModes) this is related to just the first half of nonzero 
    % frequencies since the second half is a mirror image of the first;
    % no need of factor 2 as it will be used for symmetrical conjugate
    % in ifft
    for k=2:NumModes
        Qn(k) = Fn(k)/N;   
    end
    
end