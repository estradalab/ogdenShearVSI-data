function [noisy_sig, sqrtsigma2] = gNoise( I , SNR, amp )
%{
This function adds gaussian noise to the noiseless input signal to achieve 
desired SNR
sig  : noise-less input signal
SNR: desired SNR in db (ratio of peaks)
%}

IPwr = norm(I,2) / length(I);


% noise power is equal to sigma^2
%sigma2 = IPwr / 10^(SNR/10) ;
sigma2 = 1/SNR;%10^(SNR/10);

%sigma2 = IPwr / SNR ;
%sqrtsigma2 = sqrt(sigma2);

%noisy_sig = I + sqrtsigma2*randn(size(I));
noisy_sig = I + sigma2*amp*randn(size(I));

end
