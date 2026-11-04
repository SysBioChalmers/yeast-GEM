function [rhat, ess] = mcmcDiag(x)
% mcmcDiag  Split-R-hat and effective sample size for one quantity.
%   x is draws-by-chains. Each chain is split in half (Gelman et al.,
%   Bayesian Data Analysis, 3rd ed.), and the effective sample size uses the
%   initial positive sequence of the combined autocorrelation (Geyer 1992).
%   A quantity that does not vary returns rhat = 1 and ess = number of draws.

[n, m] = size(x);
h = floor(n/2);
x = [x(1:h,:), x(h+1:2*h,:)];
[n, m] = size(x);
if all(abs(x(:) - x(1)) < 1e-12)
    rhat = 1; ess = n*m; return
end
chainMean = mean(x, 1);
W = mean(var(x, 0, 1));
B = n * var(chainMean);
varPlus = (n-1)/n*W + B/n;
rhat = sqrt(varPlus / W);

% Autocorrelation averaged over chains, via FFT.
xc = x - chainMean;
nfft = 2^nextpow2(2*n);
f = fft(xc, nfft);
acov = real(ifft(f .* conj(f)));
acov = acov(1:n, :) ./ n;
rho = 1 - (W - mean(acov, 2)) / varPlus;
% Geyer: sum consecutive pairs while positive.
tau = -1;
for t = 1:2:n-1
    pairSum = rho(t) + rho(t+1);
    if pairSum < 0, break; end
    tau = tau + 2*pairSum;
end
ess = n*m / tau;
end
