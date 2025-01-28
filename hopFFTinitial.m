function X = hopFFTinitial(N,P)
% function X = hopFFTinitial(N,P)

% TKM Jan 2025

global subFFTs


M = N/P;
X = zeros(1,N,'like',1+1i);

klist = (0:N-1);
for i=0:M-1
   phasefactor = exp(-1i*2*pi*i*P*klist/N);
   X = X + phasefactor .* subFFTs(i+1,:);
end
