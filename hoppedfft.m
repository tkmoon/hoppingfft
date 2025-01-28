
% Implement and test the hopped FFT
%
% TKM Jan 2025

clear subFFTs subFFTptr;

global subFFTs     % the set of Sub-FFTs
global subFFTptr   % index which is the oldest sub-FFT
global phasesx

N = 128;   % block size
P = 32;    % size of skip
M = N/P;

extrablocks = 40;
nblocks = M + extrablocks;
x = complex(1:nblocks*P, 1:nblocks*P);   % test data

% Initialize -- make subFFTs and subFFTptr, and make the initial transform
X = makesubFFTsandInitial(x(1:N),N,P);

% and compute the FFT to compare:
Xffttrue = fft(x(1:N));
err = norm(X - Xffttrue);
fprintf('initial: error=%g\n',err);

% initial:   x(1:N)     (N points)
% first hop  New data:  x(N+1:N+P)        (P points)   Overall: x(P+1:N + P)  % N points
% second hop New data:  x(N+P+1:N + 2P)   (P points)   Overall: x(2P + 1:2P + N )% N points
% etc.

for i=1:extrablocks
   X = FFThop(X,x(N + (i-1)*P + 1: N + i*P),N,P);
   Xffttrue = fft(x(i*P + 1: N + i*P));
   err = norm(X - Xffttrue);
   fprintf('i=%d  error=%g\n',i,err);
end

return

