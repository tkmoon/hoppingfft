function X = makesubFFTs(x,N,P)
% 
% x = initial data, x(1:N)
% N: length of initital data (in case whole vector is used
% P: length of hop

global subFFTs
global subFFTptr

M = N/P;
subFFTs = zeros(M,N,'like', 1 + 1i);

subFFTptr = 1;  % initial state

for i=1:M   % build the array of subFFTs
   xi = zeros(1,N,'like',1+1i);
   xi(1:P) = x((i-1)*P + (1:P) );
   subFFTs(i,:) = fftP(xi,P,N);
   % check
   ffttrue = fft(xi);
   fprintf('makesub i=%d error=%g\n',i,norm(subFFTs(i,:) - ffttrue));
end % for i
