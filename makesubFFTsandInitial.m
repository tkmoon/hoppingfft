function X = makesubFFTsandInitial(x,N,P)
% 
% x = initial data, x(1:N)
% N: length of initital data (in case whole vector is used
% P: length of hop

global subFFTs
global subFFTptr
global phases

M = N/P;
subFFTs = zeros(M,N,'like', 1 + 1i);

subFFTptr = 1;  % initial state

for i=1:M   % build the array of subFFTs
   xi = zeros(1,N,'like',1+1i);
   xi(1:P) = x((i-1)*P + (1:P) );
   subFFTs(i,:) = fftP(xi,P,N);
end % for i


X = zeros(1,N,'like',1+1i);

phases = -1i*2*pi*P*(0:N-1)/N;
for i=0:M-1
   X = X + exp(phases*i) .* subFFTs(i+1,:);
end
phases = exp(-phases);