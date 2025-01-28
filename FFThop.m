% FFThop.m
function Xout = FFThop(X,x,N,P)
% Compute an update to the length-N FFT.  
% X is the last FFT piece
% x consists of P points which are used in the update
% N is the length of the FFT
% P is the length of the hop (the new piece)

global subFFTs
global subFFTptr;

M = N/P;
xN = [x zeros(1,N-P)];  % create an N vector to transform
keyboard
XP = fftP(xN,P,N);      % transform the new block
phases = exp(1j*2*pi*(0:N-1)*P/N);
Xout = X - subFFTs(subFFTptr,:);
Xout = Xout + XP;
%  Xout = (X - subFFTs(subFFTptr,:) + XP);
%                ^         ^                 ^
%             current  downdate            update

Xout = Xout .* phases;
subFFTS(subFFTptr,:) = XP;        % replace the oldest sub-FFT
subFFTptr = mod(subFFTptr,M) + 1;
