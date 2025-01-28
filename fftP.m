function [x,numCmult,numCadd] = fftP(x,P,Pout)
% function [x,numCmult,numCadd] = fftP(x,P,Pout)
% given x (length N), compute the N-point FFT
%  of the first P samples, as if the P+1:N 
% points were 0.
% Return Pout frequency values.  Default: Pout = P

% TKM 1/21/25

   if(nargin == 2)
      Pout = P;
   end

   N = length(x);
   M = log2(N);

   log2P = log2(P);
   NV2 = N/2;
   NM1 = N-1;
   Q = N/P;

   % Bit reverse shuffle
   J = 1;
   for I=1:NM1
      if(I < J)
         T = x(J);
         x(J) = x(I);
         x(I) = T;
      end
      K = NV2;   % label 5
      while(1)   % label 6
         if(K < J)
            J = J-K;
            K = K/2;
            continue;
         end
         break;
      end
      J = J + K;  % label 7
   end %--------------------------------------------------

   Lfirstnoskip = log2(N/P);  % the first value of L that does not skip any blocks

   numCmult = 0; numCadd = 0;

   L = Lfirstnoskip;
   G = 2^L;

   % Input stages: distribute the inputs across the points at stage L
   Ilist = 1:2^L:N;
   for J=1:min([G,Pout])
      x(J-1 + Ilist) = x(Ilist);
   end % for J
   
   % set decision thresholds
   minGPlist = zeros(1,M);
   for I=1:M
      minGPlist(I) = min([Pout,2.^(I-1)]);
   end


   %--------------------------------------------------
   for L = Lfirstnoskip+1:M   % stages with no skipped groups and with up-branches
      G = 2^L;   % number in group
      G1 = G/2;
      U = 1 + 0i;
      W = cos(pi/G1) - sin(pi/G1)*i;
      numtoskip = G;
      if(Pout < G)  % do only upper branches
         for J=1:minGPlist(L)    %   min([G1,P]);
         for I = J:numtoskip:N
            IP = I + G1;
            T = x(IP) * U;   % input from lower branch
            numCmult = numCmult + 1;
            x(I) = x(I) + T;   % upper branch
            numCadd = numCadd + 1;
         end % for I
         U = U*W;
      end % for J
   else           % P out> G1 do both lower and upper branches of butterfly
      NG = floor(Pout/G);   % number of full groups
      NGleft = Pout - G1*NG;
      
      for J=1:minGPlist(L)    % min([G1,P]);   % do bother lower and upper branches
         for I = J:numtoskip:N
            IP = I + G1;
            T = x(IP) * U;   % input from lower branch
            numCmult = numCmult + 1;
            x(IP) = x(I) - T;  % lower branch
            x(I) = x(I) + T;   % upper branch
            numCadd = numCadd + 2;
         end % for I
         U = U*W;
      end % for J
   end % if P <= G1
end % for L
