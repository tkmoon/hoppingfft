% Make tables with results for hopfft

% TKM Jan 2025

Mlist = 5:20;
Nlist = 2.^Mlist;

str1 = '';   % first row
str2 = '';
str3 = '';
str4 = '';

table1 = {};
table2 = {};

Pmaxtable1 = 2^10;
Nmaxtable1 = 2^11;
for N = Nlist
   M = log2(N);
   str1 = sprintf('$2^{%d}$: & \\#mult ',M);
   str2 = sprintf('          & ratio ');
   str3 = sprintf('          & ratio2 ');
   str4 = sprintf('          & ratio3 ');
   if(N > Nmaxtable1)  % start the second table as well
      str5 = sprintf('$2^{%d}$: & \\#mult ',M);
      str6 = sprintf('          & ratio ');
      str7 = sprintf('          & ratio2 ');
      str8 = sprintf('          & ratio3 ');
   end

   Plist = 2.^(1:(M-1));
   for P = Plist
      Pout = N/2;

      x = zeros(1,N);
      xsave = x;

      x(1:P) = 1:P;  % complex(1:P, P+1:2*P);
      doprint = 0;

      xwithrand = x;
      xwithrand(P+1:N) = randn(1,N-P);

      [XFFTfull,numCmult,numCadd] = FFTfull(x,doprint);
      % FFT with no no noise added
      Xfft = fft(x);  % FFT with no noise added

      xsavewithrand = xwithrand;
      x = xwithrand;   % set for the the reduced computation 

      [X,numCmult,numCadd] = fftP(x,P,Pout);

      updatecount = numCmult + N;
      ratio = N/2*log2(N)/numCmult;
      ratio2 = N/2*log2(N)/updatecount;
      ratio3 = P*N/updatecount;
      
      if(P <= Pmaxtable1)
         str1 = [str1 sprintf('&  %d',numCmult)];
         str2 = [str2 sprintf('& %.1f',ratio)];
         str3 = [str3 sprintf('& %.1f',ratio2)];
         str4 = [str4 sprintf('& %.1f',ratio3)];
      else
         str5 = [str5 sprintf('&  %d',numCmult)];
         str6 = [str6 sprintf('& %.1f',ratio)];
         str7 = [str7 sprintf('& %.1f',ratio2)];
         str8 = [str8 sprintf('& %.1f',ratio3)];
      end
      
      if(P == Plist(end)) % end of row
         if(P <= Pmaxtable1)  % save to first table
            str1 = [str1 '\\'];
            str2 = [str2 '\\'];
            str3 = [str3 '\\'];
            str4 = [str4 '\\[.3em]'];

            tl = length(table1) + 1;
            table1{tl} = str1;
            table1{tl+1} = str2;
            table1{tl+2} = str3;
            table1{tl+3} = str4;
         end
         if(P > Pmaxtable1) % save to second table
            str5 = [str5 '\\'];
            str6 = [str6 '\\'];
            str7 = [str7 '\\'];
            str8 = [str8 '\\[.3em]'];
            tl = length(table2) + 1;
            table2{tl} = str5;
            table2{tl+1} = str6;
            table2{tl+2} = str7;
            table2{tl+3} = str8;
         end
      end
   end % for P
   fprintf('Table 1\n');
   for i=1:length(table1)
      fprintf('%s\n',table1{i});
   end
   if(length(table2))
      fprintf('Table 2\n');
      for i=1:length(table2)
         fprintf('%s\n',table2{i});
      end
   end
end % for N



