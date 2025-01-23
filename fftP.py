# fftP.py
import numpy as np
import math

def fftP(x,P,Pout):
    # Compute the fft of the complex input vector x
    #    using only P inputs and P outputs
    # The length of x is must be a power of 2.
    # P is a power of 2, P <= N/2.
    # The FFT is computed IN PLACE, so x is modified
    # Returns: x, numCmult, numCadd
    n = len(x);  n_2 = n//2;   m = fastlog2(n);
    # Bit reverse shuffle (following Rabiner/Gold) --
    j = 0
    for i in range(n-1):
        if(i < j):  (x[j], x[i]) = (x[i], x[j]) # swap
        k = n_2
        while(1):
            if(k <= j):  j -= k; k //= 2; continue
            break
        # end while
        j += k
    # end for i
    #-------------------------------------------------
    Lfirstnoskip = m - fastlog2(P)
         # first L that does not skip blocks
    numCmult = 0; numCadd = 0;
                       # number of multiplies and adds
    ell = Lfirstnoskip; G = 1<<ell
    # Initial Stage: Distribute the data x[0:P-1]
    #   (bit reverse indexed) across the points in x
    ilist = np.array(range(0,n,G))
    for j in range(min(G,Pout)): x[j + ilist] = x[ilist]

    minGPlist = np.zeros(m+1,dtype=int)
    for i in range(1,m+1):
        minGPlist[i] = min(Pout,  1<<(i-1) )
    #-------------------------------------------------
    for ell in range(Lfirstnoskip+1,m+1): # rest the stages
        G1 = G; G <<= 1  # G1 = 2^(L-1);  G = 2^L
        U = complex(1,0)   # twiddle factor
        W = complex(math.cos(math.pi/G1),
            -math.sin(math.pi/G1)) # twiddle factor scale
        if(Pout < G):  # do only upper branches of butterflies
            for j in range(minGPlist[ell]):
                # i.e.: for j in range(min(Pout,(1<<(ell-1))))
                # j counts butterflies within a group
                for i in range(j,n,G):
                    # i loops over the groups
                    ip = i + G1
                    t = x[ip] * U # take input from lower branch
                    numCmult += 1; numCadd += 1
                    x[i] += t
                # end for i
                U *= W
            # end for j
        else: # Do both lower and upper branches of butterflies
            for j in range(minGPlist[ell]):
                # i.e., for j in range(min(Pout,(1<<(ell-1))))
                for i in range(j,n,G):
                    ip = i + G1
                    t = x[ip] * U # input from lower branch
                    x[ip] = x[i] - t # lower branch
                    x[i] = x[i] + t # upper branch
                    numCmult += 1;  numCadd += 2
                # end for i
                U *= W
            # end for j
        # end if(Pout<G) / else
    # end for ell
    return x,numCmult,numCadd
# end def fftP

def fastlog2(x):
    # compute log2(x), where x is a power of 2
    return x.bit_length() - 1

def fftPnumCmult(N,Pout):
    M = fastlog2(N)
    Lfirstnoskip = fastlog2(N//P) # first L that does not skip blocks
    L1 = Lfirstnoskip + 1  # first block with multiplies
    Lequal = 1 + fastlog2(Pout) # first L where Pout = 2^(L-1)
    numCmult = 0
    if(Lequal >= L1):
        numCmult = (min(Lequal,M) - L1 + 1)*N//2
        if(M >= Lequal + 1):
            numCmult += Pout*((1 << (M - Lequal)) - 1)
    else:
        numCmult += numCmult + Pout*( (1 <<(M - L1 + 1)) - 1)
    return numCmult

#------------------ test code ------------------------
if __name__ == '__main__':  
    # test case
    N = (1<<19);  m = fastlog2(N)
    Plist = 1<< np.array(range(1,m))
    for P in Plist:
        Pout = int(P)
        P = int(P);  x = np.zeros(N,dtype=complex)
        x[range(P)] = (np.arange(P)+1)+1j*(np.arange(P)+P+1)
        xfft_true = np.fft.fft(x)   # used for testing
        x,numCmult, numCadd = fftP(x,P,Pout)
        print('N=',N,'P=',P,'Pout=',Pout,'numMult=',numCmult,
            'numAdd=',numCadd,'ratio=',(N/2)*math.log2(N)/numCmult,
        'error=',np.linalg.norm(x[0:Pout-1]-xfft_true[0:Pout-1]))
    # end for P
