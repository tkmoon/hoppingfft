
# hopped fft
import numpy as np
import math
from fftP import *

def makesubFFTsandInitial(x,N,P):
    # Make the subFFTs table and compute the initial FFT
    
    # global variables used by all these functions
    global subFFTs     # table of subFFTs
    global subFFTptr   # index into table of subFFTs
    global phases      # precompute the phases for speed
    
    M = N//P
    subFFTs = np.zeros((M,N),dtype=complex) # table of subFFTs
    subFFTptr = 0                           # index to subFFTs table    
    
    for i in range(M): # build the array of subFFTS
        #xi = subFFTs[i]    # get the ith row
        subFFTs[i,:P] = x[i*P:(i+1)*P]
        fftP(subFFTs[i],P,N)  # compute the FFT, in place (write into row i)
    # end for i
    phases = -1j*2*math.pi*np.array(range(N))*P/N
    X = np.zeros((N,), dtype=complex)
    
    # Build the first FFT outputs
    for i in range(M):
        X += np.exp(phases*i) * subFFTs[i]
    phases = np.exp(-phases)   # prepare for hopping updates
    return X


def FFThop(X,x,N,P):
    # Compute an update to the length-N FFT.  
    # X is the last FFT
    # X is modified IN PLACE, not returned
    # x consists of P points which are used in the update
    # N is the length of the FFT
    # P is the length of the hop (the new piece)

    global subFFTs
    global subFFTptr
    global phases


    M = N//P;
    X -= subFFTs[subFFTptr]        # downdate:remove the oldest
    subFFTs[subFFTptr,:P] = x[:P]  # copy the new data
    fftP(subFFTs[subFFTptr],P,N)   # create the transform in place of the old data
    X += subFFTs[subFFTptr]        # update: add in the newest
    subFFTptr = (subFFTptr + 1) % M
    X *= phases     

if __name__ == '__main__':  
    N = 128   # block size
    P = 32    # size of skip
    M = N // P

    extrablocks = 40;
    nblocks = M + extrablocks;
    x = (1 + 1j)*np.arange(1,nblocks*P + 1) # simple test data
    X = makesubFFTsandInitial(x[0:N],N,P)  # make the subDFTs and initial FFT
    Xffttrue = np.fft.fft(x[0:N])       # compute the FFT to compare:
    err = np.linalg.norm(X - Xffttrue);
    print('initial: error',err)
    

    # test the hopping for the extra steps
    for i in range(extrablocks):
        FFThop(X,x[N + i*P:N + (i+1)*P],N,P)
        Xffttrue = np.fft.fft(x[(i+1)*P: N + (i+1)*P])
        err = np.linalg.norm(X - Xffttrue)
        print('i=',i,'err=',err)
            

