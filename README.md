README for github

This repository contains files related to H-FFT, computing fast
updates for a hopped FFT.

hoppedFFT.py:
	makesubFFTsandInitial -- makes the initial sub-FFTs and computes
	the initial FFT using them

	FFThop -- does the hopped FFT transform.

	Test functions: generates data, computes the initial transform,
	then hops for a number of steps, comparing the results with an FFT
	computed using the library fft function

fftP.py
	fftP --- function to compute the fft of a block of data that is
	mostly zeros


	fastlog2 -- computes log_2(x), where x is a power of 2

	Test code: generates data and computes fftP, then compares with
	the fft computed using the library fft function


fftP.m -- computes the fftP function (Matlab)

hopFFTinitial.m -- computes the initial subFFTs

makesubFFTsandInitial.m -- computes the initial subFFTS, and computes
						the initial FFT

FFThop.m  --- does the hop (after initialization) 

hoppedfft.m -- tests the H-FFT function.

maketable.m --- computes the long table of results

FFTtable.tex --- makes the table of results.  The table contents are
			 obtained by cut-and-paste from output of maketable.m


