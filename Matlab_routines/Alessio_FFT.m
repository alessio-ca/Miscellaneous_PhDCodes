function [f,spectrum]=Alessio_FFT(x,fs)
%This function computes the FFT of a real-valued series x, returning the
%spectrum corresponding to positive frequencies as well as the frequency
%vector assuming a sampling rate fs

N=2*floor(length(x)/2);
df=fs/N;
X=fft(x,N);

sampleIndex=-N/2 : N/2 - 1;

spectrum=fftshift(abs(X));
spectrum=spectrum(N/2+1:end);

f=df*sampleIndex(N/2+1:end);

end

