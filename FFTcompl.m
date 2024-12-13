function [fomega,omega]=FFTcompl(acc,t,varargin)
% FFT of the accelerograms in order to obtain accelerogram  in 
% frequency domain
dt=t(2)-t(1);
Fs=1/dt;
if nargin==2
    NFFT= 2^nextpow2(length(acc));
else
    NFFT=varargin{1};
end
if mod(NFFT,2)==0
   freq = Fs/2*linspace(-1,1,NFFT+1);
   freq(end) = [ ];
else
    freq = Fs/2*linspace(-1,1,NFFT);
end
y=fft(acc,NFFT);
fomega = fftshift(y);

omega=2*pi*freq';

% -------------------------------------------------------------------------
