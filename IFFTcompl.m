function [Ut]=IFFTcompl(U)
Ut=ifft(ifftshift(U),'symmetric');
