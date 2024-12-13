function [PSD,freq]=acc2PSD(Matr)
% Determinatiof of the average Power Spectral Density from realizations
%-----------------------------------------------------------
% Matr is the matrix of the realizations
% row= time 
%column = number of accelerograms
%first column is the time
%
[N m]=size(Matr);
nA=m-1;
t=Matr(:,1);
Fs=1/(1*(t(2)-t(1)));
if mod(N,2)>0
    Matr=[Matr; zeros(1,m)];
    Matr(end,1)=t(end,1)+(t(2)-t(1));
    [N m]=size(Matr);
end
for index=2:nA+1
    acc=Matr(:,index);
    xdft = fft(acc);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)).*abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/N:Fs/2;
%     PSDi(:,index)=psdx./(sqrt(2)*pi);
    PSDi(:,index)=psdx./(2*pi);
end
PSD=mean(PSDi(:,2:end),2);