function [G,npeak]=SpCompPSD(RSA,wc,Ts)
%P. Cacciola, P. Colajanni, G. Muscolino, Combination of Modal Responses
%Consistent with Seismic Input Representation, Journal of Structural 
%Engineering (ASCE), 130(1),  pp. 47-55, 2004. (2004). 
nw=length(wc);
xi=0.05;
dw=wc(2)-wc(1);
p=0.5;
G=0;
for hh=2:nw
    Nu=Ts/2/pi*wc(hh)/(-log(p));
    du=(1-1/(1-xi^2)*(1-2/pi*atan(xi/sqrt(1-xi^2)))^2)^0.5;
    npeak(hh)=sqrt(2*log(2*Nu*(1-exp(-du^1.2*sqrt(pi*log(2*Nu))))));
    G(hh)=4*xi/(wc(hh)*pi-4*xi*wc(hh-1))*(RSA(hh)^2/npeak(hh)^2-dw*sum(G));
end
    G(G<0)=0;
    npeak(G<0)=0;