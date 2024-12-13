function [pm,npeak,lb0]=peakmean(PSD,wc,p,Ts)
% p =0.5 is the not-exceeding probability
% wc is the vector of circular frequencies
% Ts is the time observing window (es 20 s)
lb0=trapz(wc,PSD);
lb1=trapz(wc,wc.*PSD);
lb2=trapz(wc,(wc.^2).*PSD);
niplus=sqrt(lb2/lb0)/2/pi;
q=sqrt(1-lb1^2/(lb0*lb2));
npeak=sqrt(2.*log(2.*niplus.*Ts./(-log(p)).*(1-exp(-q.^1.2.*sqrt(pi.*log(2.*niplus.*Ts./(-log(p))))))));
pm=sqrt(lb0).*npeak;