function [Matr,omegatot] = QSTHgeneration(nU,t,ww,G0,phi)
% SHINOZUKA GENERATION
polin=0;
nt=length(t);
dt=t(2)-t(1);
wN=1/(2*dt)*2*pi;
omegatot=linspace(-wN,wN,nt);
df=1/(nt*dt);
w=zeros(nU,nt);
phir=0+(2*pi-0)*rand(length(omegatot),nU);
G0=G0(:);
ww=ww(:);
S0=[flip(G0(2:end))/2; G0/2];
wwq=[-flip(ww(2:end)); ww];
S0i = interp1(wwq,S0,omegatot);
S0i(isnan(S0i))=0;
for jj=1:nU
    for hh=1:length(omegatot)   
    w(jj,:)=w(jj,:)+sqrt(1*2*S0i(hh)*2*pi*df)*(cos(omegatot(hh).*t+phir(hh,jj)));
    end
w(jj,:)=w(jj,:).*phi;
end
Matr=[t.' w.'];