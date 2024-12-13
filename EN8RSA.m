function [RSA,T]=EN8RSA(ag,soiltype,varargin)
eta=1;
if nargin ==2
    T1=0.02;
    T2=4;
    T=linspace(T1,T2,100);
else
    T=varargin{1};
end
switch soiltype
    case 'A'
        S=1;
        TB=0.15;
        TC=0.4;
        TD=2;
    case 'B'
        S=1.2;
        TB=0.15;
        TC=0.5;
        TD=2;
    case 'C'
        S=1.15;
        TB=0.2;
        TC=0.6;
        TD=2;
    case 'D'
        S=1.35;
        TB=0.2;
        TC=0.8;
        TD=2;  
    case 'E'
        S=1.4;
        TB=0.15;
        TC=0.5;
        TD=2;
end
        
nT=length(T);
for hh=1:nT
    if T(hh)<=TB
        RSA(hh)=ag*S*(1+T(hh)/TB*(eta*2.5-1));
    elseif TB<T(hh) && T(hh)<=TC
        RSA(hh)=ag*S*eta*2.5;
    elseif TC<T(hh) && T(hh)<=TD
        RSA(hh)=ag*S*eta*2.5*TC/T(hh);
    elseif TD<T(hh)&& T(hh)<=4
        RSA(hh)=ag*S*eta*2.5*TC*TD/T(hh)^2;
    else
        RSA(hh)=0;
    end
end
% plot(T,RSA);


