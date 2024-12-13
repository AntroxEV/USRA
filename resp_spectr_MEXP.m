function [Sa,T] = resp_spectr_MEXP(acc,dt,varargin)
%% Data
T1=0.02;
T2=4;
csi=0.05;
if nargin==3
    n=varargin{1};
    T=linspace(T1,T2,n);
elseif nargin==4
        n=varargin{1};
        csi=varargin{2};
        T=linspace(T1,T2,n);
else
    step=0.02;
    T=T1:step:T2;
end
%
f=1./T;
nf=length(T);
M=1;
for hh=1:nf
    wo=2*pi*f(hh);
    K=wo^2*M;
    C=2*csi*wo*M;
    [Y]=MatrixExponentialIntegrator(K,M,C,acc,dt);
    Sa(hh)=wo^2*max(abs(Y(1,:)));
end