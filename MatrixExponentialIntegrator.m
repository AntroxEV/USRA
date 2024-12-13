function [Y]=MatrixExponentialIntegrator(K,M,C,f,dt)
%% Code adapted from Laura's code 30/04/2014
%% MODAL ANALYSIS
n=size(K,1);
[PHI,L]=eig(K,M);
omega=sqrt(diag(L));%------------------------------natural frequencies [rad/sec]
norm_=PHI'*M*PHI;
norm=diag(norm_);
for i=1:n
PHI_norm(:,i)=PHI(:,i)/sqrt(norm(i));
end
% O_norm1=PHI_norm'*M*PHI_norm; %-------------------------------must check identity matrix I
% O_norm2=PHI_norm'*K*PHI_norm; %-------------------------------must check
OM=diag(omega);
CC=PHI_norm'*C*PHI_norm;
%% INPUT - SEISMIC FORCE
n_step=length(f); %-------------------------------------------------time history number of intervals
% dt=0.01;%------------------------------------------duration of each interval
duration=n_step*dt;
time=linspace(0,duration,n_step);
data = f;
%% STATE SPACE MODAL ANALYSIS 
D_M=[zeros(n) eye(n); -OM^2 -CC];
teta=expm(D_M*dt);
L_M=(teta-eye(2*n))/(D_M);
gamma0=(teta-L_M/dt)/(D_M);
gamma1=(L_M/dt- eye(2*n))/(D_M);
PI=[PHI_norm zeros(n) ; zeros(n) PHI_norm ];
% tao=ones(n,1)*(-1);
% p=PHI_norm'*M*tao;
% p=PHI_norm'*eye(n,n)*tao;
p=PHI_norm'*eye(n,n);
Vm=[zeros(n,n);p];
% participation factor check %
Mtot_check=p'*p;

%% STEP BY STEP INTEGRATION
Z(:,1)=zeros(2*n,1);
Y(:,1)=PI*Z(:,1);
  
for i1=1:n_step-1
    Z(:,i1+1)=teta*Z(:,i1)+gamma0*Vm*data(i1)+gamma1*Vm*data(i1+1);
    Y(:,i1+1)=PI*Z(:,i1+1);
end