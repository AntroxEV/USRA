function [Afft,t2,H,omega,z]=LinConvTH(layers,signal,flag,G,ro,csi,acc,t)
% Site response analysis 1D on shear column. Equivalent linear
% viscous-elastic soil

%-------------M O D I F I C H E -------------------------------------------
shake91=0;                  %1  shake91 (not tested), 0  shake
load('solversettings.mat','sdata');
npoint=sdata.zeropad;

%% Initialization
z=[0; cumsum(layers)];
% hbed=sum(layers);                    
n=length(z);
if flag==0
    ilayer = n;
elseif flag==1
    ilayer=1;
else
    %error
end
% nlayer=length(layers);
% G=[G; G(end)]; %initial values
% csi=[csi; csi(end)]; %initial values
% ro=[ro; ro(end)]; %initial values
% mat=[materiale; materiale(end)]; %initial values
Gmax=G; %initial values
Csi=csi; %initial values
Ro=ro; %initial values
%% SOLVER TH     
dim=1;
dt=t(2)-t(1);
nt=length(t);
t2=t;
if npoint ~= 0
xspace=logspace(0,1,npoint)-1;
coda=acc(end)/(9)*(9-xspace);
acc=[acc;coda'];
t=[t;transpose(t(end)+dt:dt:(npoint+length(t)-1)*dt)];
clear coda xspace
end
fs=1/dt;
wb=20*2*pi;
[B, A]=butter(2,wb/2/fs,'low');
y=filter(B,A,acc);
y=baseline(y,dt,2);
[fomega,omega]=FFTcompl(y(:,4),t);
indw=find(omega==0);
fomega = fomega(indw:end).';  
omega = omega(indw:end).';    
numfreq=length(omega);      
Gx=zeros(n,1);
%------Matrice dei moduli di taglio complessi-------------------------
    switch shake91
            case 0
                for i=1:n
                    Gx(i)=G(i)*(1+2*(-1)^.5*csi(i));
                end
            case 1
                for i=1:n
                    Gx(i)=G(i)*((1-2*(csi(i))^2)+2*(-1)^.5*csi(i)*(1-(csi(i))^2)^.5);
                end            
    end    
    %------Matrice dei numeri d'onda complessi dei layer------------------
    for h=1:numfreq
        for i=1:n
            k=((Ro(i)*(omega(h))^2)/Gx(i))^0.5 ; 
            K(i,h)=k;                   %matrice dei numeri d'onda complessi
        end
    end
    %-------Imposizione condizione contorno delle amplificazioni A=B=1 ----
    amplA=1;
    amplB=1;
    %-------matrice delle amplificazioni A e B (condizioni di continuita)-
    for h=1:numfreq
        A(1,h)=1;
        B(1,h)=1;
        for i=2:n
            alpha=((Ro(i-1)*Gx(i-1))/(Ro(i)*Gx(i)))^0.5;
            A(i,h) = 0.5*A(i-1,h)*(1+alpha)*exp((-1)^.5*K(i-1,h)*layers(i-1))+ 0.5*B(i-1,h)*(1-alpha)*exp(-(-1)^.5*K(i-1,h)*layers(i-1));
            B(i,h) = 0.5*A(i-1,h)*(1-alpha)*exp((-1)^.5*K(i-1,h)*layers(i-1))+ 0.5*B(i-1,h)*(1+alpha)*exp(-(-1)^.5*K(i-1,h)*layers(i-1));
        end
    end
    %-----Funzione di trasferimento----------------------------------------
    for i=1:n
        H(i,:,dim)=(A(i,:)+B(i,:))./(A(ilayer,:)+B(ilayer,:));        %funzione di trasferimento
    end
    %----Segnale di input sismico-----------------------------------------
    signal=lower(signal);
    switch signal
        case 'outcrop'
                Ubed=(A(n,:)+B(n,:))./(2*A(n,:)).*fomega;
        case 'inside'
            Ubed=fomega;
    end
    %----Spostamento free field--------------------------------------------
    for i=1:n
        Uff(i,:,dim)=Ubed.*H(i,:,dim);
    end
%-----------------FINE CICLO-----------------------------
Afft=zeros(n,nt,dim);
for hh=1:n
    U=Uff(hh,:,dim);
    U=[fliplr(conj(U(2:end))) U];
    Ut=ifft(ifftshift(U),'symmetric');
    Afft(hh,1:nt,dim)=diffsec(Ut(1:nt),dt);
end
    