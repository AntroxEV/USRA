function [Afft,t2,H,omega,z]=LinEqConvTH(layers,signal,flag,G,ro,csi,ndata,mat,acc,t)
% Site response analysis 1D on shear column. Equivalent linear
% viscous-elastic soil

%-------------M O D I F I C H E -------------------------------------------
shake91=0;                  %1  shake91 (not tested), 0  shake
load('solversettings.mat','sdata');
R=sdata.R;                     % ration gamma eff. e gamma max
R2=sdata.R2;                     %SEE KRAMER - (Silva 1988) - R2=0.87 for deconv
shake91=0;                  %1  shake91 (not tested), 0  shake
loops=sdata.maxloops;                   %# cicli massimo; integer>0 o 'inf' fino a convergenza.
polin=sdata.polin;                    %grado del polinomio per baseline correction
tolmax=sdata.tolmax;              %toleranza percentuale per convergenza
npoint=sdata.zeropad;

%% Initialization
z=[0; cumsum(layers)];
% hbed=sum(layers);
zi=[layers./2;0]; %relative height for each layer (for shear deform)                       
n=length(z);
if flag==0
    ilayer = n;
elseif flag==1
    ilayer=1;
else
    %error
end

% nlayer=length(layers);
% G=[G(1); G]; %initial values
% csi=[csi(1); csi]; %initial values
% ro=[ro(1); ro]; %initial values
% mat=[materiale(1); materiale]; %initial values
Gmax=G; %initial values
Csi=csi; %initial values
Ro=ro; %initial values
Mat=mat; %initial values
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
%------Linear Equivalent model-----------------------------------------
ciclo=0;
tol=inf;
Gx=zeros(n,1);
while tol >= tolmax                %inizio dei cicli
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
        for h=1:numfreq
            for i=1:n
            Gamma(i,h)= (-1)^.5*K(i,h)*(A(i,h)*exp((-1)^.5*K(i,h)*zi(i))-B(i,h)*exp(-(-1)^.5*K(i,h)*zi(i)));
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
    %----Deformazioni free field--------------------------------------------
    for i=1:n
        Gammaff=Gamma(i,:).*(Uff(1,:,dim)./2);
        Gammaff=[fliplr(conj(Gammaff(2:end))) Gammaff];
        Gammat(i,:)=ifft(ifftshift(Gammaff),'symmetric');
    end 
    %---Calcolo nuovi parametri lineari equivalenti-------------------------
    for i=1:length(G)
        mat=Mat(i);
        if mat == 0
            Gnew(i,1)=Gmax(i);
            csii(i,1)=csi(i)*100;
        else
            gammaeff=(max(abs(Gammat(i,:,dim))).*R*R2)*100;      %Deformazione massima efficacie percentuale
            Gammaeffmax(i)=gammaeff;                            %vettore delle deformazioni efficaci per strato
            [rappGi, Csii]=Equivalent(ndata,gammaeff,mat);             %funzione per il calcolo dei parametri equivalenti
            Gnew(i,1)=rappGi.*Gmax(i);
            csii(i,1)=Csii;
        end
    end
    warning off 'MATLAB:divideByZero'
    Cconv=(abs(csi-csii./100))./csi.*100;     %differenza perc. tra due iterazioni successive di csi
    Gconv=(abs(G-Gnew))./G.*100;          %differenza perc. tra due iterazioni successive di G 
    warning on 'MATLAB:divideByZero'
    conver=[Cconv;Gconv];
    tol=max(conver(:));
    G=Gnew;
    csi=csii./100;
    if ciclo > loops
        warning('Convergence not reached, change parameters')
        break
    end
    ciclo=ciclo+1;
end     %fine ciclo
%-----------------FINE CICLO-----------------------------
Afft=zeros(n,nt,dim);
for hh=1:n
    U=Uff(hh,:,dim);
    U=[fliplr(conj(U(2:end))) U];
    Ut=ifft(ifftshift(U),'symmetric');
    Afft(hh,1:nt,dim)=diffsec(Ut(1:nt),dt);
end
    