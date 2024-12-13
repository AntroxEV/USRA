function [Hsoil,PSD]=LinEqConvPSD(spessorilay,passo,G,ro,csi,materiale,PSDbed,omega)
% LINEAR EQUIVALENT SITE RESPONSE ANALYSIS OF A HOMOGENEOUR SOIL DEPOSIT 
% BY MEANS OF LINEAR CONVOLUTION
% TO OBTAIN SOIL DISPLACEMENTS AND TF
%Input:
%...v_dz VECTOR OF LAYERS THICKNESS (must be constant)
%...G,ro,csi VECTORS OF INITIAL LINEAR ELASTIC SOIL PROPERTIES
%...mat Type of soil
%.. PSDbed,omega Power spectral density and circular frequencies
% -------------------------------------------------------------------------

%% Solver Settings
R=0.65;                     % ration gamma eff. e gamma max
% R2=1;                     %SEE KRAMER - (Silva 1988) - R2=0.87 for deconv
shake91=0;                  %1  shake91 (not tested), 0  shake
loops=80;                   %# cicli massimo; integer>0 o 'inf' fino a convergenza.
% polin=0;                    %grado del polinomio per baseline correction
tolmax=0.0001;              %tolleranza percentuale per convergenza
segnale='inside';           %signal (outcrop or inside)
controllo=exist('npoint','var');
if controllo==0
    npoint=0;
end
clear controllo
%% Initialization
times=spessorilay./passo;                  %# discretizzazioni layer(sublayer)
z=zeros(sum(times)+1,1);                         
n=size(spessorilay,2);                     %# strati del deposito (layer)
ind=1;
for i=1:size(spessorilay,2)
    for h=1:times(i)
    z(ind+1)=z(ind)+passo(i);           % vettore posizione top sublayer
    ind=ind+1;
    end
end
spessori=diff(z);                       %spessori di ogni sublayer
z3=[spessori./2;0];                         %quota di calcolo deformazioni
discr=size(z,1);
ilayer=discr;                   %#layer where signal is applied(only for inside)
if strcmp(segnale,'outcrop')==1;
    ilayer=discr;
end
segnale=lower(segnale);
    if strcmp(segnale,'inside')
    G=[G; G(end)];
    ro=[ro; ro(end)];
    csi=[csi; csi(end)];
    materiale=[materiale; materiale(end)];
    end
Gmax=G;
Csi=csi;
Ro=ro;
Mat=materiale;

%% SOLVER - EQUIVALENT LINEAR MODEL
numfreq=length(omega);              %# frequenze
ciclo=0;
tol=inf;
while tol >= tolmax                %inizio dei cicli
    %------Matrice dei moduli di taglio complessi-------------------------
    switch shake91
            case 0
                for i=1:discr
                    Gx(i)=G(i)*(1+2*(-1)^.5*csi(i));
                end
            case 1
                for i=1:discr
                    Gx(i)=G(i)*((1-2*(csi(i))^2)+2*(-1)^.5*csi(i)*(1-(csi(i))^2)^.5);
                end            
    end    
    %------Matrice dei numeri d'onda complessi dei layer------------------
    for h=1:numfreq
        for i=1:discr
            k=((Ro(i)*(omega(h))^2)/Gx(i))^0.5 ; 
            K(i,h)=k;                   %matrice dei numeri d'onda complessi
        end
    end
    %-------Imposizione condizione contorno delle amplificazioni A=B=1 ----
    amplA=1;
    amplB=1;
    %-------matrice delle amplificazioni A e B (condizioni di continuit�)-
    for h=1:numfreq
        A(1,h)=1;
        B(1,h)=1;
        for i=2:discr
            alpha=((Ro(i-1)*Gx(i-1))/(Ro(i)*Gx(i)))^0.5;
            A(i,h) = 0.5*A(i-1,h)*(1+alpha)*exp((-1)^.5*K(i-1,h)*spessori(i-1))+ 0.5*B(i-1,h)*(1-alpha)*exp(-(-1)^.5*K(i-1,h)*spessori(i-1));
            B(i,h) = 0.5*A(i-1,h)*(1-alpha)*exp((-1)^.5*K(i-1,h)*spessori(i-1))+ 0.5*B(i-1,h)*(1+alpha)*exp(-(-1)^.5*K(i-1,h)*spessori(i-1));
        end
    end
        for h=1:numfreq
            for i=1:discr
            Gamma(i,h)= (-1)^.5*K(i,h)*(A(i,h)*exp((-1)^.5*K(i,h)*z3(i))-B(i,h)*exp(-(-1)^.5*K(i,h)*z3(i)));
            end
        end
    %-----Funzione di trasferimento----------------------------------------
    for i=1:discr
        H(i,:)=(A(i,:)+B(i,:))./(A(ilayer,:)+B(ilayer,:));        %funzione di trasferimento
    end
%----Segnale di input sismico-----------------------------------------
    segnale=lower(segnale);
    switch segnale
        case 'outcrop'
            for h=1:size(omega,1);
%                 Ubed(1,h)=(A(discr,h)+B(discr,h))/(2*A(discr,h))*fomega(h);
                diplay('TO DO - ERROR IN OUTCROP TRANSFER FUNCTION')
            end
        case 'inside'
            Ubed=PSDbed;
    end
%     Ubed(:,:)=1
    %----Spostamento free field--------------------------------------------
    for i=1:discr
        PSD(i,:)=abs(H(i,:)).^2.*Ubed.';
    end
    Hsoil=H;
    %----Deformazioni free field--------------------------------------------
    for i=1:discr
        Gammaff(i,:)=abs(Gamma(i,:)).^2.*(PSD(i,:)./2);
    end
    % ----Calcolo del picco massimo - valore p=0.5------------------------
        for i=1:length(G)
            Gammamax(i)=peakmean(Gammaff(i,:),omega.',0.5,20);
        end
    %---Calcolo nuovi parametri lineari equivalenti-------------------------
    for i=1:length(G)
        mat=Mat(i);
        if mat == 0
            Gnew(i,1)=Gmax(i);
            csii(i,1)=csi(i)*100;
        else
            gammaeff=(Gammamax(i).*R)*100;      %Deformazione massima efficace percentuale
            Gammaeffmax(i)=gammaeff;                            %vettore delle deformazioni efficaci per strato
            [rappGi Csii]=Equivalent(gammaeff,mat);             %funzione per il calcolo dei parametri equivalenti
            Gnew(i,1)=rappGi.*Gmax(i);
            csii(i,1)=Csii;
        end
    end
    warning off 'MATLAB:divideByZero'
    Cconv=(abs(csi-csii./100))./csi.*100;     %differenza perc. tra due iterazioni successive di csi
    Gconv=(abs(G-Gnew))./G.*100;          %differenza perc. tra due iterazioni successive di G 
    warning on 'MATLAB:divideByZero'
    conver=[Cconv';Gconv'];
    tol=max(conver);
    G=Gnew;
    csi=csii./100;
    if ciclo > loops
        warning('Convergenza non raggiunta')
        break
    end
    ciclo=ciclo+1;
end     %fine ciclo


    