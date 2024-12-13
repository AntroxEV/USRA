function [Hsoil,PSD,z]=LinConvPSD(layers,signal,flag,G,ro,csi,PSDbed,omega)
% LINEAR SITE RESPONSE ANALYSIS BY MEANS OF LINEAR CONVOLUTION
% TO OBTAIN SOIL Displacements
% INSIDE INPUT
%Input:
%...G,ro,csi VECTORS OF LINEAR ELASTIC SOIL PROPERTIES
% -------------------------------------------------------------------------
% Settings
shake91=0;                  %1  shake91 (not tested), 0  shake
%% Initialization
z=[0; cumsum(layers)];
hbed=sum(layers);                    
n=length(z);
if flag==0
    ilayer = n;
elseif flag==1
    ilayer=1;
else
    %error
end
nlayer=length(layers);
Gmax=G; %initial values
Csi=csi; %initial values
Ro=ro; %initial values
% -------------------------------------------------------------------------
% Soil deposit discretization
numfreq=length(omega);              %# frequenze
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
    %-------matrice delle amplificazioni A e B (condizioni di continuit�)-
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
        H(i,:)=(A(i,:)+B(i,:))./(A(ilayer,:)+B(ilayer,:));        %funzione di trasferimento
    end
%----Segnale di input sismico-----------------------------------------
    signal=lower(signal);
    switch signal
        case 'outcrop'
                Hbed=(A(nlayer,:)+B(nlayer,:))./(2*A(nlayer,:));
                Ubed=abs(Hbed).^2.*PSDbed.';
        case 'inside'
            Ubed=PSDbed.';
    end
    %----Spostamento free field--------------------------------------------
    for i=1:n
        PSD(i,:)=abs(H(i,:)).^2.*Ubed;
    end
    Hsoil=H;
    
    