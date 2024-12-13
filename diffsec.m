function [der]=diffsec(funz,dh)
%Second-order Derivative
%Differenziale secondo. Metodo alla differenze finite
dim=length(funz);
for h=2:dim-1
    der(h)=(funz(h+1)-2*funz(h)+funz(h-1))/dh^2;
end
der=[der 0];