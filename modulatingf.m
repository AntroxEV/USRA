function [ fi ] =modulatingf(t,Ts)
% Cacciola and Deodatis (2011) 
beta=9/Ts; 
t1=2.5/beta;
t2=11.5/beta;
% beta=0.6;
% t1=3;
% t2=23;
for hh=1:length(t)
     if t(hh)<=t1 
        fi(hh)=(t(hh)/t1)^2;
     end
        
    if t(hh)>=t1 && t(hh)<=t2
               
        fi(hh)=1;
    end
        
    if t(hh)>=t2 
        
        fi(hh)=exp(-beta*(t(hh)-t2));
    end
    
end  
    
end


