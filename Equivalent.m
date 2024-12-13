function [rappGi, csii]=Equivalent(NLcurve,gammaeff,mat)
Xg=NLcurve(mat).Xg;
Yg=NLcurve(mat).Yg;
Xcsi=NLcurve(mat).Xcsi;
Ycsi=NLcurve(mat).Ycsi;
if gammaeff == 0
    rappGi=1;
    csii = Ycsi(1);
else
    rappGi=interp1(log10(Xg),Yg,log10(gammaeff),'linear','extrap');
    csii=interp1(log10(Xcsi),Ycsi,log10(gammaeff),'linear','extrap');
end