function [PSD,wg,zg,wf,zf] = PSDstatCP(w,typesoil,Sw,varargin)
%Power Spectral Density - Clough and Penzien
% PSD = PSDstatCP(w,typesoil)
% PSD = PSDstatCP(w,wg,zf)
%type soil: 1 firm, 2: medium, 3 soft :::: after De Kiureghian and Neuenhofer
%type soil: 4-8 EC8 compatible :::after Giaralis and Spanos
%After Yeh and When 1990:
%Event wg zg wf zf Sw
% EL CENTRO 
%TAFT 20 0.65 1 0.5 0.0028
%MEXICO CITY 4.2 0.1 0.23 0.1 0.0033
switch typesoil
    case 1
%% Firm soil
wg=15.0;
zg=0.6;
zf=0.6;
wf=0.1*wg;
    case 2
%% Medium soil
wg=10.0;
zg=0.4;
zf=0.6;
wf=0.1*wg;
    case 3
%% Soft soil
wg=5.0;
zg=0.2;
zf=0.6;
wf=0.1*wg;
    case 4
%% EC8 - Giaralis Spanos (Soil A)
wg=17.57;
zg=0.54;
zf=0.78;
wf=2.22;
    case 5
%% EC8 - Giaralis Spanos (Soil B)
wg=10.73;
zg=0.78;
zf=0.90;
wf=2.33;
    case 6
%% EC8 - Giaralis Spanos (Soil C)
wg=7.49;
zg=0.84;
zf=1.15;
wf=2.14;
    case 7
%% EC8 - Giaralis Spanos (Soil D)
wg=5.34;
zg=0.88;
zf=1.17;
wf=2.12;
    case 8
%% EC8 - Giaralis Spanos (Soil E)
wg=10.76;
zg=0.77;
zf=1.07;
wf=2.03;
end
if nargin==7
    wg=varargin{1};
    zg=varargin{2};
    wf=varargin{3};
    zf=varargin{4};
end

S0 = Sw^2/(pi*wg*(2*zg+1/(2*zg)));
% S = A^2*S0*((1+4*zg^2*(w/wg)^2)/((1-(w/wg)^2)^2+4*zg^2*(w/wg)^2))*((w/wf)^4/((1-(w/wf)^2)^2+4*zf^2*(w/wf)^2));
PSD=((S0*((1+4*zg^2*(w/wg).^2)))./((1-w.^2/wg^2).^2+4*zg^2*w.^2/wg^2)).*(((w/wf).^4)./((1-(w/wf).^2).^2+4*zf^2*(w/wf).^2));