% Introduzione di R2 per effettuare deconvoluzione (evita rumore)
%% Solver Settings
R=0.65;                     % ration gamma eff. e gamma max
R2=1;                     %SEE KRAMER - (Silva 1988) - R2=0.87 for deconv
shake91=0;                  %1  shake91 (not tested), 0  shake
loops=80;                   %# cicli massimo; integer>0 o 'inf' fino a convergenza.
polin=3;                    %grado del polinomio per baseline correction
tolmax=0.0001;              %toleranza percentuale per convergenza
% segnale='inside';           %signal (outcrop or inside)
% ilayer=41;                   %#layer where signal is applied(only for inside)
% controllo=exist('npoint','var');
% if controllo==0
%     npoint=0;
% end
% clear controllo
npoint=500;
