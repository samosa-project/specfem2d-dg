% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [] = compare_models(atm_file_1,atm_file_2)
  if(not(exist(atm_file_1)==2))
    error(['[',mfilename,', ERROR] File ''',atm_file_1,''' does not exist.']);
  end
  if(not(exist(atm_file_2)==2))
    error(['[',mfilename,', ERROR] File ''',atm_file_2,''' does not exist.']);
  end
  
  set(0, 'DefaultTextInterpreter', 'latex'); set(0, 'DefaultLegendInterpreter', 'latex');
  addpath('/usr/local/matlab/r2018a/toolbox/tightfig'); % tightfig
  addpath('/home/l.martire/Documents/work/mars/mars_is'); % prettyAxes
  set(0, 'defaultTextFontSize', 10);
  set(0, 'defaultAxesFontSize', 10);
  
  dn1 = '1';
  dn2 = '2';
  
  [Z1, RHO1, T1, C1, P1, H1, G1, NBVSQ1, KAP1, MU1, MUvol1, Wnorth1, Weast1, W1, Cp1, Cv1, GAM1, FR1, SVIB1] = extract_atmos_model(atm_file_1, 3, 0, -1);
  [Z2, RHO2, T2, C2, P2, H2, G2, NBVSQ2, KAP2, MU2, MUvol2, Wnorth2, Weast2, W2, Cp2, Cv2, GAM2, FR2, SVIB2] = extract_atmos_model(atm_file_2, 3, 0, -1);
  
  f=figure('units','normalized','outerposition',[0 0 1 1]);
  
  nl=3;nc=6;i=1;
  Q1=RHO1; Q2=RHO2; subplot(nl,nc,i); semilogx(Q1, Z1, 'displayname', dn1); hold on; semilogx(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$\rho$ (max. diff. ',diffStr(Q1,Q2),' \%)']); legend('location', 'best');
  Q1=T1; Q2=T2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$T$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=C1; Q2=C2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$c$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=P1; Q2=P2; subplot(nl,nc,i); semilogx(Q1, Z1, 'displayname', dn1); hold on; semilogx(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$p$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=H1; Q2=H2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$H$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=G1; Q2=G2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$g$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=NBVSQ1; Q2=NBVSQ2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$N^2$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=KAP1; Q2=KAP2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$\kappa$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=MU1; Q2=MU2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$\mu$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=MUvol1; Q2=MUvol2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$\mu_{vol}$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=Wnorth1; Q2=Wnorth2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$W_N$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=Weast1; Q2=Weast2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$W_E$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=W1; Q2=W2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$W_P$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=Cp1; Q2=Cp2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$c_p$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=Cv1; Q2=Cv2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$c_v$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=GAM1; Q2=GAM2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$\gamma$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=FR1; Q2=FR2; subplot(nl,nc,i); semilogx(Q1, Z1, 'displayname', dn1); hold on; semilogx(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$f_r$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  Q1=SVIB1; Q2=SVIB2; subplot(nl,nc,i); plot(Q1, Z1, 'displayname', dn1); hold on; plot(Q2, Z2, 'displayname', dn2); hold on; i=i+1; xlabel(['$S_{vib}$ (max. diff. ',diffStr(Q1,Q2),' \%)']);
  
  tightfig;
  
  prettyAxes(f);
end

function str=diffStr(Q1, Q2)
  str = num2str(max(diffFormula(Q1, Q2))*100);
end

function d = diffFormula(Q1, Q2)
  d = abs( (Q1-Q2) ./ (0.5*(Q1+Q2)));
end
