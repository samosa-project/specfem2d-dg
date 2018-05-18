% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Estimate new run time based on previous runs' saved data.
% Last modified: See file metadata.
% Usage:         Fill section "Parameters" according to the calculation to
%                be performed.
% Notes:         The method (N-D nearest point search) is approximate and
%                does not take into account the fact that some parameters
%                have more impact than others. All in all, it is a very
%                rough and very approximate method.

% clear all
% close all
clc
format longG;

[data, t]=load(); % Load data (see function below).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nstations   = 114;
nstepseismo = 50;

neltot      = 748360;
neldg       = 748360;

nstepsnap   = 250;
nsteptot    = 8000;

nproc       = 32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of estimate.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
point=[nstations,neldg/neltot,nstepsnap/nsteptot,nstepseismo/nsteptot,neltot/nproc]; % Format point as data format. Currently [% elements as DG, % timesteps as snapshots, elements per proc].

disp(strcat("Current point:            [ ", sprintf("%.3e ", point), "]."));
idp=dsearchn(data,point);
disp(strcat("Closest data point found: [ ", sprintf("%.3e ", data(idp,:)), "]."));
time=t(idp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(" ");

hms=fix(mod(time*neltot*nsteptot,[0,3600,60])./[3600,60,1]); cputimestr=string([num2str(hms(1)),'h ',num2str(hms(2)),'m ',num2str(hms(3)),'s']);
hms=fix(mod(time*neltot*nsteptot/nproc,[0,3600,60])./[3600,60,1]); realtimestr=string([num2str(hms(1)),'h ',num2str(hms(2)),'m ',num2str(hms(3)),'s']);

disp(strcat("Expected time per element, per iteration: ",sprintf("%.3e", time), " s."));
disp(strcat("Expected run time:                        ",cputimestr, " (CPU)."));
disp(strcat("                                          ",realtimestr, " (real)."));

disp(" ");

disp("[WARNING] Recall the method used for estimation is very rough and approximate. Do not take the estimation for granted.");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function containing data.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d, t]=load()
  d=[
-1.0000000000000000E+000,1.0000000000000000E+000,2.5000000000000000E-002,-1.0000000000000000E+000,625;
-1.0000000000000000E+000,1.0000000000000000E+000,2.5000000000000000E-002,-1.0000000000000000E+000,156;
-1.0000000000000000E+000,1.0000000000000000E+000,2.5000000000000000E-002,-1.0000000000000000E+000,10;
-1.0000000000000000E+000,1.0000000000000000E+000,1.2500000000000000E-002,-1.0000000000000000E+000,2500;
-1.0000000000000000E+000,1.0000000000000000E+000,1.2500000000000000E-002,-1.0000000000000000E+000,625;
-1.0000000000000000E+000,1.0000000000000000E+000,1.2500000000000000E-002,-1.0000000000000000E+000,39;
-1.0000000000000000E+000,1.0000000000000000E+000,6.2500000000000000E-003,-1.0000000000000000E+000,10000;
-1.0000000000000000E+000,1.0000000000000000E+000,6.2500000000000000E-003,-1.0000000000000000E+000,2500;
-1.0000000000000000E+000,1.0000000000000000E+000,6.2500000000000000E-003,-1.0000000000000000E+000,156;
-1.0000000000000000E+000,1.0000000000000000E+000,6.2500000000000000E-003,-1.0000000000000000E+000,78;
-1.0000000000000000E+000,1.0000000000000000E+000,6.2500000000000000E-003,-1.0000000000000000E+000,39;
-1.0000000000000000E+000,1.0000000000000000E+000,3.1250000000000000E-003,-1.0000000000000000E+000,625;
-1.0000000000000000E+000,1.0000000000000000E+000,3.1250000000000000E-003,-1.0000000000000000E+000,313;
-1.0000000000000000E+000,1.0000000000000000E+000,3.1250000000000000E-003,-1.0000000000000000E+000,156;
-1.0000000000000000E+000,1.0000000000000000E+000,1.5625000000000000E-003,-1.0000000000000000E+000,2500;
-1.0000000000000000E+000,1.0000000000000000E+000,1.5625000000000000E-003,-1.0000000000000000E+000,1250;
-1.0000000000000000E+000,1.0000000000000000E+000,1.5625000000000000E-003,-1.0000000000000000E+000,625;
-1.0000000000000000E+000,1.0000000000000000E+000,0.0000000000000000E+000,-1.0000000000000000E+000,10000;
-1.0000000000000000E+000,1.0000000000000000E+000,0.0000000000000000E+000,-1.0000000000000000E+000,5000;
-1.0000000000000000E+000,1.0000000000000000E+000,0.0000000000000000E+000,-1.0000000000000000E+000,2500;
-1.0000000000000000E+000,9.0773809523809500E-001,1.6666666666666700E-003,-1.0000000000000000E+000,7875;
-1.0000000000000000E+000,9.0773809523809500E-001,1.6666666666666700E-003,-1.0000000000000000E+000,7875;
-1.0000000000000000E+000,6.6079295154185000E-001,1.6666666666666700E-002,-1.0000000000000000E+000,887;
-1.0000000000000000E+000,6.6079295154185000E-001,1.6666666666666700E-002,-1.0000000000000000E+000,887;
-1.0000000000000000E+000,9.7925311203319500E-001,7.4074074074074100E-003,-1.0000000000000000E+000,5648;
-1.0000000000000000E+000,1.0000000000000000E+000,6.8027210884353700E-003,-1.0000000000000000E+000,816;
-1.0000000000000000E+000,9.7925311203319500E-001,7.4719800747198000E-003,-1.0000000000000000E+000,5648;
-1.0000000000000000E+000,9.8090388287714800E-001,2.1413276231263400E-003,-1.0000000000000000E+000,5611;
2.0000000000000000E+000,6.6079295154185000E-001,3.1746031746031700E-002,3.1746031746031700E-003,887;
1.0200000000000000E+002,6.6079295154185000E-001,1.6666666666666700E-002,1.6666666666666700E-003,887;
1.0200000000000000E+002,0.0000000000000000E+000,1.6666666666666700E-002,1.6666666666666700E-003,305;
  ];
  t=[
6.5600000000000000E-005;
7.6800000000000000E-005;
5.6320000000000000E-004;
9.2000000000000000E-005;
1.0280000000000000E-004;
1.7920000000000000E-004;
1.0358750000000000E-004;
1.1630000000000000E-004;
1.2560000000000000E-004;
2.0160000000000000E-004;
3.9360000000000000E-004;
1.3610000000000000E-004;
1.6560000000000000E-004;
2.3040000000000000E-004;
1.4953750000000000E-004;
1.7362500000000000E-004;
2.3640000000000000E-004;
1.3353125000000000E-004;
1.2979375000000000E-004;
1.1783750000000000E-004;
1.3217142857142900E-004;
1.3019470899470900E-004;
1.0976798825257000E-004;
1.1202349486049900E-004;
1.5238563598176300E-004;
2.0939463074663300E-004;
1.5199433659048300E-004;
1.3997385699311800E-004;
1.1041241871197800E-004;
1.1642173274596200E-004;
8.6974358974359000E-005;
  ];
end