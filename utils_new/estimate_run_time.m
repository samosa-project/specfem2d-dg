% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Estimate new run time based on previous runs' saved data.
% Last modified: See file metadata.
% Usage:         Fill section "Parameters" according to the calculation to
%                be performed.
%                Use the scripts in './utils_new/extract_run_information'
%                to fill data used for estimation (function 'load' below).
% Notes:         The method (N-D nearest point search) is approximate and
%                does not take into account the fact that some parameters
%                have more impact than others. All in all, it is a very
%                rough and very approximate method.

% clear all
% close all
clc
format longG;

[data, t, inf]=load(); % Load data (see function below).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nstations   = 4*90;
nstepseismo = 25;

neltot      = 1794000;
neldg       = 1794000;

nstepsnap   = 500;
nsteptot    = 120000;

nproc       = 224;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of estimate.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weigth=[1,100,5,1,50]; % Importance of each parameter.
meandata=mean(data,1);

point=[nstations,neldg/neltot,nstepsnap/nsteptot,nstepseismo/nsteptot,neltot/nproc]; % Format point as data format. Currently [% elements as DG, % timesteps as snapshots, elements per proc].

disp(['[',mfilename,'] ',num2str(neltot),' elements including ',num2str(neldg),' DG elements. ',num2str(nsteptot),' time steps. ',num2str(nstations),' stations sampling every ',num2str(nstepseismo),' iterations. Snapshots taken every ',num2str(nstepsnap),' iterations. ',num2str(nproc),' CPUs.']);
disp(" ");
disp(['[',mfilename,']                     [        n_stations       percent_DG     percent_snap    percent_synth n_elems_per_proc]']);
disp(strcat("[",mfilename,"] Current point:      [ ", sprintf("%17.3e", point), "]."));
% idp=dsearchn(data,point);
idp=dsearchn((data-meandata)./weigth,(point-meandata)./weigth); % Optimised search.
disp(strcat("[",mfilename,"] Closest data point: [ ", sprintf("%17.3e", data(idp,:)), "] (",inf{idp}{2},")."));
disp(strcat("[",mfilename,"] Weights:            [ ", sprintf("%17.f", weigth), "]."));
time=t(idp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(" ");

hms=fix(mod(time*neltot*nsteptot,[0,3600,60])./[3600,60,1]); cputimestr=strcat(sprintf('%5.f',hms(1)),"h ",sprintf('%2.f',hms(2)),"m ",sprintf('%2.f',hms(3)),'s');
hms=fix(mod(time*neltot*nsteptot/nproc,[0,3600,60])./[3600,60,1]); realtimestr=strcat(sprintf('%5.f',hms(1)),"h ",sprintf('%2.f',hms(2)),"m ",sprintf('%2.f',hms(3)),'s');

disp(strcat("[",mfilename,"] Expected time per element, per iteration:    ",sprintf("%.3e", time), " s."));
disp(strcat("[",mfilename,"] Expected run time:                        ",cputimestr, " (CPU)."));
disp(strcat("[",mfilename,"]                                           ",realtimestr, " (real)."));

disp(" ");

disp(['[',mfilename,', WARNING] Recall the method used for estimation is very rough and approximate. Do not take the estimation for granted.']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function containing data.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,t,RUNINFO]=load()
  col_nbelts=1;
  col_nbeltsdg=2;
  col_cfl=3;
  col_nsteps=4;
  col_nprocs=5;
  col_snapfreq=6;
  col_nstats=7;
  col_synthfreq=8;
  col_realtime=9;
  i=1;
  RUN_RAWDATA(i,:)=[   2500    2500 0.404   2000    4   50   0  0     82]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[   2500    2500 0.404   2000   16   50   0  0     24]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[   2500    2500 0.404   2000  256   50   0  0     11]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  10000   10000 0.404   4000    4   50   0  0    920]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  10000   10000 0.404   4000   16   50   0  0    257]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  10000   10000 0.404   4000  256   50   0  0     28]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  40000   40000 0.404   8000    4   50   0  0   8287]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  40000   40000 0.404   8000   16   50   0  0   2326]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  40000   40000 0.404   8000  256   50   0  0    157]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  40000   40000 0.404   8000  512   50   0  0    126]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[  40000   40000 0.404   8000 1024   50   0  0    123]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 160000  160000 0.404  16000  256   50   0  0   1361]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 160000  160000 0.404  16000  512   50   0  0    828]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 160000  160000 0.404  16000 1024   50   0  0    576]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 640000  640000 0.404  32000  256   50   0  0  11963]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 640000  640000 0.404  32000  512   50   0  0   6945]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 640000  640000 0.404  32000 1024   50   0  0   4728]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[2560000 2560000 0.404  64000  256    0   0  0  85460]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[2560000 2560000 0.404  64000  512    0   0  0  41534]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[2560000 2560000 0.404  64000 1024    0   0  0  18854]; RUNINFO{i}={-1,'test on Curie'}; i=i+1;
  RUN_RAWDATA(i,:)=[2016000 1830000 0.471 120000  256  200 116 50 124902]; RUNINFO{i}={583041,'OKQ45'}; i=i+1;
  RUN_RAWDATA(i,:)=[2016000 1830000 0.471 120000  256  200 116 50 123034]; RUNINFO{i}={586984,'OKQ0'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 227000  150000 0.443  30000  256  500   0  0   2920]; RUNINFO{i}={593959,'SH soft final'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 227000  150000 0.443  30000  256  500   0  0   2980]; RUNINFO{i}={593960,'SH hard final'}; i=i+1;
  RUN_RAWDATA(i,:)=[2892000 2832000 0.376  27000  512  200   0  0  23240]; RUNINFO{i}={594536,'StratoExplo66June1200 - not finished'}; i=i+1;
  RUN_RAWDATA(i,:)=[  13060   13060 0.156  29400   16  200   0  0   5025]; RUNINFO{i}={595104,'StratoExplo66June1200 - test'}; i=i+1;
  RUN_RAWDATA(i,:)=[2892000 2832000 0.376  40150  512  300   0  0  34470]; RUNINFO{i}={595500,'StratoExplo66June1200 - not finished'}; i=i+1;
  RUN_RAWDATA(i,:)=[3142000 3082000 0.376 140100  560  300   0  0 110028]; RUNINFO{i}={597316,'StratoExplo66June1200 - til n=140100'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 227000  150000 0.443  15750  256  500   2 50   1542]; RUNINFO{i}={610736,'SH soft final redone for data comparison'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 227000  150000 0.443  30000  256  500 102 50   3097]; RUNINFO{i}={616368,'SH soft final redone for data comparison'}; i=i+1;
  RUN_RAWDATA(i,:)=[  78000       0 0.443  30000  256  500 102 50    795]; RUNINFO{i}={623195,'SH soft final redone w/o fluid'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 748360  748360 0.345   8000   32  250 114 50  27994]; RUNINFO{i}={624650,'StratoBaro full with IBF'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 748360  748360 0.345   5750   64  250 114 50   9225]; RUNINFO{i}={637450,'StratoBaro re-run with EBF'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 748360  748360 0.345  16000   96  250 114 50  17354]; RUNINFO{i}={641616,'StratoBaro full with EBF'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 266162  266162 0.345  16000  128  250 114 50   5904]; RUNINFO{i}={656505,'microbaroms periodic with EBF'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 266162  266162 0.345  16000  256  250 114 50   4076]; RUNINFO{i}={655513,'microbaroms periodic with EBF'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 102150       0 0.443  30000   32  500  80 50    757]; RUNINFO{i}={660223,'SH soft axisym'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 102150       0 0.443  30000   32  500  80 50    744]; RUNINFO{i}={661601,'SH hard axisym'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 225000  150000 0.443  12250   32  500  82 50  10802]; RUNINFO{i}={668888,'SH soft first layer tweaked stopped 12kit'}; i=i+1;
  RUN_RAWDATA(i,:)=[ 224000  150000 0.443  20000   64  500  82 50   7802]; RUNINFO{i}={668888,'SH soft first layer retweaked'}; i=i+1;
  RUN_RAWDATA(i,:)=[1866000 1680000 0.471 110000  256 5000 116 50 102012]; RUNINFO{i}={668833,'OKQ0 redone'}; i=i+1;
  RUN_RAWDATA(i,:)=[1866000 1680000 0.471 110000  256 5000 116 50 101666]; RUNINFO{i}={668844,'OKQ45 redone'}; i=i+1;
  col_dgpercent=1;
  col_snappercent=2;
  col_synthpercent=3;
  col_nbeltspproc=4;
  col_hcpu=5;
  col_cputimepelpit=6;
  RUN_LV1DATA(:,col_dgpercent)=RUN_RAWDATA(:,col_nbeltsdg)./RUN_RAWDATA(:,col_nbelts);
  RUN_LV1DATA(:,col_snappercent)=RUN_RAWDATA(:,col_snapfreq)./RUN_RAWDATA(:,col_nsteps);
  RUN_LV1DATA(:,col_synthpercent)=RUN_RAWDATA(:,col_synthfreq)./RUN_RAWDATA(:,col_nsteps);
  RUN_LV1DATA(:,col_nbeltspproc)=RUN_RAWDATA(:,col_nbelts)./RUN_RAWDATA(:,col_nprocs);
  RUN_LV1DATA(:,col_hcpu)=RUN_RAWDATA(:,col_realtime).*RUN_RAWDATA(:,col_nprocs)/3600;
  RUN_LV1DATA(:,col_cputimepelpit)=(RUN_LV1DATA(:,col_hcpu)*3600./RUN_RAWDATA(:,col_nbelts))./RUN_RAWDATA(:,col_nsteps);

  x=[RUN_RAWDATA(:,col_nstats),RUN_LV1DATA(:,col_dgpercent),RUN_LV1DATA(:,col_snappercent),RUN_LV1DATA(:,col_synthpercent),RUN_LV1DATA(:,col_nbeltspproc)];
  t=RUN_LV1DATA(:,col_cputimepelpit);
end