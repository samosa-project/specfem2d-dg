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

clear all;
close all;
clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 10); % Default at 6.
set(0, 'defaultTextFontSize', 16);
set(0, 'defaultAxesFontSize', 16); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run-specific.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xminmax, zminmax, interface_z, Xsource, debfin, d, name]=mars_insight();
[xminmax, zminmax, interface, Xsource, debfin, d, name]=tir_de_mine();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automatic.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=get_domain_and_stations(xminmax, zminmax, interface, Xsource, debfin, d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ask if ok.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ok=-1;
while(not(numel(ok)==1 & ismember(ok,[0,1])))
  ok=input(['[', mfilename, '] Are those ',num2str(sum(N)),' stations ok (0 for no, 1 for yes)? > ']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display copy-pastable.      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ok)
  clc;
  disp(['nreceiversets = ',num2str(numel(d))]);
  disp(['# Orientation']);
  disp(['anglerec              = 0.0d0   # Angle to rotate components at receivers.']);
  disp(['rec_normal_to_surface = .false. # Base anglerec normal to surface (external mesh and curve file needed).']);
  for i=1:numel(d)
    disp(['# ', name{i}]);
    disp(['nrec = ',num2str(N(i))]);
    disp(['xdeb = ',char(sprintf('%17.16g',debfin(i,1,1)))]);
    disp(['zdeb = ',char(sprintf('%17.16g',debfin(i,2,1)))]);
    disp(['xfin = ',char(sprintf('%17.16g',debfin(i,1,2)))]);
    disp(['zfin = ',char(sprintf('%17.16g',debfin(i,2,2)))]);
    disp(['record_at_surface_same_vertical = .false.']);
  end

  disp(' ');
else
  disp(['[', mfilename, '] Stations not ok, stopping script.']);
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N]=get_domain_and_stations(xminmax, zminmax, interface, Xsource, debfin, delta)
  if(not( numel(xminmax)==numel(zminmax) & numel(zminmax)==2))
    error('kek');
  end
  if(not(size(debfin,1)==numel(delta)))
    error('kek');
  end
  lid=numel(delta);
  figure('units','normalized','outerposition',[0 0 1 1]);
  set(gca,'Color','k');
  set(gca,'GridColor','white');
  rectangle('Position',[xminmax(1), zminmax(1), diff(xminmax), diff(zminmax)],'edgecolor','white','linewidth',1); hold on;
  line(interface(1,:),interface(2,:),'color','white','linewidth',1);
  pct=0.05; margin=min(pct*diff(xminmax),pct*diff(zminmax));
  axis([xminmax(1)-margin, xminmax(2)+margin, zminmax(1)-margin, zminmax(2)+margin]);
  grid on; box on;
  plot(Xsource(1),Xsource(2),'r+');
  for i=1:lid
    Lline=norm(diff(debfin(i,:,:),1,3));
    if(Lline==0 || delta(i)==0)
      N(i)=1;
    else
      N(i)=Lline/delta(i)+1;
    end
    xz=[linspace(debfin(i,1,1),debfin(i,1,2),N(i)); ...
        linspace(debfin(i,2,1),debfin(i,2,2),N(i))]
    plot(xz(1,:),xz(2,:),'g.');
  end
  daspect([1 1 1]);
  xlabel(['$x$ [m]']); ylabel(['$z$ [m]']);
  set(gca,'ticklabelinterpreter','latex');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run configurations.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xminmax, zminmax, interface, Xsource, debfin, d, name]=mars_insight()
  Xsource=[0,40e3];
  xminmax=[-100,100]*1e3;
  zminmax=[-5,60]*1e3;
  spacingstations=4e3;
  verticalaway_x=30e3;
  interface=[-1e9,0;1e9,0];

  lid=0;
  lid=lid+1;
  % d(lid)=2.5e3;
  d(lid)=spacingstations;
  debfin(lid,1,:)=[Xsource(1),Xsource(1)]; % xdeb xfin
  debfin(lid,2,:)=[0,Xsource(2)]; % zdeb zfin
  name{lid} = ['vertical under source'];

  lid=lid+1;
  % d(lid)=2.5e3;
  d(lid)=spacingstations;
  debfin(lid,1,:)=[Xsource(1)-verticalaway_x,Xsource(1)-verticalaway_x]; % xdeb xfin
  debfin(lid,2,:)=[0,Xsource(2)]; % zdeb zfin
  name{lid} = ['vertical ',num2str(-verticalaway_x),' m away from source'];

  lid=lid+1;
  % d(lid)=2.5e3;
  d(lid)=spacingstations;
  debfin(lid,1,:)=[Xsource(1)+verticalaway_x,Xsource(1)+verticalaway_x]; % xdeb xfin
  debfin(lid,2,:)=[0,Xsource(2)]; % zdeb zfin
  name{lid} = ['vertical ',num2str(verticalaway_x),' m away from source'];

  lid=lid+1;
  % d(lid)=2.5e3;
  d(lid)=spacingstations;
  debfin(lid,1,:)=[xminmax(1)+d(lid),xminmax(2)-d(lid)]; % xdeb xfin
  debfin(lid,2,:)=[1,1]; % zdeb zfin
  name{lid} = ['horizontal over ground'];
  lid=lid+1;
  d(lid)=d(lid-1);
  debfin(lid,1,:)=debfin(lid-1,1,:);
  debfin(lid,2,:)=-1*[1,1]; % zdeb zfin
  name{lid} = ['horizontal under ground'];
end

function [xminmax, zminmax, interface, Xsource, debfin, d, name]=tir_de_mine()
  Xsource=[0,-25];
  xminmax=[-440,240];
  zminmax=[-190,150];
  ISAEIS2_z=65;
  ISAEIS3_z=104;
  to_gla08_z=(195.-167.);
  to_gla08_r=187.;
  to_gla08_x=(to_gla08_r^2.-to_gla08_z^2.)^0.5;
  xr0 = -to_gla08_x/to_gla08_r; zr0 = -to_gla08_z/to_gla08_r; zx0 = -to_gla08_z/to_gla08_x;
  interface=[-400,10;400*zx0,-10*zx0];

  lid=0;
  lid=lid+1; d(lid)=0; r=0;   name{lid} = ['above source (@r=0)',num2str(r)]; debfin(lid,1,:)=[1,1]*r*xr0; debfin(lid,2,:)=[1,1]*r*zr0;
  lid=lid+1; d(lid)=0; r=100; name{lid} = ['GLA04 @r=',num2str(r)]; debfin(lid,1,:)=[1,1]*r*xr0; debfin(lid,2,:)=[1,1]*r*zr0;
  lid=lid+1; d(lid)=0; r=193; name{lid} = ['GLA05 @r=',num2str(r)]; debfin(lid,1,:)=[1,1]*r*xr0; debfin(lid,2,:)=[1,1]*r*zr0;
  lid=lid+1; d(lid)=0; r=150; name{lid} = ['GLA06 @r=',num2str(r)]; debfin(lid,1,:)=[1,1]*r*xr0; debfin(lid,2,:)=[1,1]*r*zr0;
  lid=lid+1; d(lid)=0; r=5;   name{lid} = ['GLA07 @r=',num2str(r)]; debfin(lid,1,:)=[1,1]*r*xr0; debfin(lid,2,:)=[1,1]*r*zr0;
  lid=lid+1; d(lid)=0; r=187; name{lid} = ['GLA08 @r=',num2str(r)]; debfin(lid,1,:)=[1,1]*r*xr0; debfin(lid,2,:)=[1,1]*r*zr0; idGLA08=lid;
  lid=lid+1; d(lid)=0; r=185; name{lid} = ['GLA09 @r=',num2str(r)]; debfin(lid,1,:)=[1,1]*r*xr0; debfin(lid,2,:)=[1,1]*r*zr0;
  lid=lid+1; d(lid)=0; r=50;  name{lid} = ['GLA11 @r=',num2str(r)]; debfin(lid,1,:)=[1,1]*r*xr0; debfin(lid,2,:)=[1,1]*r*zr0;
  debfin(:,2,1)=debfin(:,2,1)+1; % beginning above ground
  debfin(:,2,2)=debfin(:,2,2)-1; % end under ground
  d=d+2; % separate correctly
  lid=lid+1; d(lid)=0; r=187; name{lid} = ['ISAE IS Sensor 2, above GLA08 (@r=',num2str(r),'), z=+',num2str(ISAEIS2_z),' (ISAEIS ground is GLA08 z=+1m)']; debfin(lid,:,:)=debfin(idGLA08,:,:); debfin(lid,2,:)=debfin(lid,2,:)+ISAEIS2_z;
  lid=lid+1; d(lid)=0; r=187; name{lid} = ['ISAE IS Sensor 3, above GLA08 (@r=',num2str(r),'), z=+',num2str(ISAEIS3_z),' (ISAEIS ground is GLA08 z=+1m)']; debfin(lid,:,:)=debfin(idGLA08,:,:); debfin(lid,2,:)=debfin(lid,2,:)+ISAEIS3_z;
end