% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

% clear all;
% close all;
clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 18); set(0, 'defaultAxesFontSize', 18);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

f0=100;
% wind=0;
wind=102;
% wind=204;
c=340;

M=wind/c;

% According to Ostashev D:
% dx=dy=dr<lambda*(1-M)/Npts=c*(1-M)/(Npts*f0)
% dt<(1-M)/((1+M)*Npts*f0)
Npts=2;
dxmax=c*(1-M)/(Npts*f0);
dtmax=(1-M)/((1+M)*Npts*f0);

NSTATIONS=15;

alpha=linspace(0,pi,NSTATIONS);

k=2*pi*f0/c;
r=20/k;

% A=0.98;
A=1;

p_th= A*(sqrt(1-M^2*sin(alpha).^2) - M*cos(alpha)) ./ (sqrt(2*pi*k*r)*(1-M^2)*(1-M^2*sin(alpha).^2).^0.75);
p0_th= A ./ (sqrt(2*pi*k*r));

figure();
plot(alpha*180/pi,p_th/p0_th,'k','displayname','Ostashev''s (72)');
set(gca,'ticklabelinterpreter','latex');
set(gca,'xtick',180*[0,0.25,0.5,0.75,1]);
xlim([0,180]);
hold on;
title(['$f_0=',num2str(f0),'$ Hz, $c=',num2str(c),'$, $w=',num2str(wind),'$, $M=',num2str(M),'$, $k=',num2str(k),'$, $r=',num2str(r),'$, $kr=',num2str(k*r),'$']);
ylabel(['normalised sound pressure amplitude']);
xlabel(['azimuth $\alpha$ [deg]']);
legend('location','northwest');
grid on;

if(0)
  p0=Zamp; % Use this after having ran synth_load on Mach=0 simulation and selected the right station.
  p0_ampl=max(p0)-min(p0)
  
  p=Zamp; % Use this after having ran synth_load on relevant simulation.
  p_ampl=(max(p')-min(p'))
  angglle=atan(ystattab./xstattab)*180/pi;
  angglle(angglle<0)=angglle(angglle<0)+180;
  angglle(end)=180;
  
  plot(angglle,p_ampl/p0_ampl,'displayname','LNS');
  plot(angglle,p_ampl/p0_ampl,'displayname','FNS');
end


clc;
disp(['nreceiversets = ',num2str(numel(alpha))]);
disp(['# Orientation']);
disp(['anglerec                        = 0.0d0          # Angle to rotate components at receivers.']);
disp(['rec_normal_to_surface           = .false.        # Base anglerec normal to surface (external mesh and curve file needed).']);
for ang=alpha
  disp(['# az ', num2str(ang*180/pi),'°']);
  disp(['nrec = 1']);
  disp(['xdeb = ',char(sprintf('%.16g',r*cos(ang)))]);
  disp(['zdeb = ',char(sprintf('%.16g',r*sin(ang)))]);
  disp(['xfin = 0.']);
  disp(['zfin = 0.']);
  disp(['record_at_surface_same_vertical = .false.']);
end

disp(' ');

disp(['[',mfilename,'] Check if dx<',num2str(dxmax),' and dt<',num2str(dtmax),'']);

% GMSH point definition
if(0)
  for angi=1:length(alpha)
    disp(['Point(',num2str(4+angi),') = {',char(sprintf('%.18g',r*cos(alpha(angi)))),',',char(sprintf('%.18g',r*sin(alpha(angi)))),',0};']);
  end
  for angi=1:length(alpha)
    disp(['Point{',num2str(4+angi),'} In Surface{1};']);
  end
end
