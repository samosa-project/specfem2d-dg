% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         This script can be used to provide illustration of the various bottom forcings implemented in the DG extension.

clear all;
% close all
clc;
format compact;
set(0, 'DefaultLineLineWidth', 3); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 8); % Default at 6.
set(0, 'defaultTextFontSize', 22);
set(0, 'defaultAxesFontSize', 22); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

%% Microbarom
perio=7;
lambdo=200;
MICROBAROM_MAXTIME=8*perio;
MICROBAROM_RANGE=15e3;
dt=1.5d-2;
% t=max(0,MICROBAROM_MAXTIME-4*perio):dt:MICROBAROM_MAXTIME+perio/4;
t=0:dt:3*perio;
x=linspace(MICROBAROM_RANGE-15*lambdo,MICROBAROM_RANGE+lambdo/4,1000);
% x=linspace(-MICROBAROM_RANGE-lambo/4,MICROBAROM_RANGE+lambo/4,1000);

% Build spatial and temporal phases based on normal random walks.
Phix=zeros(length(t),1); Phix(1)=0;
Phit=zeros(length(t),1); Phit(1)=0;
for i=2:length(t)
  % Phase must at most shift pi/2 over one period (otherwise it would be
  % too chaotic). Thus, if we denote the step s, and if the random walk is
  % the most unlucky possible (is all steps done in only one direction),
  % (period/dt)*s must be < pi/2. In other words, d must be <
  % pi*dt/(2*period). If the walk is normal, we thus choose n*sigma = 
  % pi*dt/(2*period), that is sigma = pi*dt/(n*2*period). n=1 is generally
  % a good choice.
  n=0.2;
  sigma=pi*dt/(n*2*perio);
  Phix(i)=Phix(i-1)+normrnd(0,sigma);
  Phit(i)=Phit(i-1)+normrnd(0,sigma);
end

% Apodisation
% n=10; apox=0.25.*(1.-erf((x-MICROBAROM_RANGE+n*lambdo)/(0.5*n*lambdo))).*(1+erf((x+MICROBAROM_RANGE-n*lambdo)/(0.5*n*lambdo)));
n=10; apox=(1-((-abs(x)+MICROBAROM_RANGE)/(n*lambdo)-1).^2 .*(abs(x)>MICROBAROM_RANGE-n*lambdo)).*(abs(x)<MICROBAROM_RANGE);
% n=1; apot=0.5.*(1.-erf((t-MICROBAROM_MAXTIME+n*perio)/(0.5*n*perio)));
n=2; apot=(t<MICROBAROM_MAXTIME-n*perio)+ (1-((t-MICROBAROM_MAXTIME)/(n*perio)+1).^2) .*(t>MICROBAROM_MAXTIME-n*perio).*(t<MICROBAROM_MAXTIME);
% n=1; apot0=0.5*(t/(n/2*perio)).^2 .* (t<n/2*perio) + (0.5+0.5*(1-(2-t/(n/2*perio)).^2)) .*(t>n/2*perio).*(t<n*perio);
n=1; apot0=2*(t/(n*perio)).^2 .* (t<n/2*perio) + ((4*t)/(n*perio)-2*(t/(n^2*perio)).^2-1) .*(t>n/2*perio).*(t<n*perio) + (t>n*perio);
apot=apot.*apot0;

[X,T]=meshgrid(x,t);

Fx=sin(2.*pi*X/lambdo+Phix);
Ft=sin(2.*pi*T/perio+Phit);
Fx=Fx.*(abs(X)<MICROBAROM_RANGE);
Ft=Ft.*(T<MICROBAROM_MAXTIME);
Fx=Fx.*apox;
Ft=Ft.*apot';
F=Fx.*Ft;

figure();
surf(X,T,F,'EdgeColor','none');
view([0,0,1]);
colormap('jet');
xlim([min(x),max(x)]);
ylim([min(t),max(t)]);
xlabel('$x$');
ylabel('$t$');