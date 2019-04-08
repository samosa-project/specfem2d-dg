% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
% clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 6); % Default at 6.
set(0, 'defaultTextFontSize', 18);
set(0, 'defaultAxesFontSize', 16); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

energyfile='energy.dat';

OFDIR = input('  Folder containing the energy.dat file? > ', 's');
if(not(strcmp(OFDIR(end),'/'))); OFDIR=[OFDIR,'/']; end;
if(not(exist(OFDIR)==7))
  error(['not a directory']);
end
energypath = [OFDIR,energyfile];
if(not(exist(energypath)==2))
  error(['energy file ''',energypath,'''does not exist']);
end

% logscale=-1;
% while(not(logscale==0 || logscale==1))
%   logscale=input('  Normal y-scale (0) or log y-scale (1)? > ');
% end
logscale=1;

spl=split(OFDIR,filesep);
simulationname=regexprep(spl{end-2},'_','\\_');

% Load the data.
EF = importdata(energypath);
t  = EF.data(:, 1);
ke = EF.data(:, 2);
pe = EF.data(:, 3);
te = EF.data(:, 4);

f = figure('units','normalized','outerposition',[0 0 0.5 1]);
if(logscale)
  semilogy(t, te, 'displayname',simulationname);
else
  plot(t, te, 'displayname',simulationname);
end
xlabel('t (s)');
ylabel('total energy (J)');
title({['Total Simulation Energy'],['(',simulationname,')']});
xlim([min(t),max(t)]);
prettyAxes(f);
customSaveFig([OFDIR,'total_energy']);