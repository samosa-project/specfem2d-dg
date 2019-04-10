% Author:        Léo Martire.
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.


function f = plot_total_energy(OFDIRs, switchE_0K_1P_2T, logscale, titlefig, outputfigpath)
  addpath('/home/l.martire/Documents/work/mars/mars_is'); % prettyAxes
  if(not(exist('OFDIRs')))
    OFDIRsProvided=0;
  else
    OFDIRsProvided=1;
  end
  if(not(exist('switchE_0K_1P_2T')))
    switchE_0K_1P_2T = 0;
  else
    if(not(ismember(switchE_0K_1P_2T,[0,1,2])))
      disp('switchE_0K_1P_2T not in [0,1,2], setting it to 0 (kinetic energy).');
      switchE_0K_1P_2T = 0;
    end
  end
  if(not(exist('logscale')))
    logscale=1;
  end
  if(not(exist('titlefig')))
    titlefig_provided=0;
  else
    titlefig_provided=1;
  end
  if(not(exist('outputfigpath')))
    outputfigpath_provided=0;
  else
    outputfigpath_provided=1;
  end
  
  % clear all;
  % clc;
  format compact;
  set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
  set(0, 'DefaultLineMarkerSize', 6); % Default at 6.
  set(0, 'defaultTextFontSize', 18);
  set(0, 'defaultAxesFontSize', 16); % Default at 10.
  set(0, 'DefaultTextInterpreter', 'latex');
  set(0, 'DefaultLegendInterpreter', 'latex');

  energyfilename='energy.dat';
  
  if(OFDIRsProvided)
    % nothing
  else
    % get OFDIRs
    many0orone1=-1;
    while(not(numel(many0orone1) & ismember(many0orone1,[0,1])))
      many0orone1 = input(['[',mfilename,'] Plot many (0) or one (1) energy curves? > ']);
    end

    OFDIRs={};
    if(many0orone1)
      in = input('  Folder containing the energy.dat file? > ', 's');
      if(not(exist(in)==7))
        error(['this is not a directory, try again']);
      end
      OFDIRs{1} = in;
    else
      in='-1';
      i=1;
      stop=0;
      while(not(stop))
        in = input(['[',mfilename,'] Next folder containing the energy.dat file (or 0 to stop)? > '], 's');
        if(numel(in)==1)
          if(str2num(in)==0)
            stop=1;
            continue;
          end
        end
        OFDIRs{i}=in;
        i=i+1;
      end
      OFDIRs
    end
  end

  f = figure('units','normalized','outerposition',[0 0 0.5 1]);
  for i=1:numel(OFDIRs)
    OFDIR = OFDIRs{i};
    if(not(strcmp(OFDIR(end),'/'))); OFDIR=[OFDIR,'/']; end;
    if(not(exist(OFDIR)==7))
      disp(['[',mfilename,', INFO] ''',OFDIR,''' is not a directory. Skipping.']);
      continue;
    end
    energypath = [OFDIR,energyfilename]
    if(not(exist(energypath)==2))
      disp(['[',mfilename,', INFO] Energy file ''',energypath,'''does not exist. Skipping.']);
      continue;
    end

    % Prepare a displayname.
    spl=split(OFDIR,filesep);
    simulationname = [spl{end-2},'__',spl{end-1}];
    simulationname = regexprep(simulationname,'_','\\_');

    % Load the data.
    EF = importdata(energypath);
    t  = EF.data(:, 1);
    ke = EF.data(:, 2);
    pe = EF.data(:, 3);
    te = EF.data(:, 4);
    
    switch(switchE_0K_1P_2T)
      case 0
        y=ke;
        titlefig_local='Simulation Kinetic Energy';
        ylab = 'kinetic energy (J)';
      case 1
        y=pe;
        titlefig_local='Simulation Potential Energy';
        ylab = 'potential energy (J)';
      case 2
        y=te;
        titlefig_local='Simulation Total (K+P) Energy';
        ylab = 'total energy (J)';
      otherwise
        error('Set switchE_0K_1P_2T only to either 0, 1, or 2');
    end
    if(not(titlefig_provided))
      titlefig=titlefig_local;
    end

    if(logscale)
      semilogy(t, y, 'displayname',simulationname);
    else
      plot(t, y, 'displayname',simulationname);
    end
    hold on;
  end
  xlabel('t (s)');
  ylabel(ylab);
%   title({['Total Simulation Energy'],['(',titlefig,')']});
  title(titlefig);
  
  sel = (t>=0);
  kek=find(sel==1);
  sel(kek(1)-1) = 1; % for safety, add the previous time step to include 0
  if(logscale)
    % create nice logarithmic ylims
    curMinMax = [min(y(sel)),max(y(sel))];
    if(curMinMax(1)<=0)
      suby = sort(y(sel));
      curMinMax(1) = suby(find(suby>0,1,'first'));
    end
    powTen = 10.^(floor(log10(curMinMax)));
    mM_multiplier = curMinMax./powTen;
    mM_multiplier = [floor(mM_multiplier(1)), ceil(mM_multiplier(2))];
    newYlim = mM_multiplier.* powTen;
    ylim(newYlim);
  end
  xlim([min(t(sel)),max(t(sel))]);
  
  legend('location','best');
  prettyAxes(f);
  if(outputfigpath_provided)
    customSaveFig(outputfigpath);
  else
    % customSaveFig([OFDIR,'total_energy']);
    customSaveFig(['~/TMP/total_energy']);
  end
end