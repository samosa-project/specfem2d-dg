% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   f = plot_total_energy(OFDIRs, swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP, logscale, titlefig, outputfigpath)
% with:
%   TODO.
% yields:
%   TODO.

function [figureHandle, ke, pe, te] = plot_total_energy(OFDIRs, swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP, logscale, titlefig, outputfigpath)
  addpath('/home/l.martire/Documents/work/mars/mars_is'); % prettyAxes
  if(not(exist('OFDIRs')))
    OFDIRsProvided=0;
  else
    OFDIRsProvided=1;
  end
  if(not(exist('swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP')))
    swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP = 0;
  else
    if(not(ismember(swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP,[0,1,2,4])))
      disp(['[',mfilename,', INFO] swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP not in [0,1,2,4], setting it to 0 (kinetic energy).']);
      swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP = 0;
    end
  end
  if(not(exist('logscale')))
    logscale=1;
  end
  if(not(exist('titlefig')) | strcmp(titlefig,'auto'))
    titlefig_provided = 0;
  else
    titlefig_provided = 1;
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

  
  % prepare labels
  switch(swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP)
    case 0
      titlefig_local='Simulation Kinetic Energy';
      ylab = 'kinetic energy [J]';
    case 1
      titlefig_local='Simulation Potential Energy';
      ylab = 'potential energy [J]';
    case 2
      titlefig_local='Simulation Total (K+P) Energy';
      ylab = 'total energy [J]';
    case {3,4}
      titlefig_local='Simulation Energies';
      ylabL_base = ['potential energy [J]'];
      ylabR_base = ['kinetic energy [J]'];
      if(swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP==3)
        LS_left='-';
        LS_right='--';
        prefixL = LS_left;
        prefixR = LS_right;
        ylabL_base = [ylabL_base, ' (\texttt{',prefixL,'})'];
        ylabR_base = [ylabR_base, ' (\texttt{',prefixR,'})'];
      else
        LS_left='-';
        LS_right='-';
      end
      ylabL=ylabL_base;
      ylabR=ylabR_base;
    otherwise
      error('Set swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP only to either 0, 1, 2, or 3.');
  end

  if(not(titlefig_provided))
    titlefig = titlefig_local;
  end

  figureHandle = figure('units','normalized','outerposition',[0 0 0.5 1]);

  for i=1:numel(OFDIRs)
    % set loading directory, and check it
    OFDIR = OFDIRs{i};
    if(not(strcmp(OFDIR(end),'/'))); OFDIR=[OFDIR,'/']; end;
    if(not(exist(OFDIR)==7))
      disp(['[',mfilename,', INFO] ''',OFDIR,''' is not a directory. Skipping.']);
      continue;
    end
    energypath = [OFDIR,energyfilename];
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
    
    % If LNS, update values given the file format (see energy file, iterate_time.f90, and compute_energy.f90).
    if(readExampleFiles_extractParam([OFDIR,'input_parfile'],'USE_LNS','bool'))
%       ke(2:end) = ke(1)+ke(2:end);
%       pe(2:end) = pe(1)+pe(2:end);
%       te(2:end) = te(1)+te(2:end);
      ke(1)=0;
      pe(1)=0;
      te(1)=0;
    end
    
    % define what goes where
    switch(swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP)
      case 0
        y=ke;
      case 1
        y=pe;
      case 2
        y=te;
      case {3,4}
        yL = pe;
        yR = ke;
      otherwise
        error('ddddd');
    end
    
    % do the plot
    switch(swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP)
      case {3,4}
        % classical but two yaxis or two subplots
        
        if(swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP==4)
          % prepare axes for subplots
          if(i==1)
            topscale = 1;
            if(isa(titlefig,'cell'))
              topscale = numel(titlefig);
            end
            axxx = tight_subplot(2, 1, 0.03, [0.07,topscale*0.03], [0.15, 0.011]); % not mandatory, but prettier
%           else
%             f=gcf; axxx = f.Children;
          end
        end
        
        % switch to good axis
        if(swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP==3)
          yyaxis left;
        else
          axes(axxx(1));
        end
        % plot
        if(logscale && peak2peak(yL)>0)
          if(sign(min(yL))==sign(max(yL)))
            % if no change of sign, plot normal log
            semilogy(t, yL, 'displayname',simulationname, 'linestyle', LS_left); hold on;
            y=yL; sel=(y>0 & (not(max(t)>0)|(t>=0))); cMinMax=[min(yL(sel)),max(yL(sel))]; ylim(niceLogYLim(cMinMax)); sprintf('%.1e ',cMinMax);
          else
            % if change of sign, symlog (https://fr.mathworks.com/matlabcentral/fileexchange/57902-symlog)
            plot(t, yL, 'displayname',simulationname, 'linestyle', LS_left); hold on;
%             autoSymLog(axxx(1), yL); ylabL={ylabL_base,'(symlog scale)'};
          end
        else
          plot(t, yL, 'displayname',simulationname, 'linestyle', LS_left); hold on;
        end
        ylabel(ylabL); set(gca, 'ycolor', 'k');
        % switch to good axis
        if(swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP==3)
          yyaxis right;
        else
          set(axxx(1),'xticklabels',[]);
          axes(axxx(2));
        end
        % plot
        if(logscale && peak2peak(yR)>0)
          if(sign(min(yL))==sign(max(yL)))
            semilogy(t, yR, 'displayname',simulationname, 'linestyle', LS_right); hold on;
            y=yR; sel=(y>0 & (not(max(t)>0)|(t>=0))); cMinMax=[min(yR(sel)),max(yR(sel))]; ylim(niceLogYLim(cMinMax)); sprintf('%.1e ',cMinMax);
          else
            % if change of sign, symlog (https://fr.mathworks.com/matlabcentral/fileexchange/57902-symlog)
            plot(t, yR, 'displayname',simulationname, 'linestyle', LS_left); hold on;
%             autoSymLog(axxx(2), yR); ylabR={ylabR_base,'(symlog scale)'};
          end
        else
          plot(t, yR, 'displayname',simulationname, 'linestyle', LS_right); hold on;
        end
        ylabel(ylabR); set(gca,'ycolor','k');
        % link axes t
        linkaxes(axxx, 'x');
        
      case {0,1,2}
        % if classical plot
        if(logscale && peak2peak(y)>0)
          semilogy(t, y, 'displayname',simulationname); hold on;
        else
          plot(t, y, 'displayname',simulationname); hold on;
        end
      otherwise
        error('kuk');
    end
  end
  
  xlabel('time $t$ [s]');
  switch(swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP)
    case {0,1,2}
      ylabel(ylab);
  end
  
  % update xlim (and eventually ylim)
%   sel = (y>0 & (not(max(t)>0)|(t>=0)));
  if(max(t)>=0)
    sel = (t>=0);
  else
    sel = 1:numel(t);
  end
  kek = find(sel==1);
  if(kek(1)>1)
    sel(kek(1)-1) = 1; % for safety, add the previous time step to include 0
  end
  xlim([min(t(sel)),max(t(sel))]);
%   if(logscale)
%     % create nice logarithmic ylims
%     if(swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP~=3)
%       ylim(niceLogYLim([min(y(sel)),max(y(sel))]));
%     end
%   end

  if(swtchE_0K_1P_2T_3yyaxisKP_4sbpltKP==4)
    axes(axxx(1));
%     niceFormatForYTickLabels(gcf);
  end
%   title({['Total Simulation Energy'],['(',titlefig,')']});
  title(titlefig);
  legend('location', 'best');
  prettyAxes(figureHandle);
  if(outputfigpath_provided)
    customSaveFig(outputfigpath);
  else
    % customSaveFig([OFDIR,'total_energy']);
    customSaveFig(['~/TMP/total_energy']);
  end
end

function newYlim = niceLogYLim(curMinMax)
  if(curMinMax(1)<=0)
    suby = sort(y(sel));
    curMinMax(1) = suby(find(suby>0,1,'first'));
  end
  powTen = 10.^(floor(log10(curMinMax)));
  mM_multiplier = curMinMax./powTen;
  mM_multiplier = [floor(mM_multiplier(1)), ceil(mM_multiplier(2))];
  newYlim = mM_multiplier.* powTen;
end