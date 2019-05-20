% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   f = plot_total_energy(OFDIRs, swtchE_0K_1P_2T_4sbpltKP, logscale, titlefig, outputfigpath)
% with:
%   TODO.
% yields:
%   TODO.

function [figureHandle, ke, pe, te] = plot_total_energy(OFDIRs, swtchE_0K_1P_2T_4sbpltKP, logscale, titlefig, outputfigpath, linestyles)
  addpath('/home/l.martire/Documents/work/mars/mars_is'); % prettyAxes
  energyfilename='energy.dat';
  if(not(exist('OFDIRs')))
    OFDIRsProvided=0;
  else
    OFDIRsProvided=1;
  end
  if(not(exist('swtchE_0K_1P_2T_4sbpltKP')))
    swtchE_0K_1P_2T_4sbpltKP = 0;
  else
    if(not(ismember(swtchE_0K_1P_2T_4sbpltKP,[0,1,2,4])))
      disp(['[',mfilename,', INFO] swtchE_0K_1P_2T_4sbpltKP not in [0,1,2,4], setting it to 0 (kinetic energy).']);
      swtchE_0K_1P_2T_4sbpltKP = 0;
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
  if(not(exist('outputfigpath')) | strcmp(outputfigpath,'no') | isempty(outputfigpath))
    outputfigpath_provided=0;
  else
    outputfigpath_provided=1;
  end
  if(not(exist('linestyles')))
    linestyles_provided = 0;
  else
    linestyles_provided = 1;
    if(not(iscell(linestyles)))
      error(['linestyles input must be a cell array']);
    end
    if(not(numel(OFDIRs)==numel(linestyles)))
      error(['linestyles input must be the same length as OFDIRs input']);
    end
  end
  
  if(OFDIRsProvided)
    % nothing
  else
    format compact;
    set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
    set(0, 'DefaultLineMarkerSize', 6); % Default at 6.
    set(0, 'defaultTextFontSize', 18);
    set(0, 'defaultAxesFontSize', 16); % Default at 10.
    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
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
      of_id=1;
      stop=0;
      while(not(stop))
        in = input(['[',mfilename,'] Next folder containing the energy.dat file (or 0 to stop)? > '], 's');
        if(numel(in)==1)
          if(str2num(in)==0)
            stop=1;
            continue;
          end
        end
        OFDIRs{of_id}=in;
        of_id=of_id+1;
      end
      OFDIRs
    end
  end
  
  if(not(linestyles_provided))
    for of=1:numel(OFDIRs)
      linestyles{of} = '-';
    end
  end
  
  % prepare labels
  switch(swtchE_0K_1P_2T_4sbpltKP)
    case 0
      titlefig_local='Simulation Kinetic Energy';
      ylab = 'kinetic energy [J]';
    case 1
      titlefig_local='Simulation Potential Energy';
      ylab = 'potential energy [J]';
    case 2
      titlefig_local='Simulation Total (K+P) Energy';
      ylab = 'total energy [J]';
%     case {3,4}
    case 4
      titlefig_local='Simulation Energies';
      ylabL_base = ['potential energy [J]'];
      ylabR_base = ['kinetic energy [J]'];
      if(swtchE_0K_1P_2T_4sbpltKP==3)
        LS_left='-';
        LS_right='--';
        prefixL = LS_left;
        prefixR = LS_right;
        ylabL_base = [ylabL_base, ' (\texttt{',prefixL,'})'];
        ylabR_base = [ylabR_base, ' (\texttt{',prefixR,'})'];
      end
      ylabL=ylabL_base;
      ylabR=ylabR_base;
    otherwise
      error('Set swtchE_0K_1P_2T_4sbpltKP only to either 0, 1, 2, or 4.');
  end

  if(not(titlefig_provided))
    titlefig = titlefig_local;
  end

  figureHandle = figure('units','normalized','outerposition',[0 0 0.5 1]);

  for i=1:numel(OFDIRs)
    % set loading directory, and check it
    OFDIR = OFDIRs{i};
    localLS = linestyles{i};
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
    switch(swtchE_0K_1P_2T_4sbpltKP)
      case 0
        y=ke;
      case 1
        y=pe;
      case 2
        y=te;
      case 4
        y_top = pe;
        y_bottom = ke;
      otherwise
        error('ddddd');
    end
    
    % do the plot
    switch(swtchE_0K_1P_2T_4sbpltKP)
%       case {3,4}
      case 4
        % classical but two yaxis or two subplots
        
%         if(swtchE_0K_1P_2T_4sbpltKP==4)
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
%         end
        
        % switch to good axis
%         if(swtchE_0K_1P_2T_4sbpltKP==3)
%           yyaxis left;
%         else
        axes(axxx(1));
%         end
        % plot
        if(logscale && peak2peak(y_top)>0)
          if(sign(min(y_top))==sign(max(y_top)))
            % if no change of sign, plot normal log
            semilogy(t, y_top, 'displayname',simulationname, 'linestyle', localLS); hold on;
            y=y_top; sel=(y>0 & (not(max(t)>0)|(t>=0))); cMinMax=[min(y_top(sel)),max(y_top(sel))]; ylim(niceLogYLim(cMinMax)); sprintf('%.1e ',cMinMax);
          else
            % if change of sign, symlog (https://fr.mathworks.com/matlabcentral/fileexchange/57902-symlog)
            plot(t, y_top, 'displayname',simulationname, 'linestyle', localLS); hold on;
%             autoSymLog(axxx(1), yL); ylabL={ylabL_base,'(symlog scale)'};
          end
        else
          plot(t, y_top, 'displayname',simulationname, 'linestyle', localLS); hold on;
        end
        ylabel(ylabL);
%         set(gca, 'ycolor', 'k');
        % switch to good axis
%         if(swtchE_0K_1P_2T_4sbpltKP==3)
%           yyaxis right;
%         else
        set(axxx(1),'xticklabels',[]);
        axes(axxx(2));
%         end
        % plot
        if(logscale && peak2peak(y_bottom)>0)
          if(sign(min(y_top))==sign(max(y_top)))
            semilogy(t, y_bottom, 'displayname',simulationname, 'linestyle', localLS); hold on;
            y=y_bottom; sel=(y>0 & (not(max(t)>0)|(t>=0))); cMinMax=[min(y_bottom(sel)),max(y_bottom(sel))]; ylim(niceLogYLim(cMinMax)); sprintf('%.1e ',cMinMax);
          else
            % if change of sign, symlog (https://fr.mathworks.com/matlabcentral/fileexchange/57902-symlog)
            plot(t, y_bottom, 'displayname',simulationname, 'linestyle', localLS); hold on;
%             autoSymLog(axxx(2), yR); ylabR={ylabR_base,'(symlog scale)'};
          end
        else
          plot(t, y_bottom, 'displayname',simulationname, 'linestyle', localLS); hold on;
        end
        ylabel(ylabR);
%         set(gca,'ycolor','k');
        % link axes t
        linkaxes(axxx, 'x');
        
      case {0,1,2}
        % if classical plot
        if(logscale && peak2peak(y)>0)
          semilogy(t, y, 'displayname',simulationname, 'linestyle', localLS); hold on;
        else
          plot(t, y, 'displayname',simulationname, 'linestyle', localLS); hold on;
        end
      
      otherwise
        error('swtchE_0K_1P_2T_4sbpltKP has wrong value');
    end
  end
  
  xlabel('time $t$ [s]');
  switch(swtchE_0K_1P_2T_4sbpltKP)
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
%     if(swtchE_0K_1P_2T_4sbpltKP~=3)
%       ylim(niceLogYLim([min(y(sel)),max(y(sel))]));
%     end
%   end

  if(swtchE_0K_1P_2T_4sbpltKP==4)
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