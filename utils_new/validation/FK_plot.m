% Prepare figure name.
if(externalDGAtmosModel)
  savefigname = [savefigname, '_externalDG'];
else
  if(USE_ISOTHERMAL_MODEL)
    savefigname = [savefigname, '_isothermal'];
  else
    savefigname = [savefigname, '_isobaric'];
  end
end
spl = split(OFd,'/');
savefigname = [savefigname, regexprep(spl{end-2}, prefix, '')];

% Actually plot.
fh = figure('units','normalized','outerposition',[0 0 1 1]);
titt_nlines = 2;
axxxx = tight_subplot(nstatToPlot, 1, [0, -1], [0.1,0.02+(titt_nlines*0.11/3)], [0.07, 0.02]); % not necessary, but prettier
for j = 1:numel(stationsToPlot)
  axes(axxxx(nstatToPlot-j+1)); hold on;
  i = stationsToPlot(j);
  plot(t_forErr{i}, anal_forErr{i}, 'Color', [0 0 0], 'LineWidth', 1, 'displayname', 'analytical'); % Plot analytic value.
  plot(t_forErr{i}, synth_forErr{i}, '-.k', 'LineWidth', 2, 'displayname', 'synthetic'); % Plot synthetic.
  ylim_b4_err = get(gca,'ylim'); ylim_b4_err=max(abs(ylim_b4_err))*1.1; pow10=10^(-floor(log10(max(ylim_b4_err)))); ylim_b4_err = ceil(pow10*ylim_b4_err)/pow10*[-1,1]; % save and round ylim
  plot(t_forErr{i},err_v{i},':k','LineWidth',1.5,'displayname',['$',num2str(factor_err),'{\times}|$anal.$-$synth.$|$']); % Plot difference.
  ylim([ylim_b4_err]); % Redo ylim.
  if (j == round(nstatToPlot/2))
    ylabel(['$',vname_long,'_z$ ',vunit]);
  end
  yticklabels('auto');
  if (j == 1)
    xticklabels('auto');
    xlabel('time [s]');
  else
    xticklabels([]);
  end
  if(j==nstatToPlot)
    legend('location','east');
  end
end
linkaxes(axxxx, 'x');
if(not(USE_ISOTHERMAL_MODEL))
  linkaxes(axxxx, 'xy');
  curylim = ylim; curylim=max(abs(curylim))*1.1; pow10=10^(-floor(log10(max(curylim)))); ylim(ceil(pow10*curylim)/pow10*[-1,1]); % round ylim
end
xlim([max(min(min(Ztime)),min(t)), min(max(max(Ztime)),max(t))]);
% title(['Gravito-acoustic wave propagation. Stations at altitude : ', num2str(zstattab(locStatNumber)),'km (along z)']);
title({['Stations'' Altitudes = [', sprintf('%.1f ',zstattab(istattab)/1000),'] [km]'],['Sound Speed = ',num2str(SOUNDSPEED),' [m/s]. Max. Relative Error = ',num2str(max(max_relative_err_all_stats*100)),' \%']});
% addpath('/usr/local/matlab/r2018a/toolbox/tightfig'); tightfig; % eventually tighten fig
prettyAxes(fh);
savefigfullpath = [savefigpath, regexprep(savefigname,'\.','')]; % regexprep because latex crashed when filenames have dots in it
% customSaveFig(savefigfullpath,{'fig','eps','jpg'});
customSaveFig(fh, savefigfullpath, {'fig', 'jpg'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure of the nstat horizontal components and synthetic signals against 
% time
% figure
% for istat = 1 : nstat
%     
%     ax(istat)=subplot(nstat,1,nstat-istat+1) ; 
%     hold on
%     
%     if(seismotype ==1)
%          plot(timef,real(synf_X(istat,:)),'Color',[0 0 0],'LineWidth',1)
% %        plot(timef,real(synf_X(istat,:))/max(real(synf_X(istat,:))),'Color',[0 0 0],'LineWidth',1)
% 
% %                plot(Ztime(istat,:),real(synf_X(istat,1:length(Xtime(istat,:))))/max(real(synf_X(istat,1:length(Xtime(istat,:))))),'Color',[0 0 0],'LineWidth',1)
%     else
%     plot(Xtime(istat,:),syn(istat,:),'Color',[0 0 0],'LineWidth',1)
%     end
%   
%      plot(Xtime(istat,1:nt),Xamp(istat,1:nt),'-.k','LineWidth',2)
% %    plot(Xtime(istat,:),Xamp(istat,:)/max(Xamp(istat,:)),'-.k','LineWidth',2)
%     
%     if (istat == round(nstat/2))
%         ylabel([variable,' x-axis (m)'] ,'FontSize',14.3)
%     end
%     if (istat == 1)
%         xlabel('time (s)','FontSize',14.3)
%     end
%     
%     text(00.943*max(Xtime(istat,:)),0.95,['X position : $',num2str((istat-1)*dx_station + x_0),'km$ (along x)']) 
%     
%     legend('Analytical','Modeled','10*(Mo-An)','Location','West')
% 
% end
% 
% linkaxes(ax,'x')
% % xlim([0 20000])
% 
% title(strcat('Gravito-acoustic wave propagation. Stations at altitude : ', num2str(z_station/1000),'km (along z)'),'FontSize',24)


%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory of semd files
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4/') ;
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_LDDRK/') ;
% directory2 = strcat('/home/garcia/SATGRAVI/Brissaud/test_solution_analytique_gravi/Results.semd/') ;
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/gravi_forcing_nowind/') ; 
%directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/forcing_gravi_1Diso_windcte_noatten/') ; 
%  directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/FG1_1Diso_wind_noatten/') ; 
% directory = strcat('/home/garcia/SATGRAVI/Brissaud/files_wind_quentin/FG1_1Diso_nowind_noatten/') ; 
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4_3D/');
% directory = strcat('/home/qbrissaud/Documents/Results/LAST_GRAVI_Roland/200PROCS/');
% directory = strcat('./');
% OFd = [rootd,'OUTPUT_FILES/'];
% directory = strcat('/home/qbrissaud/Documents/FD/FD_14/');
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4_3D/');
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4/') ; 
% directory = strcat('/home/qbrissaud/Documents/Results/GJI_PAPER/1Diso_rhovar_grav_noatten_wind_wx10_graviForc_RK4/');

% OFd = [rootd,'OUTPUT_FILES_long/'];
% OFd = [rootd,'OUTPUT_FILES_1904171716_redone_velocity/'];
% OFd = [rootd,'OUTPUT_FILES_1904171808_vel_isobaric/'];
% OFd = [rootd,'OUTPUT_FILES_1904171832_vel_isobaric_LNS/'];
% OFd = [rootd,'OUTPUT_FILES_isobaric_LNS_190603_st2_morestations_corrected']; % this is p', we want v'
% OFd = [rootd,'OUTPUT_FILES_isobaric_FNS_190603_st1'];



%%%%%%%%%%%%%%%%%%%%%
% Stations parameters
% nstat      = 7 ;             % Number of station.
% dx_station = 60000; %58333.3333333;   % x-distance between stations
% z_station  = 100250.0;        % Height along z of stations
% % x_0        = 150000.0; % first station along x
% x_0        = 600250.0; % first station along x
% istattab=[6:7:48]
% xstattab=[450250.0:50000.0:750250.0]