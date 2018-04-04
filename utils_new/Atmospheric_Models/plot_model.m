% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [] = plot_model(output_file, subplotcode, marker, colour, legends)
  [Z, RHO, ~, C, ~, ~, ...
   ~, NSQ, ~, ~, ~, WN, WE, ~, ~, ~, GAMMA] = ...
   extract_data(output_file, 3, 0, 0);

  [~, posstr, ~, d, s, ~, ~, F107A, F107, AP] = extract_data_setup(output_file);
  apf107str=strcat("F10.7 avg. ", sprintf("%.1f",F107A), ", F10.7 ", sprintf("%.1f",F107), ", AP ", sprintf("%.1f",AP));
%       datestr=[num2str(y),', ',num2str(d), 'th day, ', num2str(floor(s/3600)),':',num2str(floor((s - floor(s/3600)*3600)/60)), ' UT'];
  datestr=[num2str(d), 'th day, ', num2str(floor(s/3600)),':',num2str(floor((s - floor(s/3600)*3600)/60)), ' UT'];

  figure();
  subplot(subplotcode+1); myplot(RHO, Z, marker, colour, datestr, 'sx'); xlabel('$\rho$ (kg/m$^3$)'); % xlim([0.1*min(RHO),10*max(RHO)]);
  subplot(subplotcode+2); myplot(C, Z, marker, colour, datestr, 'p'); xlabel('$c$ (m/s)'); % LIMX=C; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]);
  title({[posstr],apf107str,''});
  subplot(subplotcode+3); myplot(NSQ, Z, marker, colour, datestr, 'p'); xlabel('$N^2$ (rad$^2$/s$^2$)'); % LIMX=NSQ; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]);
  if(strcmp(legends,'no')==0)
    legend('Location', 'northwest');
  end
  subplot(subplotcode+4); myplot(WN, Z, marker, colour, datestr, 'p'); xlabel('$w_M$ (m/s)'); % LIMX=WN; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]);
  subplot(subplotcode+5); myplot(WE, Z, marker, colour, datestr, 'p'); xlabel('$w_Z$ (m/s)'); % LIMX=WE; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]);
  subplot(subplotcode+6); myplot(GAMMA, Z, marker, colour, datestr, 'p'); xlabel('$\gamma$'); % LIMX=GAMMA; xlim([1.1*min(LIMX)-0.1*max(LIMX),1.1*max(LIMX)-0.1*min(LIMX)]); legend('Location', 'best');
end

function myplot(XDATA, YDATA, MARKER, COLOUR, DISPLAYNAME, TYPE)
  if(strcmp(TYPE,'p'))
    plot(XDATA, YDATA, MARKER, 'Color', COLOUR, 'DisplayName', DISPLAYNAME); hold on;
  elseif(strcmp(TYPE,'sx'))
    semilogx(XDATA, YDATA, MARKER, 'Color', COLOUR, 'DisplayName', DISPLAYNAME); hold on;
  end
  ylim([min(YDATA), max(YDATA)]);
end