% Author:        Léo Martire.
% Description:   Plots effective sound speed of an atmospheric model file.
% Notes:         TODO.
% Usage:
%   plot_model_effective_soundspeed(DATAFILE, upwind1_downwind0)
% with:
%   TODO.
% yields:
%   TODO.
%
% /!\ /!\ /!\ SOMEWHAT DEPRECATED /!\ /!\ /!\

function [] = plot_model_ceff_one(DATAFILE, upwind1_downwind0)
  [Z, ~, ~, C, ~, ~, ...
   ~, ~, ~, ~, ~, ~, ~, W, ~, ~, ~] = ...
   extract_atmos_model(DATAFILE, 3, 0, 0);

%   [~, posstr, ~, d, s, ~, ~, F107A, F107, AP] = extract_atmos_model_setup(output_file);
  [datestr, posstr] = extract_atmos_model_setup(DATAFILE);
%   apf107str=strcat("F10.7 avg. ", sprintf("%.1f",F107A), ", F10.7 ", sprintf("%.1f",F107), ", AP ", sprintf("%.1f",AP));
%       datestr=[num2str(y),', ',num2str(d), 'th day, ', num2str(floor(s/3600)),':',num2str(floor((s - floor(s/3600)*3600)/60)), ' UT'];
%   datestr=[num2str(d), 'th day, ', num2str(floor(s/3600)),':',num2str(floor((s - floor(s/3600)*3600)/60)), ' UT'];
  tit_plus={posstr, datestr};
  
  figure();
  if(upwind1_downwind0~=0)
    semilogx(C, Z, ":", min(C+W,C-W), Z, "-", max(C+W,C-W), Z, "-");
    legend('sound speed $c$', 'upwind effective sound speed', 'downwind effective sound speed', 'Location', 'best');
  else
    semilogx(C, Z, ":", C+W, Z, "-", C-W, Z, "-");
    legend('sound speed $c$', 'rightward effective sound speed $c+w$', 'leftward effective sound speed $c-w$', 'Location', 'best');
  end
  %xlim([0.5 * min([east_Richardson; north_Richardson; proj_Richardson]), 2 * max([east_Richardson; north_Richardson; proj_Richardson])]);
  xlabel('sound speed (m/s)'); ylabel('altitude (m)');
  title(tit_plus);
end

