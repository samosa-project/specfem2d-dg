% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [] = plot_model_effective_soundspeed(output_file, upwinddownwind)
  [Z, ~, ~, C, ~, ~, ...
   ~, ~, ~, ~, ~, ~, ~, W, ~, ~, ~] = ...
   extract_data(output_file, 3, 0, 0);

  [~, posstr, ~, d, s, ~, ~, F107A, F107, AP] = extract_data_setup(output_file);
%   apf107str=strcat("F10.7 avg. ", sprintf("%.1f",F107A), ", F10.7 ", sprintf("%.1f",F107), ", AP ", sprintf("%.1f",AP));
%       datestr=[num2str(y),', ',num2str(d), 'th day, ', num2str(floor(s/3600)),':',num2str(floor((s - floor(s/3600)*3600)/60)), ' UT'];
  datestr=[num2str(d), 'th day, ', num2str(floor(s/3600)),':',num2str(floor((s - floor(s/3600)*3600)/60)), ' UT'];
  tit_plus={posstr, datestr};
  
  figure();
  if(upwinddownwind~=0)
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

