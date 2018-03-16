function [] = plot_model_effective_soundspeed(output_file)
  [Z, ~, ~, C, ~, ~, ...
   ~, ~, ~, ~, ~, ~, ~, W, ~, ~, ~] = ...
   extract_data([output_file], 3, 0, 0);

  [~, posstr, ~, d, s, ~, ~, F107A, F107, AP] = extract_data_setup(output_file);
  apf107str=strcat("F10.7 avg. ", sprintf("%.1f",F107A), ", F10.7 ", sprintf("%.1f",F107), ", AP ", sprintf("%.1f",AP));
%       datestr=[num2str(y),', ',num2str(d), 'th day, ', num2str(floor(s/3600)),':',num2str(floor((s - floor(s/3600)*3600)/60)), ' UT'];
  datestr=[num2str(d), 'th day, ', num2str(floor(s/3600)),':',num2str(floor((s - floor(s/3600)*3600)/60)), ' UT'];
  tit_plus={posstr, datestr};
  
  figure();
  semilogx(C, Z, ":", C+W, Z, "-", C-W, Z, "-");
  %xlim([0.5 * min([east_Richardson; north_Richardson; proj_Richardson]), 2 * max([east_Richardson; north_Richardson; proj_Richardson])]);
  xlabel('sound speed (m/s)'); ylabel('altitude (m)');
  legend('sound speed $c$', 'rightward effective sound speed $c+w$', 'leftward effective sound speed $c-w$', 'Location', 'best');
  title(tit_plus);
end

