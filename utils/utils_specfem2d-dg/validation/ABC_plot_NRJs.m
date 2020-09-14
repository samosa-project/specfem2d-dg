% custom function based off plot_total_energy, but for nicer ABC validation
% visualisations

function [figureHandle, tight_axes] = ABC_plot_NRJs(OFDIRS_struct, DISPLAYNAMES_struct, LS_strct, colours_runs, outputfigpath)
  
  % OFDIRS should be a struct containing the test cases
  % OFD_LIST should be a length 9 cell: {[3 PW], [3 PS], [3 WPS]}.
  i = 1;
  OFD_cell{i, 1} = OFDIRS_struct.PW.LAR; OFD_cell{i, 2} = OFDIRS_struct.PW.FAF; OFD_cell{i, 3} = OFDIRS_struct.PW.BUF; i = i + 1;
  OFD_cell{i, 1} = OFDIRS_struct.PS.LAR; OFD_cell{i, 2} = OFDIRS_struct.PS.FAF; OFD_cell{i, 3} = OFDIRS_struct.PS.BUF; i = i + 1;
  OFD_cell{i, 1} = OFDIRS_struct.WPS.LAR; OFD_cell{i, 2} = OFDIRS_struct.WPS.FAF; OFD_cell{i, 3} = OFDIRS_struct.WPS.BUF;
  DISPLAYNAMES_cell = {DISPLAYNAMES_struct.LAR, DISPLAYNAMES_struct.FAF, DISPLAYNAMES_struct.BUF};
  LS_cell = {LS_strct.LAR, LS_strct.FAF, LS_strct.BUF};
  cases_titles = {'PW', 'PS', 'WPS'};
  num_cases = 3;
  
  YLAB_top = {['normalised'], ['potential energy'], ['perturbation']};
  YLAB_bot = {['normalised'], ['kinetic energy'], ['perturbation']};
  XLAB = ['time [s]'];
  
  % load all
  t = {};
  ke = {};
  pe = {};
  te = {};
  for test_case = 1:num_cases
    for run = 1:3
      OFD = OFD_cell{test_case, run};
      [t{test_case, run}, ke{test_case, run}, pe{test_case, run}, te{test_case, run}] = load_energy(OFD);
      
      sel = (t{test_case, run}>=0);
      t{test_case, run} = t{test_case, run}(sel);
      pe{test_case, run} = pe{test_case, run}(sel);
      ke{test_case, run} = ke{test_case, run}(sel);
      te{test_case, run} = te{test_case, run}(sel);
      pe{test_case, run} = pe{test_case, run}/max(pe{test_case, run});
      ke{test_case, run} = ke{test_case, run}/max(ke{test_case, run});
      te{test_case, run} = te{test_case, run}/max(te{test_case, run});
    end
  end
  
  % do the fig
  figureHandle = figure('units','normalized','outerposition', [0 0 1 0.75]);
  gap_h = 0.01; gap_w = 0.01;
  marg_h_low = 0.14; marg_h_upp = 0.08;
  marg_w_lef = 0.11; marg_w_rig = 0.01;
  tight_axes = tight_subplot(2, 3, [gap_h, gap_w], [marg_h_low, marg_h_upp], [marg_w_lef, marg_w_rig]); % not mandatory, but prettier
  for i = 1:num_cases
    % Plot potential energy.
    axes(tight_axes(i));
    cur_y = pe;
    for j = 1:3
      plot(t{i,j}, cur_y{i,j}, 'color', colours_runs{j}, 'displayname', DISPLAYNAMES_cell{j}, 'linestyle', LS_cell{j}); hold on;
    end
%     if(min(cur_y{i,1})>0 && min(cur_y{i,2})>0 && min(cur_y{i,3})>0)
%       set(gca,'yscale','log');
%     end
    if(i==1)
      ylab_top_h = ylabel(YLAB_top);
      legend('location', 'northeast');
    else
      yticklabels({});
    end
    xticklabels({});
    title(cases_titles{i});
    
    % Plot kinetic energy.
    axes(tight_axes(i+num_cases));
    cur_y = ke;
    for j = 1:3
      plot(t{i,j}, cur_y{i,j}, 'color', colours_runs{j}, 'displayname', DISPLAYNAMES_cell{j}, 'linestyle', LS_cell{j}); hold on;
    end
%     if(min(cur_y{i,1})>0 && min(cur_y{i,2})>0 && min(cur_y{i,3})>0)
%       set(gca,'yscale','log');
%     end
    if(i==1)
      ylab_bot_h = ylabel(YLAB_bot);
    else
      yticklabels({});
    end
    xlabel(XLAB);
  end
  
  linkaxes(tight_axes, 'x'); linkprop(tight_axes, 'XTick');
  custom_xticks = linspace(0,0.6,4); custom_xticks(end) = [];
  set(tight_axes(1), 'xlim', [0, max(t{1,1})], 'xtick', custom_xticks);
  
%   maxabsval = 1.15;
%   ylimmm = [-1, 1]*maxabsval; yticksss = 'auto';
%   ylimmm = [-1.5, 1.15];
  yticksss = [-2:0.5:2];
%   % linked y, but separate top/bot
%   linkaxes(axxx(1:num_cases), 'y'); linkprop(axxx(1:num_cases), 'YTick');
%   set(axxx(1), 'ylim', [-1, 1]*maxabsval);
%   linkaxes(axxx(num_cases + (1:num_cases)), 'y'); linkprop(axxx(num_cases + (1:num_cases)), 'YTick');
%   set(axxx(num_cases+1), 'ylim', [-1, 1]*maxabsval);
  % linked y, all
  linkaxes(tight_axes(1:3), 'y'); linkprop(tight_axes(1:3), 'YTick', 'Yscale');
  linkaxes(tight_axes(4:6), 'y'); linkprop(tight_axes(4:6), 'YTick', 'Yscale');
  set(tight_axes, 'xtick', 0:0.1:0.5);
  set(tight_axes(1:3), 'ylim', [-1, 1]*1.15, 'ytick', yticksss, 'yscale', get(tight_axes(1), 'yscale'));
  set(tight_axes(4:6), 'ylim', [-1.5, 1.15], 'ytick', yticksss, 'yscale', get(tight_axes(1), 'yscale'));
  
  ylab_top_h.Units = 'normalized';
  ylab_bot_h.Units = ylab_top_h.Units;
  min_ylab_h_pos = min([ylab_top_h.Position(1), ylab_bot_h.Position(1)]);
  set(ylab_top_h, 'Position', [min_ylab_h_pos, ylab_top_h.Position(2:3)]);
  set(ylab_bot_h, 'Position', [min_ylab_h_pos, ylab_bot_h.Position(2:3)]);
  
  ll = add_labels_subplots(figureHandle, 0.85);
  for i=1:numel(ll);set(ll(i), 'Position', get(ll(i), 'Position')+[0,-0.3225,0,0]);end
  
  customSaveFig(figureHandle,outputfigpath,{'fig', 'eps'},9999);
end