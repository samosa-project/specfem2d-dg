max_relative_err_all_stats = []; % save maximum error
t_forErr = {};
synth_forErr = {};
anal_forErr = {};
err_v = {};
for i = 1:numel(allStations)
  % Produce difference.
  if(min(diff(Ztime(1,:))) < min(diff(t)))
    % if dt_sythetic < dt_analytic, interpolate analytic on synthetic t (Ztime)
    t_forErr{i} = Ztime(i,1:nt);
    synth_forErr{i} = Zamp(i,1:nt);
    anal_forErr{i} = interp1(t, synf(i,:), t_forErr{i});
  else
    % if dt_sythetic < dt_analytic, interpolate synthetic on analytic t (t)
    t_forErr{i} = t;
    synth_forErr{i} = interp1(Ztime(i,1:nt), Zamp(i,1:nt), t_forErr{i});
    anal_forErr{i} = synf(i,:);
  end
  err_v{i} = factor_err * abs(synth_forErr{i} - anal_forErr{i}); max_relative_err_all_stats(i) = (max(err_v{i})/factor_err)/peak2peak(synf(i,:));
  
  % Seek a possible lag.
  synth_forErr_forXCorr = synth_forErr{i};
  synth_forErr_forXCorr(isnan(synth_forErr{i}))=0; % remove nans for xcorr
  anal_forErr_forXCorr = anal_forErr{i};
  anal_forErr_forXCorr(isnan(anal_forErr{i}))=0; % remove nans for xcorr
  [xc, xcl] = xcorr(synth_forErr_forXCorr, anal_forErr_forXCorr);
  lag = xcl(find(abs(xc) == max(abs(xc))));
  disp(['[',mfilename,'] Lag for station ',num2str(i),' = ',num2str(lag),' indices. Shifing.']);
  anal_forErr{i}(max(1+lag,1):end+min(lag,0)) = anal_forErr{i}(max(1-lag,1):end+min(-lag,0));
  if(lag<0)
    anal_forErr{i}(end+lag:end) = nan;
  else
    anal_forErr{i}(1:1+lag) = nan;
  end
  err_v{i} = factor_err * abs(synth_forErr{i} - anal_forErr{i});
  max_relative_err_all_stats(i) = (max(err_v{i})/factor_err)/peak2peak(synf(i,:));
end