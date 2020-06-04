function [tstack, pselsave, pstacc] = get_staccs(Tda, Pda, tstacctime, kind, fcut, order)

  tstack = {};
  pstacc = {};
  pselsave = {};
  
  % design a singular filter for all time series
  switch(kind)
    case 'bp'
      theFilter = designfilt('bandpassiir', 'FilterOrder', order, ...
                             'HalfPowerFrequency1', min(fcut), ...
                             'HalfPowerFrequency2', max(fcut), ...
                             'SampleRate', 1/mean(diff(Tda(1,:))));
  end
  
  for i = 1:size(Pda, 1)
    numtstac = numel(tstacctime{i});
    for ti=1:numtstac
      twin = tstacctime{i}(ti)+[-1e-2, 0.1];
      if(ti==1)
        sel = (Tda(i, :)>=min(twin)) & (Tda(i, :)<=max(twin));
        num_ids = sum(sel);
      end
      ID_begin = find(Tda(i, :)>=min(twin), 1, 'first');
      sel = ID_begin + [0:num_ids];
      if(ti==1)
        pstacc{i} = zeros(1, numel(sel));
        pselsave{i} = zeros(1, numel(sel));
      end
    %   pselsave{i}(ti, :) = Pda(i, sel);
%       pselsave{i}(ti, :) = quick_filter(Tda(i, sel), Pda(i, sel), kind, fcut, order);
      pselsave{i}(ti, :) = filtfilt(theFilter, Pda(i, sel));
      pstacc{i} = pstacc{i} + pselsave{i}(ti, :);
      if(ti==1)
        tstack{i} = Tda(i, sel);
      end
    end
    pstacc{i}=pstacc{i}/numtstac;
  end
end