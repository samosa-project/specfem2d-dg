function customSaveFig(path, extensions)
  if(not(exist('extensions')))
    extensions={'jpg', 'fig'};
  else
    for i=1:numel(extensions)
      switch(extensions{i})
        case {'jpg', 'fig', 'png', 'eps'}
          % nothing
        otherwise
          error(['[',mfilename,', ERROR] Extension ''',extensions{i},''' not implemented.']);
      end
    end
  end
  
  fig = gcf;
  fig.InvertHardcopy = 'off';
  figsavecol = fig.Color;
  
  for i = 1:numel(extensions)
    ee = extensions{i};
    eee = ee;
    if(strcmp(ee,'eps'))
      eee = 'epsc';
    end
    
    fullpath=[path, '.', ee];
    
    if(strcmp(ee,'png') || strcmp(ee,'eps'))
      set(fig, 'color', 'none');
    elseif(strcmp(ee,'jpg') || strcmp(ee,'fig'))
      set(fig, 'color', 'w');
    else
      error('kek');
    end
    
    saveas(fig, fullpath, eee);
    
    disp(['[', mfilename, '] Saved figure to ''', fullpath, '''.']);
  end
  
  set(fig, 'color', figsavecol);
end