function [s2f1_or_f2s0, fig_title, f0] = check_test_case(parfile, sourcefile)
  zinterface = readExampleFiles_extractParam(parfile,'coord_interface','float');
  zs = readExampleFiles_extractParam(sourcefile,'zs','float');
  if(numel(zs)>1)
    % If many sources found (in the case of S2F coupling), check if all same altitude and keep only one.
    if(numel(unique(zs))==1)
      zs = zs(1);
    else
      error(['[',mfilename,'] Sources must be at same altitude.']);
    end
  end
  f0 = readExampleFiles_extractParam(sourcefile,'f0','float');
  if(numel(f0)>1)
    % If many f0 found (in the case of S2F coupling), check if all same and keep only one.
    if(numel(unique(f0))==1)
      f0 = f0(1);
    else
      error(['[',mfilename,'] Sources must have same frequency.']);
    end
  end
  % Distinguish type of coupling based on source altitude.
  if(zs>zinterface)
    % source is above ground, coupling is thus fluid to solid
    s2f1_or_f2s0 = 0;
    fig_title = ['F2S'];
  else
    if(zs==zinterface)
      error(['source must not be on interface']);
    end
    % source is above ground, coupling is thus solid to fluid
    s2f1_or_f2s0 = 1;
    fig_title = ['S2F'];
  end
end