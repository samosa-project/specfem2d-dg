function [ext, chantxt] = getUnknowns(type_displ, chan)
  % Switch on type of display.
  switch type_displ
%   if (type_display == 1) % Original SPECFEM2D's synthetic is displacement.
    case 1 % Original SPECFEM2D's synthetic is displacement.
      ext = 'semd'; % Because original SPECFEM2D's synthetic is displacement.
      % For stations in solid zones it's displacement. For stations in DG zones it's velocity.
      if (strcmp(chan, 'BXZ'))
        chantxt = 'vertical {$u_z$ (m), $v_z$ [m/s]}';
      elseif (strcmp(chan, 'BXX'))
        chantxt = 'horizontal {$u_x$ (m), $v_x$ [m/s]}';
      else
        error(['[',mfilename,', ERROR] The variable ''unknown'' has a non-standard value.']);
      end
%   elseif (type_display == 2) % Original SPECFEM2D's synthetic is velocity.
    case 2 % Original SPECFEM2D's synthetic is velocity.
      ext = 'semv'; % Because original SPECFEM2D's synthetic is velocity.
      if (strcmp(chan, 'BXZ'))
        chantxt = '{$v_z$ [m/s], $\delta P$ [Pa]}';
      elseif (strcmp(chan, 'BXX'))
        chantxt = '{$v_x$ [m/s], $\delta P$ [Pa]}';
      else
        error(['[',mfilename,', ERROR] The variable ''unknown'' has a non-standard value.']);
      end
%   elseif (type_display == 4) % Original SPECFEM2D's synthetic is pressure.
    case 4 % Original SPECFEM2D's synthetic is pressure.
      ext = 'semp'; % Because original SPECFEM2D's synthetic is pressure.
      if (strcmp(chan, 'PRE'))
        chantxt = '$\delta P$ [Pa]';
      else
        error(['[',mfilename,', ERROR] The variable ''unknown'' has a non-standard value.']);
      end
    otherwise
      error(['[',mfilename,', ERROR] This type_display (',num2str(type_displ),') is not implemented for this script.']);
  end
end