for i = 1:numel(cases)
  if(cases{i}.fts0_stf1)
    if(cases{i}.ortho0_slant1)
      istfslant = i;
    else
      istfortho = i;
    end
  else
    if(cases{i}.ortho0_slant1)
      iftsslant = i;
    else
      iftsortho = i;
    end
  end
end