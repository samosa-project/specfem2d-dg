function factor = computeScalings(stat_number, geomAtt, x_stat, z_stat, d_stat, rescale_fact, askUserToConfirm)
  if(not(exist('askUserToConfirm', 'var')))
    askUserToConfirm = 1;
  end
  % Renormalisation (global, and geometric).
  factor = 1; % Reset to default value for each station.
  if (geomAtt ~= 0)
    % If geometric_attenuation was asked by user.
    switch geomAtt
      case 1
        geom_att_fact = d_stat(stat_number) ^ 0.5; % Geometric attenuation respective to raw distance to source.
      case 2
        geom_att_fact = abs(x_stat(stat_number)) ^ 0.5; % Geometric attenuation respective to horizontal distance to source.
      case 3
        geom_att_fact = abs(z_stat(stat_number)) ^ 0.5; % Geometric attenuation respective to vertical distance to source.
      otherwise
        error(['[',mfilename,', ERROR] Geometric attenuation parameter not implemented.']);
    end
    if (geom_att_fact == 0)
      % If exactly at zero distance, consider no geometric rescaling.
      geom_att_fact = 1;
    end
    factor = factor / geom_att_fact;
  end
  renorm_statbystat = - 1;
  if(rescale_fact~=1 && askUserToConfirm==1)
    % Rescaling was asked. Check again with user.
    renorm_statbystat = - 1;
    disp(['[',mfilename,'] Specified rescale factor is ', num2str(rescale_fact), '.']);
    inputtxt = ['[',mfilename,'] Rescale data for station ', num2str(stat_number), '? (0 for no, 1 for yes) > '];
    while (not(ismember(renorm_statbystat, [0, 1])))
      renorm_statbystat = input(inputtxt);
      if (isempty(renorm_statbystat))
        renorm_statbystat = 1;
      end
    end
  else
    % Confirmation overriden, do apply factor.
    renorm_statbystat = 1;
  end
  if (renorm_statbystat == 1)
    % If rescaling is actually wanted by user, do it.
    factor = factor * rescale_fact;
  end
end