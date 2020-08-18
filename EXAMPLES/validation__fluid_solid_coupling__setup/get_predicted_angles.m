function [i1_i2_j2t_j2r] = get_predicted_angles(parfile, ic)
  if(not(strcmp(readExampleFiles_extractParam(parfile, 'MODEL', 'string'),'default')) || readExampleFiles_extractParam(parfile, 'USE_ISOTHERMAL_MODEL', 'bool'))
    error(['[',mfilename,', ERROR] Should only use isobaric models.']);
  end
  c = readExampleFiles_extractParam(parfile, 'sound_velocity', 'float');
  [~, solid_model] = readExampleFiles_extractParfileModels(parfile);
  vp = solid_model(2).vp; 
  vs = solid_model(2).vs;
  i1_i2_j2t_j2r = asin([c*sin(ic)/vp, vp*sin(ic)/c, vs*sin(ic)/c, vs*sin(ic)/vp])* 180/pi;
end

