function [rhof, c, rhos, vp, vs] = get_models(parfile)
  if(not(strcmp(readExampleFiles_extractParam(parfile, 'MODEL', 'string'),'default')) || readExampleFiles_extractParam(parfile, 'USE_ISOTHERMAL_MODEL', 'bool'))
    error(['[',mfilename,', ERROR] Should only use isobaric models.']);
  end
  c = readExampleFiles_extractParam(parfile, 'sound_velocity', 'float');
  rhof = readExampleFiles_extractParam(parfile, 'surface_density', 'float');
  [~, solid_model] = readExampleFiles_extractParfileModels(parfile);
  rhos = solid_model(2).rho;
  vp = solid_model(2).vp; 
  vs = solid_model(2).vs;
end

