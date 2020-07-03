function [rho_1, vp_1, Z_1, rho_2, vp_2, vs_2, Z_2P, Z_2S] = grab_models(parfile)
  rho_1 = readExampleFiles_extractParam(parfile,'surface_density','float');
  vp_1 = readExampleFiles_extractParam(parfile,'sound_velocity','float');
  Z_1 = rho_1*vp_1;
  [~, mas]=readExampleFiles_extractParfileModels(parfile);
  if(not(numel(mas)==2))
    error(['Must use exactly 2 models (one for air and one for solid).']);
  end
  mas(readExampleFiles_extractParam(parfile,'id_region_DG','int'))=[]; % remove dg model
  rho_2 = mas.rho;
  vp_2 = mas.vp;
  vs_2 = mas.vs;
  Z_2P = rho_2*vp_2;
  Z_2S = rho_2*vs_2;
end

