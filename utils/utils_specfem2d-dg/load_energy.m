function [t, ke, pe, te] = load_energy(OFDIR)
  energyfilename='energy.dat';
  energypath = [OFDIR,filesep,energyfilename];
  if(not(exist(energypath)==2))
    disp(['[',mfilename,', INFO] Energy file ''',energypath,'''does not exist. Skipping.']);
    return;
  end
  EF = importdata(energypath);
  t  = EF.data(:, 1);
  ke = EF.data(:, 2);
  pe = EF.data(:, 3);
  te = EF.data(:, 4);
end

