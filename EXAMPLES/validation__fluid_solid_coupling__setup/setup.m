ic = 30*pi/180; % incidence angle for slanted case [rad]
f0 = 100; % [Hz]

SPCFMEXLOC = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES',filesep];
basename = 'validation__fluid_solid_coupling__';
tagfts = 'fts';
tagstf = 'stf';
tagortho = 'ortho';
tagslant = 'slant';

cases = {};
i = 1;
cases{i}.fts0_stf1 = 0; cases{i}.ortho0_slant1 = 0; cases{i}.tlim = [0.1, 0.2]; i = i+1;
cases{i}.fts0_stf1 = 0; cases{i}.ortho0_slant1 = 1; cases{i}.tlim = [0.1, 0.2]; i = i+1;
cases{i}.fts0_stf1 = 1; cases{i}.ortho0_slant1 = 0; cases{i}.tlim = [0.04, 0.14]; i = i+1;
cases{i}.fts0_stf1 = 1; cases{i}.ortho0_slant1 = 1; cases{i}.tlim = [0.04, 0.14]; i = i+1;

folderz = {};
for i = 1:numel(cases)
  % Define working folders.
  folder = [basename];
  if(cases{i}.fts0_stf1)
    folder = [folder, tagstf];
  else
    folder = [folder, tagfts];
  end
  folder = [folder, '__'];
  if(cases{i}.ortho0_slant1)
    folder = [folder, tagslant];
  else
    folder = [folder, tagortho];
  end
%   folder
  folder = [SPCFMEXLOC, folder, filesep];
  if(not(exist(folder, 'dir')))
    mkdir(folder);
  end
  folderz{i} = folder;
end