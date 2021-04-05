ic_rad = 18*pi/180; % incidence angle for slanted case [rad]
f0 = 100; % [Hz]

thisFolder = [regexprep(mfilename('fullpath'),[filesep,mfilename],filesep)];
plotFolder = [thisFolder,filesep,'plots'];
extToSave = {'jpg', 'eps'};

SPCFMEXLOC = [thisFolder,filesep,'..',filesep];
basename = 'validation__fluid_solid_coupling__';
tagfts = 'fts';
tagstf = 'stf';
tagortho = 'ortho';
tagslant = 'slant';

ic_deg = ic_rad*180/pi;

cases = {};
i = 1;
cases{i}.fts0_stf1 = 1; cases{i}.ortho0_slant1 = 1; cases{i}.tlim = [0.035,0.095]; cases{i}.code='SF_slant'; i = i+1;
cases{i}.fts0_stf1 = 1; cases{i}.ortho0_slant1 = 0; cases{i}.tlim = [0.035, 0.095]; cases{i}.code='SF_ortho'; i = i+1;
cases{i}.fts0_stf1 = 0; cases{i}.ortho0_slant1 = 1; cases{i}.tlim = [0.08, 0.18]; cases{i}.code='FS_slant'; i = i+1;
cases{i}.fts0_stf1 = 0; cases{i}.ortho0_slant1 = 0; cases{i}.tlim = [0.08, 0.18]; cases{i}.code='FS_ortho'; i = i+1;

nvals = [50, 100, 500];

stftag = 'S2F';
ftstag = 'F2S';

folderz = {}; nelz = {}; idcasez = {};
c = 1;
for j = 1:numel(nvals)
  curnel = nvals(j);
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
    folder = [folder, '__',sprintf('%04d',curnel)];
    folder = [SPCFMEXLOC, folder, filesep];
    if(not(exist(folder, 'dir')))
      mkdir(folder);
    end
    folderz{c} = folder;
    nelz{c} = curnel;
    idcasez{c} = i;
    c = c+1;
  end
end

if(not(exist(plotFolder,'dir')))
  mkdir(plotFolder);
end