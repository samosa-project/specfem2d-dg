clear all;
close all;
clc;

thisFolder = regexprep(mfilename('fullpath'),mfilename,''); cd(thisFolder);
name_fig = ['snapshot_mosaic'];

% t = [460, 754, 1048, 1342]+800;
% IDs = t/1.4e-3;
root = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_transmission/'];
OFDs = {[root, 'OUTPUT_FILES_LNS_F2S_plotvz/'], [root, 'OUTPUT_FILES_LNS_S2F_plotvz/']};
IDs = [100, 700, 900];

parfile = [thisFolder, 'parfile_input'];
dt = readExampleFiles_extractParam(parfile, 'DT', 'float');

t = []; c = 1;
for j = 1:numel(OFDs)
  for i = 1:numel(IDs)
    id = IDs(i);
    t(c) = id * dt;
    snap = [OFDs{j}, 'image',sprintf('%07d', id),'.jpg'];

    img{c} = imread(snap);
    NX{c} = size(img, 2);
    NZ{c} = size(img, 1);

    img{c} = img{c}(:,1:1085,:); % crop right
    img{c} = img{c}(:,100:end,:); % crop left
    img{c} = img{c}(1:885,:,:); % crop bottom
    img{c} = img{c}(135:end,:,:); % crop top
    c = c+1;
  end
end

% Plot.
fh = figure('units','normalized','outerposition',[0,0,1,1]);
tightAxes = tight_subplot(2, 3, [0.075, 0.01], [0.01, 0.01], [1, 1]*0.02);

for i = 1:numel(img)
  axes(tightAxes(i));
  TIT = ['$t=',sprintf('%.0f', 1e3*t(i)),'$~ms'];
  image(img{i});
  daspect([1,1,1]);
  title(TIT);
  xticks({}); yticks({});
end

ll = add_labels_subplots(gcf,0.9,0,[0,-0.07]);
% set(tightAxes, 'xtick', {}, 'ytick', {});

figpath = [thisFolder, name_fig];
customSaveFig(fh, figpath, {'fig', 'eps'}, 9999);

% move produced figures to thesis
disp(['[] Starting to move Figures to thesis folder.']);
system(['cp ', figpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/images_verif_transmission/']);
