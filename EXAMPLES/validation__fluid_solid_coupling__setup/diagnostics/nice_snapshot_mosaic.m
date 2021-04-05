clear all;
close all;
clc;

setup;
setup_find_all_cases; % produces istfortho, istfslant, iftsslant, iftsortho
itodo = [iftsortho, iftsslant, istfortho, istfslant];

figpath = [plotFolder,filesep,'snapshot_mosaic'];

% t = [460, 754, 1048, 1342]+800;
% IDs = t/1.4e-3;
% root = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_transmission/'];
% folderz = {[root, 'OUTPUT_FILES_LNS_F2S_plotvz/'], [root, 'OUTPUT_FILES_LNS_S2F_plotvz/']};
IDsFTS = [500, 1600];
IDsSTF = [300, 1100];

% ncol=2; nlin=4;
% map_ita_to_itodo = [iftsortho, iftsslant, iftsortho, iftsslant, istfortho, istfslant, istfortho, istfslant];
% map_ita_to_IDs = [IDsFTS(1), IDsFTS(1), IDsFTS(2), IDsFTS(2), IDsSTF(1), IDsSTF(1), IDsSTF(2), IDsSTF(2)];

ncol=4; nlin=2;
map_ita_to_itodo = [iftsortho, iftsortho, iftsslant, iftsslant, istfortho, istfortho, istfslant, istfslant];
map_ita_to_IDs = [IDsFTS, IDsFTS, IDsSTF, IDsSTF];

npan=ncol*nlin;
parfile = [folderz{1}, 'OUTPUT_FILES', filesep, 'input_parfile'];
dt = readExampleFiles_extractParam(parfile, 'DT', 'float');

% t = []; c = 1;
for ita = 1:npan
%   for i = 1:numel(IDsFTS)
    id = map_ita_to_IDs(ita);
    t(ita) = id * dt;
    snap = [folderz{map_ita_to_itodo(ita)}, 'OUTPUT_FILES', filesep, 'image',sprintf('%07d', id),'.jpg'];

    img{ita} = imread(snap);
    NX{ita} = size(img, 2);
    NZ{ita} = size(img, 1);

    img{ita} = img{ita}(:,1:787,:); % crop right
    img{ita} = img{ita}(:,45:end,:); % crop left
    img{ita} = img{ita}(1:750,:,:); % crop bottom
    img{ita} = img{ita}(80:end,:,:); % crop top
%     c = c+1;
%   end
end

% Plot.
fh = figure('units','normalized','outerposition',[0,0,1,1]);
gaph = 0.025;
tightAxes = tight_subplot(nlin, ncol, [0.05, gaph], [0.05, 0.02], [0.03, 0.005]);
fs=26;

XTCK = [1, 247, 371, 495, size(img{end},2)]; a=(30-0)/(XTCK(4)-XTCK(3)); b = -a*XTCK(3); XTCKLBL = split(sprintf('%.0f|', round(XTCK*a+b)), '|'); XTCKLBL(end)=[];
YTCK = [13, 296, 337.5, 378, size(img{end},1)-9];a=(10-0)/(YTCK(4)-YTCK(3)); b = -a*YTCK(3); YTCKLBL = split(sprintf('%.0f|', round(YTCK*a+b)), '|'); YTCKLBL(end)=[];

for i = 1:npan
  axes(tightAxes(i));
  image(img{i});
  daspect([1,1,1]);
  TIT = ['$t=',sprintf('%.0f', 1e3*t(i)),'$~ms'];
  title(TIT, 'fontsize', fs);
  xticks(XTCK); xticklabels({});
  yticks(YTCK); yticklabels({});
  if(i==1); ylabel([ftstag, ' orthogonal'], 'fontsize', fs); end
  if(i==5); ylabel([stftag, ' orthogonal'], 'fontsize', fs); end
  if(i==3); ylabel([ftstag, ' slanted'], 'fontsize', fs); end
  if(i==7); ylabel([stftag, ' slanted'], 'fontsize', fs); end
end

hshift = gaph*0.9;
for i=[2,6]; tightAxes(i).Position = tightAxes(i).Position - [hshift,0,0,0]; end
for i=[3,7]; tightAxes(i).Position = tightAxes(i).Position + [hshift,0,0,0]; end
for i=[1,2,5,6]; tightAxes(i).Position = tightAxes(i).Position + [hshift*1.15,0,0,0]; end
for i=[1:4]; tightAxes(i).Position = tightAxes(i).Position - [0,0.01,0,0]; end
for i=[5:8]; tightAxes(i).Position = tightAxes(i).Position - [0,0.01,0,0]; end

set(tightAxes(5), 'xticklabel', XTCKLBL);
set(tightAxes([1,5]), 'yticklabel', YTCKLBL);

ll = add_labels_subplots(gcf,0.9,0,[0,-0.02]);

coeffont = 0.8;
i = 1; x=300; h = annotateAxes(tightAxes(i), [x, size(img{i},1)-200], [x, size(img{i},1)-450], ['incoming P'], fs*coeffont); set(h, 'color', [1,1,1]);
i = 5; x=300; h = annotateAxes(tightAxes(i), [x, size(img{i},1)-470], [x, size(img{i},1)-230], ['incoming P'], fs*coeffont); set(h, 'color', [1,1,1]);

i = 4; x=320; h = annotateAxes(tightAxes(i), [x, size(img{i},1)-290], [x, size(img{i},1)-150], ['reflected P'], fs*coeffont); set(h, 'color', [1,1,1]);
i = 4; x=260; h = annotateAxes(tightAxes(i), [x, size(img{i},1)-400], [x, size(img{i},1)-500], ['transmitted S'], fs*coeffont); set(h, 'color', [1,1,1]);
i = 4; x=400; h = annotateAxes(tightAxes(i), [x, size(img{i},1)-390], [x, size(img{i},1)-450], ['$\quad$','transmitted P'], fs*coeffont); set(h, 'color', [1,1,1]);

i = 7; x=230; h = annotateAxes(tightAxes(i), [x, size(img{i},1)-425], [x, size(img{i},1)-260], ['incoming P'], fs*coeffont); set(h, 'color', [1,1,1]);
i = 7; x=450; h = annotateAxes(tightAxes(i), [x, size(img{i},1)-525], [x, size(img{i},1)-180], ['incoming S (not studied)'], fs*coeffont); set(h, 'color', [1,1,1]);

i = 8; x=260; h = annotateAxes(tightAxes(i), [x, size(img{i},1)-265], [x, size(img{i},1)-150], ['transmitted P'], fs*coeffont); set(h, 'color', [1,1,1]);
i = 8; x=460; h = annotateAxes(tightAxes(i), [x, size(img{i},1)-375], [x, size(img{i},1)-190], ['reflected S'], fs*coeffont); set(h, 'color', [1,1,1]);
i = 8; x=610; h = annotateAxes(tightAxes(i), [x, size(img{i},1)-385], [x, size(img{i},1)-230], ['reflected P'], fs*coeffont); set(h, 'color', [1,1,1]);
i = 8; x=310; h = annotateAxes(tightAxes(i), [x, size(img{i},1)-358], [x, size(img{i},1)-560], ['incoming S (not studied)'], fs*coeffont); set(h, 'color', [1,1,1]);

customSaveFig(fh, figpath, {'fig', 'jpg', 'eps'}, 9999);

% move produced figures to thesis
% disp(['[] Starting to move Figures to thesis folder.']);
% system(['cp ', figpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/images_verif_transmission/']);
