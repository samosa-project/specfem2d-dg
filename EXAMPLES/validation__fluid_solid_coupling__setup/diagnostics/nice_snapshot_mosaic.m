clear all;
close all;
clc;

thisFolder=regexprep(mfilename('fullpath'),mfilename,'');
cd(thisFolder);
addpath(genpath([thisFolder,'/../']));
setupLocal;
setup_find_all_cases; % produces istfortho, istfslant, iftsslant, iftsortho
itodo = [iftsortho, iftsslant, istfortho, istfslant];
folderz = folderz(end-3:end); % consider only the one with highest N (lowest dx)

figpath = [plotFolder,filesep,'snapshot_mosaic'];

% t = [460, 754, 1048, 1342]+800;
% IDs = t/1.4e-3;
% root = ['/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_transmission/'];
% folderz = {[root, 'OUTPUT_FILES_LNS_F2S_plotvz/'], [root, 'OUTPUT_FILES_LNS_S2F_plotvz/']};
IDsFTS = [5, 7000];
IDsSTF = [1000, 4000];

% ncol=2; nlin=4;
% map_ita_to_itodo = [iftsortho, iftsslant, iftsortho, iftsslant, istfortho, istfslant, istfortho, istfslant];
% map_ita_to_IDs = [IDsFTS(1), IDsFTS(1), IDsFTS(2), IDsFTS(2), IDsSTF(1), IDsSTF(1), IDsSTF(2), IDsSTF(2)];

ncol=4; nlin=2;
map_ita_to_itodo = [iftsortho, iftsortho, iftsslant, iftsslant, istfortho, istfortho, istfslant, istfslant];
map_ita_to_IDs = [IDsFTS, IDsFTS, IDsSTF, IDsSTF];

npan=ncol*nlin;

% t = []; c = 1;
for ita = 1:npan
%   for i = 1:numel(IDsFTS)
  id = map_ita_to_IDs(ita);
  OFD=dir([folderz{map_ita_to_itodo(ita)}, filesep, 'OUTPUT_FILES_*']);
  switch(numel(OFD))
    case 0
      error(['found no output files folder']);
    case 1
      % ok
    otherwise
      error(['more than one output files folder']);
  end
  OFD = [OFD.folder,filesep,OFD.name,filesep];

  parfile = [OFD, filesep, 'input_parfile'];
  dt = readExampleFiles_extractParam(parfile, 'DT', 'float');
  t(ita) = id * dt;
  snap = [OFD, filesep, 'image',sprintf('%07d', id),'.jpg'];

  img{ita} = imread(snap);
  NX{ita} = size(img, 2);
  NZ{ita} = size(img, 1);
  
  switch(ita)
    case {1,2}
      cl = 128;
      cr = 476;
      ct = 38;
      cb = 304;
    case {3,4}
      cl = 248;
      cr = 715;
      ct = 98;
      cb = 368;
    case {5,6}
      cl = 154;
      cr = 402;
      ct = 31;
      cb = 177;
    case {7,8}
      cl = 250;
      cr = 660;
      ct = 51;
      cb = 268;
  end

  img{ita} = img{ita}(:,cl:cr,:);
  img{ita} = img{ita}(ct:cb,:,:);
%     c = c+1;
%   end
end

% Plot.
fh = figure('units','normalized','outerposition',[0,0,1,1]);
gaph = 0.035;
tax = tight_subplot(nlin, ncol, [0.05, gaph], [0.05, 0.02], [0.03, 0.005]);
fs=26;

XTCK = [1, size(img{end},2)]; %a=(30-0)/(XTCK(4)-XTCK(3)); b = -a*XTCK(3); XTCKLBL = split(sprintf('%.0f|', round(XTCK*a+b)), '|'); XTCKLBL(end)=[];
YTCK = [1, size(img{end},1)]; %a=(10-0)/(YTCK(4)-YTCK(3)); b = -a*YTCK(3); YTCKLBL = split(sprintf('%.0f|', round(YTCK*a+b)), '|'); YTCKLBL(end)=[];

for i = 1:npan
  axes(tax(i));
  image(img{i});
  daspect([1,1,1]);
  TIT = ['$t=',sprintf('%.0f', 1e3*t(i)),'$~ms'];
  title(TIT, 'fontsize', fs);
  
  switch(i)
    case {1,2}
      XTCK = [166];
      YTCK = [87, 115.5, 145];
    case {3,4}
      XTCK = [342];
      YTCK = [169, 198.5, 258];
    case {5,6}
      XTCK = [122];
      YTCK = [32, 59.5, 87];
    case {7,8}
      XTCK = [260];
      YTCK = [36, 63.5, 93];
  end
  pix2m = 10/diff(YTCK(1:2));
  XTCK = XTCK + [-1,0,1]*30/pix2m;
  switch(i)
    case {1,5,7}
      YTCKLBL = {'', '10', '0', '-10', ''};
    case {3}
      YTCKLBL = {'', '10', '0', '-20', ''};
    otherwise
      YTCKLBL = {};
  end
%   switch(i)
%     case {1,3,5,7}
%       XTCKLBL = {'', '-30', '0', '30', ''};
%     otherwise
%       XTCKLBL = {};
%   end
  XTCKLBL = {'', '-30', '0', '30', ''};
  XTCK = [1, XTCK, size(img{i},2)];
  YTCK = [1, YTCK, size(img{i},1)];

  xticks(XTCK); xticklabels(XTCKLBL);
  yticks(YTCK); yticklabels(YTCKLBL);
  switch(i)
    case 2
      annotation('textbox',[tax(i).Position(1), min(sum(tax(i).Position([2,4]))+.01,1),0,0],'String',[ftstag, ' orthogonal'], 'fontsize', fs, 'FitBoxToText','on','horizontalalignment','center', 'interpreter', 'latex');
    case 4
      annotation('textbox',[tax(i).Position(1), min(sum(tax(i).Position([2,4]))+.01,1),0,0],'String',[ftstag, ' slanted'], 'fontsize', fs, 'FitBoxToText','on','horizontalalignment','center', 'interpreter', 'latex');
    case 6
      annotation('textbox',[tax(i).Position(1), min(sum(tax(i).Position([2,4]))+.01,1),0,0],'String',[stftag, ' orthogonal'], 'fontsize', fs, 'FitBoxToText','on','horizontalalignment','center', 'interpreter', 'latex');
    case 8
      annotation('textbox',[tax(i).Position(1), min(sum(tax(i).Position([2,4]))+.01,1),0,0],'String',[stftag, ' slanted'], 'fontsize', fs, 'FitBoxToText','on','horizontalalignment','center', 'interpreter', 'latex');
  end
  switch(i)
    case {1,3,5,7}
      ylabel('$z$ [m]');
  end
  xlabel('$x$ [m]');
end

hshift = gaph*0.9;
for i=[2,6]; tax(i).Position = tax(i).Position - [hshift,0,0,0]; end
for i=[3,7]; tax(i).Position = tax(i).Position + [hshift,0,0,0]; end
for i=[1,2,5,6]; tax(i).Position = tax(i).Position + [hshift*1.15,0,0,0]; end
for i=[1:4]; tax(i).Position = tax(i).Position - [0,0.01,0,0]; end
for i=[5:8]; tax(i).Position = tax(i).Position - [0,0.01,0,0]; end

% set(tightAxes(5), 'xticklabel', XTCKLBL);
% set(tightAxes([1,5]), 'yticklabel', YTCKLBL);

ll = add_labels_subplots(gcf,0.9,0,[0,-0.02]);

coeffont = 0.8;
% i = 1; x=300; h = annotateAxes(tax(i), [x, size(img{i},1)-200], [x, size(img{i},1)-450], ['incoming P'], fs*coeffont); set(h, 'color', [1,1,1]);
% i = 5; x=300; h = annotateAxes(tax(i), [x, size(img{i},1)-470], [x, size(img{i},1)-230], ['incoming P'], fs*coeffont); set(h, 'color', [1,1,1]);
% 
% i = 4; x=320; h = annotateAxes(tax(i), [x, size(img{i},1)-290], [x, size(img{i},1)-150], ['reflected P'], fs*coeffont); set(h, 'color', [1,1,1]);
% i = 4; x=260; h = annotateAxes(tax(i), [x, size(img{i},1)-400], [x, size(img{i},1)-500], ['transmitted S'], fs*coeffont); set(h, 'color', [1,1,1]);
% i = 4; x=400; h = annotateAxes(tax(i), [x, size(img{i},1)-390], [x, size(img{i},1)-450], ['$\quad$','transmitted P'], fs*coeffont); set(h, 'color', [1,1,1]);
% 
% i = 7; x=230; h = annotateAxes(tax(i), [x, size(img{i},1)-425], [x, size(img{i},1)-260], ['incoming P'], fs*coeffont); set(h, 'color', [1,1,1]);
% i = 7; x=450; h = annotateAxes(tax(i), [x, size(img{i},1)-525], [x, size(img{i},1)-180], ['incoming S (not studied)'], fs*coeffont); set(h, 'color', [1,1,1]);
% 
% i = 8; x=260; h = annotateAxes(tax(i), [x, size(img{i},1)-265], [x, size(img{i},1)-150], ['transmitted P'], fs*coeffont); set(h, 'color', [1,1,1]);
% i = 8; x=460; h = annotateAxes(tax(i), [x, size(img{i},1)-375], [x, size(img{i},1)-190], ['reflected S'], fs*coeffont); set(h, 'color', [1,1,1]);
% i = 8; x=610; h = annotateAxes(tax(i), [x, size(img{i},1)-385], [x, size(img{i},1)-230], ['reflected P'], fs*coeffont); set(h, 'color', [1,1,1]);
% i = 8; x=310; h = annotateAxes(tax(i), [x, size(img{i},1)-358], [x, size(img{i},1)-560], ['incoming S (not studied)'], fs*coeffont); set(h, 'color', [1,1,1]);

customSaveFig(fh, figpath, {'fig', 'jpg', 'eps'}, 9999);

% move produced figures to thesis
% disp(['[] Starting to move Figures to thesis folder.']);
% system(['cp ', figpath, '.* /home/l.martire/Documents/work/THESE/PHD_THESIS/images/chap2/images_verif_transmission/']);
