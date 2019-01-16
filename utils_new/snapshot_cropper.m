% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
close all;
clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 8); % Default at 6.
set(0, 'defaultTextFontSize', 24);
set(0, 'defaultAxesFontSize', 24); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

addpath('/usr/local/matlab/r2017b/toolbox/tightfig');

suffix="_crop";

% FOLDER = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/SH_final/snapshots_edits/kushida';
% FOLDER = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200/snapshots_edits';
% FOLDER='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/OKQ/snapshots_edits/cropping';
% FOLDER='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/SH_final/snapshots_edits/cropping';
% FOLDER='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/OKQ/snapshots_edits/redone';
% FOLDER='/home/l.martire/Downloads/gifs/exp_ballons_helium/shhard';
% FOLDER='/home/l.martire/Downloads/gifs/exp_ballons_helium/shsoft';
% FOLDER='/home/l.martire/Documents/MATLAB/snapshot_beautifier_tests/';
% FOLDER='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/tir_de_mine/OUTPUT_FILES_74752/illustration';
% FOLDER='/home/l.martire/Documents/SPECFEM/Ongoing_Work/18_microbaroms/microbaroms_patch/OUTPUT_FILES_668482_disp7_isp6_full/cropppp';
% FOLDER='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/demo_lns/OUTPUT_FILES_826213/croppp';
% FOLDER='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/demo_lns/OUTPUT_FILES_lnsf2s_local/croppp';
FOLDER='/home/l.martire/Documents/Ongoing_Work/1811_glanes/190116_presentation_JPL/simulation';

if(not(strcmp(FOLDER(end),'/'))); FOLDER=[FOLDER,'/']; end;
list=dir(strcat(FOLDER, '*.jpg'));
rem=[];
for i=1:size(list,1)
  if(not(isempty(regexp(list(i).name,suffix,'ONCE'))))
    rem=[rem,i]; % Remove already cropped snapshots.
  end
end
list(rem)=[];
disp(strcat("[",mfilename,"] Treating folder '",FOLDER,"', containing ",num2str(size(list,1))," snapshot files."));

save_raw=-1;
while(not(ismember(save_raw,[0,1])))
  save_raw=input(['[',mfilename,'] > Which type of cropped snapshots (0: formatted, 1: raw)? > ']);
end

ext="KEK";
if(save_raw)
  while(not(ismember(ext,["jpg","png"])))
    ext=input(['[',mfilename,'] > Extension (jpg, png)? > '],'s');
  end
else
  while(not(ismember(ext,["fig","jpg","png","eps"])))
    ext=input(['[',mfilename,'] > Extension (fig, jpg, png, eps)? > '],'s');
  end
  if(strcmp(ext,"eps")); ext="epsc"; end
  ext=char(ext);
end

if(not(save_raw))
  fs=input(['[',mfilename,'] > Font size (default 24)? > ']);
  set(0, 'defaultTextFontSize', fs);
  set(0, 'defaultAxesFontSize', fs);
end

newsnap={};
for s=1:size(list,1)
  snap = strcat(FOLDER,list(s).name);
  
  %   splsnap=split(snap,"/"); splsnapend=split(splsnap(end),"."); splsnapend{1}=[splsnapend{1},char(suffix)]; splsnapend{2}=''; splsnap{end}=char(join(splsnapend,".")); newsnap{s}=char(join(splsnap,"/"));
  splsnap=split(snap,"/"); splsnapend=split(splsnap(end),"."); splsnap{end}=[splsnapend{1},char(suffix)]; newsnap{s}=char(join(splsnap,"/"));
  
  img=imread(snap);
  NX=size(img,2);
  NZ=size(img,1);
  
  if(not(save_raw))
    figure(s);
  end
  
  if(s==1)
    image(img); set(gca,'DataAspectRatio',[1,1,1]);
    if(not(save_raw))
      disp(['[',mfilename,'] Run size.']);
      xunit=input(['[',mfilename,'] > X unit? > '], 's');
      zunit=input(['[',mfilename,'] > Z unit? > '], 's');
      XZunit=input(['[',mfilename,'] > (X,Z) min/max values, x_source, and z_interface, in unit? Format [min(x), max(x), min(z), max(z), x_source, z_interface].\n    > ']);
      xmin=XZunit(1); xmax=XZunit(2); zmin=XZunit(3); zmax=XZunit(4); xs=XZunit(5); zint=XZunit(6);
      xspan=xmax-xmin; zspan=zmax-zmin;
      nx_s=(xs-xmin)*NX/xspan; nz_int=(zmax-zint)*NZ/zspan;
    end
    disp(['[',mfilename,'] Cropping.']);
    satisfied=0;
    while(satisfied~=1)
      XZcrop=input(['[',mfilename,'] > (X,Z) crop, in pixels? Format [left, right, bottom, top]. > ']);
      crop_left=XZcrop(1); crop_right=XZcrop(2); crop_bottom=XZcrop(3); crop_top=XZcrop(4);
      cropped_img=img(1+crop_top:end-crop_bottom,1+crop_left:end-crop_right,:);
      close(s);f=figure(s);image(cropped_img);set(gca,'DataAspectRatio',[1,1,1]);
      satisfied=input(['[',mfilename,'] > Satisfied (0 for no, 1 for yes)? > ']);
    end
    
    fconfig=fopen(strcat(FOLDER,"/config"),"w");
    if(not(save_raw))
      fprintf(fconfig,num2str(fs));
      fprintf(fconfig,"\n");
      fprintf(fconfig,strcat(xunit," ",zunit));
      fprintf(fconfig,"\n");
      fprintf(fconfig,num2str(XZunit));
      fprintf(fconfig,"\n");
    end
    fprintf(fconfig,num2str(XZcrop));
    fclose("all");
  end

  cropped_img=img(1+crop_top:end-crop_bottom,1+crop_left:end-crop_right,:);
  img=cropped_img;
  if(save_raw)
    imwrite(img, strcat(newsnap{s}, '.', ext));
  else
    nNX=size(img,2);nNZ=size(img,1);
    close(s); f=figure(s); image(img); set(gca,'DataAspectRatio',[1,1,1]);

    if(s==1)
      left_xval=xmin+crop_left*xspan/NX; right_xval=xmax-crop_right*xspan/NX;
      bottom_zval=zmin+crop_bottom*zspan/NZ; top_zval=zmax-crop_top*zspan/NZ;
      source_nx_from_left=nx_s-crop_left; interface_nz_from_top=nz_int-crop_top;
    end

    set(gca,'TickLabelInterpreter', 'latex');
    xlabel(strcat("$x$ (",xunit,")"));
    ylabel(strcat("$z$ (",zunit,")"));
    if(interface_nz_from_top<nNZ)
      yticks([1,interface_nz_from_top,nNZ]);
    else
      yticks([1,nNZ]);
    end
    yticklabels({strcat(" $",sprintf("%.1f",top_zval-zint),"$")," $0$",strcat(" $",sprintf("%.1f",bottom_zval-zint),"$")});
    xticks([1,source_nx_from_left,nNX]);
    xticklabels({strcat(" $",sprintf("%.1f",left_xval-xs),"$")," $0$",strcat("$",sprintf("%.1f",right_xval-xs),"$ ")});

    xlp = get(get(gca, 'XLabel'), 'Position'); zlp = get(get(gca, 'YLabel'), 'Position');
    set(get(gca, 'XLabel'), 'Position', [ xlp(1),    2.02*zlp(2), 0]);
    set(get(gca, 'YLabel'), 'Position', [-xlp(1)/18, zlp(2),      0]);
    set(gca,'TickDir','out');

    f.set('pos',[10 10 1.1*nNX 1.1*nNZ]);
    figure(s);
    gcf;
%     tightfig;

    saveas(gcf,newsnap{s},ext);
  end
end

if(0==1) % One-liners to be ran to correct axis label positionning on all cropped snapshots.
  for s=1:size(list,1); figure(s); xlp = get(get(gca, 'XLabel'), 'Position'); set(get(gca, 'XLabel'), 'Position', xlp+[0,10,0]); saveas(gcf,newsnap{s},ext); end % Shift XLabel vertically.
  for s=1:size(list,1); figure(s); xlp = get(get(gca, 'XLabel'), 'Position'); set(get(gca, 'XLabel'), 'Position', xlp+[10,0,0]); saveas(gcf,newsnap{s},ext); end % Shift XLabel horizontally.
  for s=1:size(list,1); figure(s); zlp = get(get(gca, 'YLabel'), 'Position'); set(get(gca, 'YLabel'), 'Position', zlp+[0,-40,0]); saveas(gcf,newsnap{s},ext); end % Shift YLabel vertically.
  for s=1:size(list,1); figure(s); zlp = get(get(gca, 'YLabel'), 'Position'); set(get(gca, 'YLabel'), 'Position', zlp+[-40,0,0]); saveas(gcf,newsnap{s},ext); end % Shift YLabel horizontally.
end
