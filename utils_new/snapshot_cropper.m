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

FOLDER = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/SH_final/snapshots_edits/kushida';
% FOLDER = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200/snapshots_edits';
% FOLDER='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/OKQ/snapshots_edits/cropping';
% FOLDER='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/SH_final/snapshots_edits/cropping';
% FOLDER='/home/l.martire/Documents/MATLAB/snapshot_beautifier_tests/';

if(not(strcmp(FOLDER(end),'/'))); FOLDER=[FOLDER,'/']; end;
list=dir(strcat(FOLDER, '*.jpg'));
rem=[];
for i=1:size(list,1)
  if(not(isempty(regexp(list(i).name,suffix,'ONCE'))))
    rem=[rem,i]; % Remove already cropped snapshots.
  end
end
list(rem)=[];
disp(strcat("Treating folder '",FOLDER,"', containing ",num2str(size(list,1))," snapshot files."));

fs=input('Font size (default 24)? > '); set(0, 'defaultTextFontSize', fs); set(0, 'defaultAxesFontSize', fs);

newsnap={};
for s=1:size(list,1)
  snap = strcat(FOLDER,list(s).name);

  img=imread(snap); NX=size(img,2);NZ=size(img,1);
  figure(s);
  
  if(s==1)
    image(img); set(gca,'DataAspectRatio',[1,1,1]);
    disp("Run size.");
    xunit=input('  X unit? > ', 's'); zunit=input('  Z unit? > ', 's');
    XZunit=input('  (X,Z) min/max values, x_source, and z_interface, in unit? Format [min(x), max(x), min(z), max(z), x_source, z_interface].\n    > ');
    xmin=XZunit(1); xmax=XZunit(2); zmin=XZunit(3); zmax=XZunit(4); xs=XZunit(5); zint=XZunit(6);
    xspan=xmax-xmin; zspan=zmax-zmin;
    nx_s=(xs-xmin)*NX/xspan; nz_int=(zmax-zint)*NZ/zspan;
    disp("Cropping.");
    satisfied=0;
    while(satisfied~=1)
      XZcrop=input('  (X,Z) crop, in pixels? Format [left, right, bottom, top]. > ');
      crop_left=XZcrop(1); crop_right=XZcrop(2); crop_bottom=XZcrop(3); crop_top=XZcrop(4);
      cropped_img=img(1+crop_top:end-crop_bottom,1+crop_left:end-crop_right,:);
      close(s);f=figure(s);image(cropped_img);set(gca,'DataAspectRatio',[1,1,1]);
      satisfied=input('  Satisfied (0 for no, 1 for yes)? > ');
    end
    
    fconfig=fopen(strcat(FOLDER,"/config"),"w");
    fprintf(fconfig,strcat(xunit," ",zunit));
    fprintf(fconfig,"\n");
    fprintf(fconfig,num2str(XZunit));
    fprintf(fconfig,"\n");
    fprintf(fconfig,num2str(XZcrop));
    fclose("all");
  end

  cropped_img=img(1+crop_top:end-crop_bottom,1+crop_left:end-crop_right,:);
  img=cropped_img;
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
  yticks([1,interface_nz_from_top,nNZ]);
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
  tightfig;

  splsnap=split(snap,"/"); splsnapend=split(splsnap(end),"."); splsnapend{1}=[splsnapend{1},char(suffix)];splsnap{end}=char(join(splsnapend,"."));newsnap{s}=char(join(splsnap,"/"));
  saveas(gcf,newsnap{s});
end

if(0==1) % One-liner to be ran to correct axis label positionning on all cropped snapshots.
  for s=1:size(list,1); figure(s); xlp = get(get(gca, 'XLabel'), 'Position'); set(get(gca, 'XLabel'), 'Position', xlp+[-90,20,0]); saveas(gcf,newsnap{s}); end
end