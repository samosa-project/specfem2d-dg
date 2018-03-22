% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
% close all
clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 8); % Default at 6.
set(0, 'defaultTextFontSize', 20);
set(0, 'defaultAxesFontSize', 20); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

FOLDER='/home/l.martire/Documents/MATLAB/snapshot_beautifier_tests/';

if(not(strcmp(FOLDER(end),'/'))); FOLDER=[FOLDER,'/']; end;
list=dir(strcat(FOLDER, '*.jpg'));

i=1;
snap = strcat(FOLDER,list(i).name);

img=imread(snap);
NX=size(img,2);NZ=size(img,1);

figure(i);image(img);set(gca,'DataAspectRatio',[1,1,1]);

crop_left=0;
crop_right=NX-1400;
crop_top=0;
crop_bottom=0;

cropped_img=img(1+crop_top:end-crop_bottom,1+crop_left:end-crop_right,:);
img=cropped_img;
NX=size(img,2);NZ=size(img,1);

close(i);figure(i);image(img);set(gca,'DataAspectRatio',[1,1,1]);

interface_nz_from_top = 580;
zunit="m";
bottom_zval=-300;
top_zval=600;
xunit="m";
left_xval=-300;
right_xval=600;

set(gca,'TickLabelInterpreter', 'latex');
xlabel(strcat("$x$ (",xunit,")"));
ylabel(strcat("$z$ (",zunit,")"));
yticks([1,interface_nz_from_top,NZ]);
yticklabels({strcat("$",num2str(top_zval),"$"),"$0$",strcat("$",num2str(bottom_zval),"$")});
xticks([1,NX]);
xticklabels({strcat("$",num2str(left_xval),"$"),strcat("$",num2str(right_xval),"$")});

xlp = get(get(gca, 'XLabel'), 'Position');set(get(gca, 'XLabel'), 'Position', xlp+[0,-100,0]);
zlp = get(get(gca, 'YLabel'), 'Position');set(get(gca, 'YLabel'), 'Position', zlp+[0,0,-100]);