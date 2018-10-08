% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Use this script to plot partitions of a SPECFEM run.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
% close all;
% clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

boundarythreshold=1;
ptsize=10;
dofill=0;
fillalpha=0.25;
figsize=1000;

FOLDER=input('Folder containing Database files? > ','s');
if(not(FOLDER(end)=='/')) FOLDER=[FOLDER,'/']; end;

bound=-1;
while(not(ismember(bound,[0,1])))
  bound=input('For each partition, plot all points (0, might be very heavy) or only boundary (1)? > ');
end

database_list=dir([FOLDER,'Database*']);
numdatabases=numel(database_list);
for i=1:numdatabases
  X{i}=scan_database_file([FOLDER,database_list(i).name]);
  if(bound)
    j=boundary(X{i},boundarythreshold);
    X{i}=X{i}(j,:);
  end
end
arrcpus=1:numdatabases;

f=figure();
f.set('pos',[10 10 figsize, figsize]);
% colorvect=jet(numdatabases);
linescm=lines; colorvect=linescm(1:numdatabases,:); clear('linescm');
for i=arrcpus
  cpuname=['P',num2str(i-1)];
  if(bound)
    j=1:size(X{i},1);
    h=plot(X{i}(:,1),X{i}(:,2),'.-','markersize',ptsize,'color',colorvect(i,:)); hold on;
  else
    j=boundary(X{i},boundarythreshold);
    h=scatter(X{i}(:,1),X{i}(:,2),ptsize,repmat(colorvect(i,:),size(X{i},1),1),'filled'); hold on;
  end
  if(dofill)
    hfill=fill(X{i}(j,1),X{i}(j,2),colorvect(i,:)); hold on;
    set(hfill,'facealpha',fillalpha,'displayname',cpuname);
  else
    scatter(X{i}(j(1),1),X{i}(j(1),2),ptsize,colorvect(i,:),'filled','displayname',cpuname); hold on;
  end
  set(h,'handlevisibility','off');
  for j=arrcpus(arrcpus~=i)
    inboth=ismember(X{i},X{j},'rows');
    h=scatter(X{i}(inboth,1),X{i}(inboth,2),ptsize*2*(j+1),repmat(colorvect(j,:),size(X{i}(inboth,1),1),1),'linewidth',1);
    set(h,'handlevisibility','off');
  end
end
legend;
set(gca, 'tickdir','both');
set(gca, 'TickLabelInterpreter','latex');
xlabel('$x$ (m)'); ylabel('$z$ (m)');
%title({['Partitions'],FOLDER});
%axis square
set(gca,'DataAspectRatio',[1,1,1]);
exts=["eps","png"];
formats=["epsc","png"];
for i=1:numel(exts)
  saveas(gcf,[FOLDER,'partitions.',char(exts(i))],char(formats(i)));
end
savefig(gcf,[FOLDER,'partitions.fig']);