% Author:        Léo Martire.
% Description:   Use this script to plot partitions of a SPECFEM run.
% Notes:         TODO.
%
% Usage:
%   Just run it.

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

PseudoGMSH0orDATABASE1=-1;
while(not(numel(PseudoGMSH0orDATABASE1)==1 & ismember(PseudoGMSH0orDATABASE1,[0,1])))
  PseudoGMSH0orDATABASE1 = input(['[',mfilename,'] Pseudo-GMSH plot (0) or SPECFEM Database plot (1)? > ']);
end

if(PseudoGMSH0orDATABASE1==0)
%   error(['[',mfilename,', ERROR] Not implemented yet.']);
%   Nodes_extMesh = input(['[',mfilename,'] Path to Nodes_extMesh? > '],'s');
%   Nodes_extMesh_f = fopen(Nodes_extMesh,'r');
%   lines = textscan(Nodes_extMesh_f, '%f %f');
% %   X = lines{1}(2:end); Z = lines{2}(2:end); P = [X,Z];
%   P = [lines{1}(2:end), lines{2}(2:end)];
%   if( not(lines{1}(1)==size(P, 1)))
%     error(['[',mfilename,', ERROR] Issue while reading Nodes_extMesh: number of points read (number of lines) does not correspond to number of points specified (first line).']);
%   end
  [P] = extMesh_loadNodes(Nodes_extMesh);
  
  doWindowX = -1;
  while (not((numel(doWindowX)==1 && doWindowX==0) || (numel(doWindowX)==2)))
    doWindowX = input(['[', mfilename, '] Window X ([x1 x2] vector for yes, 0 for no)? > ']);
  end
  doWindowZ = -1;
  while (not((numel(doWindowZ)==1 && doWindowZ==0) || (numel(doWindowZ)==2)))
    doWindowZ = input(['[', mfilename, '] Window Z ([z1 z2] vector for yes, 0 for no)? > ']);
  end
  
  if(any(doWindowX))
    sel = (P(:,1)>=min(doWindowX) & P(:,1)<=max(doWindowX));
    P = P(sel,:);
  end
  if(any(doWindowZ))
    sel = (P(:,2)>=min(doWindowZ) & P(:,2)<=max(doWindowZ));
    P = P(sel,:);
  end
  
  f = figure();
  plot(P(:,1), P(:,2), '.');
  daspect([1 1 1]);
  margin=0.05*min((max(P(:,1))-min(P(:,1))), (max(P(:,2))-min(P(:,2))));
  xlim([min(P(:,1)),max(P(:,1))]+margin*[-1,1]);
  ylim([min(P(:,2)),max(P(:,2))]+margin*[-1,1]);
  prettyAxes(f);
end

if(PseudoGMSH0orDATABASE1==1)
  FOLDER = input(['[',mfilename,'] Folder containing Database files (usually your run''s OUTPUT_FILES folder)? > '],'s');
  if(not(FOLDER(end)=='/')) FOLDER=[FOLDER,'/']; end;

  bound=-1;
  while(not(ismember(bound,[0,1])))
    bound=input(['[',mfilename,'] For each partition, plot all points (0, might be very heavy) or only boundary (1)? > ']);
  end

  database_list=dir([FOLDER,'Database*']);
  if(isempty(database_list))
    error(['[',mfilename,', ERROR] No ''Database'' file found.']);
  end
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
end