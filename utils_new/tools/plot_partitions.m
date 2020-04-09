% Author:        Léo Martire.
% Description:   Use this script to plot partitions of a SPECFEM run.
% Notes:         TODO.
%
% Usage:
%   Just run it.

function plot_partitions(PseudoGMSH0orDATABASE1, path_to_mesh, boundary_only, ids_choice)
  
  % Parameters.
  boundarythreshold = 1;
  ptsize    = 5;
  dofill    = 0;
  fillalpha = 0.25;
  figsize   = 1000;
  
  % Input parsing.
  if(not(exist('PseudoGMSH0orDATABASE1', 'var')))
    PseudoGMSH0orDATABASE1=-1;
    while(not(numel(PseudoGMSH0orDATABASE1)==1 & ismember(PseudoGMSH0orDATABASE1,[0,1])))
      PseudoGMSH0orDATABASE1 = input(['[',mfilename,'] Pseudo-GMSH plot (0) or SPECFEM Database plot (1)? > ']);
    end
  end
  if(not(exist('path_to_mesh', 'var')))
    path_to_mesh_provided=0;
  else
    path_to_mesh_provided=1;
  end
  if(not(exist('boundary_only','var')))
    boundary_only=-1;
    while(not(ismember(boundary_only,[0,1])))
      boundary_only=input(['[',mfilename,'] For each partition, plot all points (0, might be very heavy) or only boundary (1)? > ']);
    end
  end
  if(not(exist('ids_choice', 'var')))
    ids_choice_provided = 0;
  else
    ids_choice_provided = 1;
  end
  
  % Error checking.
  if(not(ismember(boundary_only,[0,1])))
    error(['[',mfilename,', ERROR] Boundary only switch must be either 0 (for all points), or 1 (for boundary only).']);
  end
  
  % Load paths to mesh.
  switch(PseudoGMSH0orDATABASE1)
    case 0
      if(not(path_to_mesh_provided))
        path_to_mesh = input(['[',mfilename,'] Path to Nodes_extMesh? > '],'s');
        if(not(exist(path_to_mesh)==2))
          error(['[',mfilename,', ERROR] ''',path_to_mesh,''' is not a file.']);
        end
      end
    case 1
      if(not(path_to_mesh_provided))
        path_to_mesh = input(['[',mfilename,'] Folder containing Database files (usually your run''s OUTPUT_FILES folder)? > '],'s');
      end
      if(not(path_to_mesh(end)=='/')) path_to_mesh=[path_to_mesh,'/']; end;
    otherwise
      error(['[',mfilename,', ERROR] First parameter must be either 0, or 1.']);
  end
  % Check existence.
  switch(PseudoGMSH0orDATABASE1)
    case 0
      if(not(exist(path_to_mesh,'file') && not(exist(path_to_mesh,'dir'))))
        % if not at least a file (must check it is not a folder, because a folder will check out also as file), error
        error(['[',mfilename,', ERROR] If first parameter is ',num2str(PseudoGMSH0orDATABASE1),', path parameter must be a file.']);
      end
    case 1
      if(not(exist(path_to_mesh,'dir')))
        % if not at least a folder, error
        error(['[',mfilename,', ERROR] If first parameter is ',num2str(PseudoGMSH0orDATABASE1),', path parameter must be a folder.']);
      end
    otherwise
      error(['[',mfilename,', ERROR] First parameter must be either 0, or 1.']);
  end
  
  % Actual treatment.
  switch(PseudoGMSH0orDATABASE1)
    case 0
    %   error(['[',mfilename,', ERROR] Not implemented yet.']);
      [P] = readExampleFiles_extmshLoadNodes(path_to_mesh);

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
    case 1
      database_list=dir([path_to_mesh,'Database*']);
      if(isempty(database_list))
        error(['[',mfilename,', ERROR] No ''Database'' file found.']);
      end
      numdatabases=numel(database_list);
      CPUs = [];
      for i = 1:numdatabases
        found_cpu = regexp(database_list(i).name,'[0-9]+','match'); 
        CPUs = [CPUs, str2num(found_cpu{1})]; % save cpu id
        X{i} = scan_database_file([path_to_mesh,database_list(i).name]);
        if(boundary_only)
          j = boundary(X{i}, boundarythreshold);
          X{i} = X{i}(j,:);
        end
      end
      arr_ids_cpus_original=1:numdatabases;

      if(not(ids_choice_provided))
        ids_choice=-10;
        while(not((numel(ids_choice)==1 && ids_choice==-1) || all(ismember(ids_choice,arr_ids_cpus_original))))
          disp(['[',mfilename,'] Database to plot ? Must be either -1 to plot all, or an array of IDs in']);
          disp(['[',mfilename,']   CPU = ',sprintf('%5d',CPUs)]);
          disp(['[',mfilename,']   ID  = ',sprintf('%5d',arr_ids_cpus_original), ' <- choose among those']);
          ids_choice = input(['[',mfilename,'] > ']);
        end
        if(ids_choice==-1)
          arr_ids_cpus = arr_ids_cpus_original;
        else
          arr_ids_cpus = ids_choice;
        end
      else
        % user provided cpu ids, find them
        arr_ids_cpus = find(CPUs==ids_choice);
      end

      f=figure();
      f.set('pos',[10 10 figsize, figsize]);
      % colorvect=jet(numdatabases);
      linescm=lines; colorvect=linescm(1:numdatabases,:); clear('linescm');
      for i=arr_ids_cpus
        cpuname=['P',num2str(CPUs(i))];
        if(boundary_only)
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
        for j=arr_ids_cpus(arr_ids_cpus~=i)
          inboth=ismember(X{i},X{j},'rows');
          h=scatter(X{i}(inboth,1),X{i}(inboth,2),ptsize*2*(j+1),repmat(colorvect(j,:),size(X{i}(inboth,1),1),1),'linewidth',1);
          set(h,'handlevisibility','off');
        end
      end
      legend('location', 'eastoutside');
%       set(gca, 'tickdir','both');
%       set(gca, 'TickLabelInterpreter','latex');
      xlabel('$x$ (m)'); ylabel('$z$ (m)');
      %title({['Partitions'],FOLDER});
      %axis square
      set(gca,'DataAspectRatio',[1,1,1]);
      
      if(numel(arr_ids_cpus)==1)
        plural = '';
      else
        plural = 's';
      end
      filename = sprintf('%05d_',sort(CPUs(arr_ids_cpus)));
      filename(end)=[];
      filename = ['partition',plural,'_',filename];
      
      customSaveFig(gcf,[path_to_mesh,filename],{'fig','jpg'});
%       exts=["eps","png"];
%       formats=["epsc","png"];
%       for i=1:numel(exts)
%         saveas(gcf,[path_to_mesh,'partitions.',char(exts(i))],char(formats(i)));
%       end
%       savefig(gcf,[path_to_mesh,'partitions.fig']);
      
    otherwise
      error(['[',mfilename,', ERROR] First parameter must be either 0, or 1.']);
  end
end