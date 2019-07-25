% Function coded quickly, maybe to be redone cleaner some day.

function [X,Y,V, imagetype_wavefield_dumps] = readDumps(OFD, IT, verbose)
  if(not(exist('verbose')))
    verbose=0;
  end
%   OFD='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/validation_lns_manufactured/OUTPUT_FILES';
%   IT = 2000;

  % safety.
  OFD = char(OFD);
  if(not(strcmp(OFD(end),filesep)))
    OFD = [OFD,filesep];
  end

  % try identifying NPROC, for safety checks later.
  NPROC = readExampleFiles_extractParam([OFD,'input_parfile'],'NPROC','int');
  
  % check dump type
  imagetype_wavefield_dumps = readExampleFiles_extractParam([OFD,'input_parfile'],'imagetype_wavefield_dumps','int');
  switch(imagetype_wavefield_dumps)
    case {1,2,3}
      nvalues_dumped = 2; % vector dump, 2 values per point
    case 4
      nvalues_dumped = 1; % scalar dump, 1 value per point
    otherwise
      error(['[] imagetype_wavefield_dumps not implemented']);
  end
  if(verbose)
    disp(['[] imagetype_wavefield_dumps = ',num2str(imagetype_wavefield_dumps), ', ',num2str(nvalues_dumped),' to be read per mesh point']);
  end

  wavefield_dump_wildcard = [OFD, 'wavefield*.txt'];

  % dir(wavefield_dump_wildcard)
  wavefield_files = dir(wavefield_dump_wildcard);


  XY = [];
  V = [];
  for i=1:numel(wavefield_files)
    wavefield_file_name = wavefield_files(i).name;
    wavefield_file_path = [wavefield_files(i).folder,filesep,wavefield_file_name];

  %   disp([' >>>>>>>> Reading ''',wavefield_file_name,'''. <<<<<<<<']); % info

    if(regexp(wavefield_file_name,'grid_value_of_nglob'))
      % nglob, do not load it (maybe usefule for checking)
  %     disp('  This is a grid file.');
    elseif(regexp(wavefield_file_name,'grid_for_dumps'))
      % grid, load it;
      inp = importdata(wavefield_file_path);
      if(verbose)
        disp(['grid  file ''',wavefield_file_name,''' needs to be loaded, will be, contains [',num2str(size(inp)),'] values.']);
      end
      XY = [XY; importdata(wavefield_file_path)];
%       size(XY)
    else
  %     disp(['''',wavefield_file_name,''' is a value file, loading it.']);
  
      itname_fills_number_span = regexp(wavefield_file_name,['d',num2str(IT),'_']);
      if(isempty(itname_fills_number_span))
        itname_fills_number_span = 0;
      end
      itname_has_zero_right_before = regexp(wavefield_file_name,['0',num2str(IT),'_']);
      if(isempty(itname_has_zero_right_before))
        itname_has_zero_right_before = 0;
      end
      itname_has_something_else_right_before_and_before_text = regexp(wavefield_file_name,['0+[1-9]0*',num2str(IT),'_']);
      if(isempty(itname_has_something_else_right_before_and_before_text))
        itname_has_something_else_right_before_and_before_text = 0;
      end
      
      if( itname_fills_number_span | (itname_has_zero_right_before & not(itname_has_something_else_right_before_and_before_text)) )
        inp = importdata(wavefield_file_path);
        if(verbose)
          disp(['value file ''',wavefield_file_name,''' has right iteration ID, needs to be loaded, will be, contains [',num2str(size(inp)),'] values.']);
        end
        V = [V;inp];
%         size(V)
      else
  %       disp(['''',wavefield_file_name,''' is a value file with wrong iteration ID, do nothing.']);
      end
    end
  end
  disp(['finished loading']);

  % check
  if(size(V,1)==0)
    error(['[, ERROR] No value has been loaded. Maybe dumps for this iteration (',num2str(IT),') do not exist.']);
  else
    if(not(size(XY,1)==size(V,1)))
      error(['DID NOT LOAD AS MANY POINTS (',num2str(size(XY,1)),') AS VALUES (',num2str(size(V,1)),')']);
    end
  end
  
%   XY = reshape(XY,max(size(XY)),nvalues_dumped);
  if(all(size(V)==[max(size(XY)), nvalues_dumped]) && all(size(XY)==[max(size(XY)), 2]))
    % ok
    if(verbose)
      disp(['[] Shapes are ok, considering dump reading done.']);
    end
  else
    error(['[, ERROR] Shapes (XY=',num2str(size(XY)),', V=',num2str(size(V)),') is wrong w.r.t. nvalues_dumped (=',num2str(nvalues_dumped),')']);
  end
%   V = reshape(V,max(size(V)),1);
  X = reshape(XY(:,1),max(size(XY)),1);
  Y = reshape(XY(:,2),max(size(XY)),1);
end