function [X,Y,V] = readDumps(OFD, IT)
%   OFD='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/validation_lns_manufactured/OUTPUT_FILES';
%   IT = 2000;

  % safety.
  OFD = char(OFD);
  if(not(strcmp(OFD(end),filesep)))
    OFD = [OFD,filesep];
  end

  % try identifying NPROC, for safety checks later.
  NPROC = readExampleFiles_extractParam([OFD,'input_parfile'],'NPROC','int');

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
      disp(['grid  file ''',wavefield_file_name,''' needs to be loaded, will be, needs to be loaded, will be, contains [',num2str(size(inp)),'] values.']);
      XY = [XY; importdata(wavefield_file_path)];
%       size(XY)
    else
  %     disp(['''',wavefield_file_name,''' is a value file, loading it.']);
      if( regexp(wavefield_file_name,['0',num2str(IT),'_']) | regexp(wavefield_file_name,['d',num2str(IT),'_']))
        inp = importdata(wavefield_file_path);
        disp(['value file ''',wavefield_file_name,''' has right iteration ID, needs to be loaded, will be, contains [',num2str(size(inp)),'] values.']);
        V = [V;inp];
%         size(V)
      else
  %       disp(['''',wavefield_file_name,''' is a value file with wrong iteration ID, do nothing.']);
      end
    end
  end
  disp(['finished loading']);

  XY = reshape(XY,max(size(XY)),2);
  V = reshape(V,max(size(V)),1);
  % check
  if(not(size(XY,1)==size(V,1)))
    error(['DID NOT LOAD AS MANY POINTS (',num2str(size(XY,1)),') AS VALUES (',num2str(size(V,1)),')']);
  end
  X = reshape(XY(:,1),max(size(XY)),1);
  Y = reshape(XY(:,2),max(size(XY)),1);
end