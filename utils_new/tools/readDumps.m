% Function coded quickly, maybe to be redone cleaner some day.

function [X, Y, V] = readDumps(OFD, IT, verbose)
  if(not(exist('verbose', 'var')))
    verbose=0;
  end
%   OFD='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/validation_lns_manufactured/OUTPUT_FILES';
%   IT = 2000;
  
  % Types of variables to query.
  tags = {'rho', 'vel', 'pre'};
  
  % Safeguard.
  OFD = char(OFD);
  if(not(strcmp(OFD(end),filesep)))
    OFD = [OFD,filesep];
  end
  
  parfile=[OFD,'input_parfile'];
  
  % are we indeed dumping dumps?
  output_wavefield_dumps = readExampleFiles_extractParam(parfile, 'output_wavefield_dumps', 'bool');
  if(not(output_wavefield_dumps))
    error(['[',mfilename,', ERROR] Required OUTPUT_FILES (''',OFD,''') is not dumping dumps (as read from parfile ''',parfile,''').']);
  end
  
% %   binary or ascii
  
  % try identifying NPROC, for safety checks later.
  NPROC = readExampleFiles_extractParam(parfile,'NPROC','int');
  
%   % check dump type
%   imagetype_wavefield_dumps = readExampleFiles_extractParam([OFD,'input_parfile'],'imagetype_wavefield_dumps','int');
%   switch(imagetype_wavefield_dumps)
%     case {1,2,3,4}
%       nvalues_dumped = 2; % vector dump, 2 values per point
%       isWholeDumps = 0;
%       current_load = '123';
%     case 4
%       nvalues_dumped = 1; % scalar dump, 1 value per point
%       isWholeDumps = 0;
%       current_load = 'pre';
%     case 10
%       disp(['[] Special case: DG full wavefields dumping. For now, read one value.']);
%       isWholeDumps = 1;
%     otherwise
%       error(['[] imagetype_wavefield_dumps not implemented']);
%   end
% %   if(verbose)
% %     disp(['[] imagetype_wavefield_dumps = ',num2str(imagetype_wavefield_dumps), ', ',num2str(nvalues_dumped),' to be read per mesh point']);
% %   end

  wvflddmp_wildcard_txt = [OFD, 'wavefield*.txt'];
  wvflddmp_wildcard_bin = [OFD, 'wavefield*.bin'];

  % dir(wavefield_dump_wildcard)
  wvflddmp_txtfiles = dir(wvflddmp_wildcard_txt);
  wvflddmp_binfiles = dir(wvflddmp_wildcard_bin);
  
  % Grab numbers of points, and setup remaining files to read.
  is_binary = readExampleFiles_extractParam(parfile, 'use_binary_for_wavefield_dumps', 'bool');
  if(is_binary)
    % If binary, they should be the only .txt files anyhow.
    if(numel(wvflddmp_txtfiles)==NPROC)
      % ok
    else
      error(['[] In this configuration, the number of .txt files should match NPROC from parfile.']);
    end
    bin_nglob = zeros(NPROC, 1);
    for f=1:numel(wvflddmp_txtfiles)
      bin_nglob(f) = importdata([wvflddmp_txtfiles(f).folder,filesep,wvflddmp_txtfiles(f).name]);
    end
    toRead = wvflddmp_binfiles;
  else
    % If ASCII, read all.
    toRead = wvflddmp_txtfiles;
  end

  % Prepare.
  XY = [];
  V = struct();
  for t=1:numel(tags)
    V.(tags{t}) = [];
  end
  
  % Run.
  for i=1:numel(toRead)
    wavefield_file_name = toRead(i).name;
    wavefield_file_path = [toRead(i).folder,filesep,wavefield_file_name];
    if(verbose)
      disp(['[',mfilename,'] Reading ''',wavefield_file_name,'''.']);
    end

    if(regexp(wavefield_file_name,'grid_value_of_nglob'))
      % If binary, files were read beforehand and should not be reached here.
      % If ASCII, we can still reach them. TODO: read nglob above, alongside binary reading.
      
    elseif(regexp(wavefield_file_name,'grid_for_dumps'))
      % This is a grid file, load it.
      if(verbose)
        disp(['[',mfilename,'] Grid file. ''',wavefield_file_name,'''. Needs to be loaded.']);
%         pause
      end
      XY = [XY; readDumps_grid(is_binary, wavefield_file_path)];
      
    else
      % This should be a file containing values. Update V with the values in this file.
      V = readDumps_valuesWrapper(parfile, wavefield_file_name, wavefield_file_path, IT, V, bin_nglob, verbose);
      
    end
  end % end on files to read
  disp(['[',mfilename,'] Finished loading dumps, and assembling them.']);

  % check
  for t=1:numel(tags)
    if(size(V.(tags{t}),1)==0)
      disp(['[',mfilename,', INFO] No value has been loaded in V.',tags{t},'. Maybe ',tags{t},' dumps for this iteration (',num2str(IT),') do not exist.']);
    else
      if(not(size(XY,1)==size(V.(tags{t}),1)))
        error(['DID NOT LOAD AS MANY POINTS (',num2str(size(XY,1)),') AS VALUES (',num2str(size(V,1)),')']);
      end
    end
  end
  
% %   XY = reshape(XY,max(size(XY)),nvalues_dumped);
%   if(all(size(V)==[max(size(XY)), nvalues_dumped]) && all(size(XY)==[max(size(XY)), 2]))
%     % ok
%     if(verbose)
%       disp(['[] Shapes are ok, considering dump reading done.']);
%     end
%   else
%     error(['[, ERROR] Shapes (XY=',num2str(size(XY)),', V=',num2str(size(V)),') is wrong w.r.t. nvalues_dumped (=',num2str(nvalues_dumped),')']);
%   end
  
%   V = reshape(V,max(size(V)),1);
  X = reshape(XY(:,1),max(size(XY)),1);
  Y = reshape(XY(:,2),max(size(XY)),1);
end