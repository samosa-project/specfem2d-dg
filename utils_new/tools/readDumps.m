% Function coded quickly, maybe to be redone cleaner some day.

function [X, Y, V, imagetype_wavefield_dumps] = readDumps(OFD, IT, verbose)
  if(not(exist('verbose', 'var')))
    verbose=0;
  end
%   OFD='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/validation_lns_manufactured/OUTPUT_FILES';
%   IT = 2000;

  % safety.
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
  
  % binary or ascii
  use_binary_for_wavefield_dumps = readExampleFiles_extractParam(parfile, 'use_binary_for_wavefield_dumps', 'bool');
  
  % try identifying NPROC, for safety checks later.
  NPROC = readExampleFiles_extractParam(parfile,'NPROC','int');
  
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

  wvflddmp_wildcard_txt = [OFD, 'wavefield*.txt'];
  wvflddmp_wildcard_bin = [OFD, 'wavefield*.bin'];

  % dir(wavefield_dump_wildcard)
  wvflddmp_txtfiles = dir(wvflddmp_wildcard_txt);
  wvflddmp_binfiles = dir(wvflddmp_wildcard_bin);
  
  if(use_binary_for_wavefield_dumps)
    % If binary, grab numbers of points. They should be the only .txt files anyhow.
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
    % If ASCII.
    toRead = wvflddmp_txtfiles;
  end

  XY = [];
  V = [];
  for i=1:numel(toRead)
    wavefield_file_name = toRead(i).name;
    wavefield_file_path = [toRead(i).folder,filesep,wavefield_file_name];

    disp(['[',mfilename,'] Reading ''',wavefield_file_name,'''.']);
%     pause

    if(regexp(wavefield_file_name,'grid_value_of_nglob'))
      % If binary, no such file exist. If ASCII, we do not care (may be useful for checking). Skip.
  %     disp('  This is a grid file.');
    elseif(regexp(wavefield_file_name,'grid_for_dumps'))
      % This is a grid file, load it;
      if(verbose)
        disp(['grid  file ''',wavefield_file_name,''' needs to be loaded.']);
%         pause
      end
      if(use_binary_for_wavefield_dumps)
        % binary loading
        binfile = fopen(wavefield_file_path, 'r');
        machinefmt = 'single'; % This is determined by the 'write' command in 'write_wavefields_dumps.f90'.
        inp = fread(binfile, machinefmt); fclose(binfile);
        inp = reshape(inp, [2, numel(inp)/2])'; % Reshape HAS TO be like this, because binary file lists points as 'x,y,x,y,x,y' and hence first loading of inp stores them as [x;y;x;y;x;y]. N.B.: 'numel/2' would crash if wrong number of singles loaded.
      else
        % ascii loading
        inp = importdata(wavefield_file_path);
      end
      XY = [XY; inp];
    else
      
      % This should be a file containing values.
      disp(['[',mfilename,'] ''',wavefield_file_name,''' should be a value file, checking it has right shape, is of right iteration, and then loading it.']);
%       pause
  
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
        if(verbose)
          disp(['[',mfilename,'] Value file ''',wavefield_file_name,''' has right iteration ID, needs to be loaded.']);
%           pause
        end
        if(use_binary_for_wavefield_dumps)
          % binary loading
          binfile = fopen(wavefield_file_path, 'r');
          machinefmt = 'single'; % This is determined by the 'write' command in 'write_wavefields_dumps.f90'.
          inp = fread(binfile, machinefmt); fclose(binfile);
          % identify which proc was that
          curProc = str2double(regexp(regexp(wavefield_file_path,'_[0-9]+\.bin','match','once'),'[0-9]+','match','once'));
          if(numel(inp)==bin_nglob(curProc+1))
            % Loaded the right number of values, all is right.
          else
            error(['[',mfilename,', ERROR] This file is (supposedly) related to PROC n°',num2str(curProc),'. Hence, should have loaded ',num2str(bin_nglob(curProc+1)),' values, but only loaded ',num2str(numel(inp)),'.']);
          end
        else
          % ascii loading
          inp = importdata(wavefield_file_path);
        end
%         if(verbose)
%           disp(['[',mfilename,'] Value file ''',wavefield_file_name,''' has right iteration ID, needs to be loaded, will be, contains [',num2str(size(inp)),'] values.']);
%         end
        V = [V;inp];
%         size(V)
      else
        if(verbose)
          disp(['[',mfilename,'] Value file ''',wavefield_file_name,''' should not be loaded. Skipping.']);
        end
      end
    end
  end
  disp(['[',mfilename,'] Finished loading dumps, and assembling them.']);

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