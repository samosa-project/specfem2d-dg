function [V_out] = readDumps_valuesWrapper(parfile, wvfld_filename, wvfld_path, IT, V_in, bin_nglob, verbose)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PARAMETERS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % RegExps to find CPU.
  rgxp_findProc_step{1} = '_[0-9]+\.bin'; rgxp_findProc_step{2}='[0-9]+';
  % RegExps to find tag, among wholeDumps.
  rgxp_wholeDumps_step{1} = 'wavefield_[a-z]+_'; rgxp_wholeDumps_step{2} = '_[a-z]+_'; rgxp_wholeDumps_step{3} = '[a-z]+';
  % This is determined by the 'write' command in 'write_wavefields_dumps.f90'.
  binaryfile_machinefmt = 'single';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ACTUAL FUNCTION
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Save V_out as is.
  V_out = V_in;
  
  % Check if current file has to be loaded.
  if(verbose)
    disp(['[',mfilename,'] ''',wvfld_filename,''' should be a value file, checking it has right shape, is of right iteration, and then loading it.']);
  end
  [fileNeedsToBeLoaded] = readDumps_valuesCheck(wvfld_filename, IT);
  
  % If needs to be loaded, do it.
  if(fileNeedsToBeLoaded)
    if(verbose)
      disp(['[',mfilename,'] Value file ''',wvfld_filename,''' has right iteration ID, needs to be loaded.']);
    end
    
    % Check whether we will need binary loading.
    is_binary = readExampleFiles_extractParam(parfile, 'use_binary_for_wavefield_dumps', 'bool');
    
    if(is_binary)
      % binary loading
      binfile = fopen(wvfld_path, 'r');
      
      VALUES = fread(binfile, binaryfile_machinefmt); fclose(binfile);
      % identify which proc was that
      curProc = wvfld_path;
      for rgxp=1:numel(rgxp_findProc_step)
        curProc = regexp(curProc, rgxp_findProc_step{rgxp},'match','once');
      end
      curProc = str2double(curProc);
      if(numel(VALUES)==bin_nglob(curProc+1))
        % Loaded the right number of values, all is right.
      else
        error(['[',mfilename,', ERROR] This file is (supposedly) related to PROC n°',num2str(curProc),'. Hence, should have loaded ',num2str(bin_nglob(curProc+1)),' values, but only loaded ',num2str(numel(VALUES)),'.']);
      end
      
    else
      % ascii loading
      VALUES = importdata(wvfld_path);
      
    end

    % Check dump type to add to right field.
    imagetype_wavefield_dumps = readExampleFiles_extractParam(parfile,'imagetype_wavefield_dumps','int');
    switch(imagetype_wavefield_dumps)
      case 1
        current_load = 'rho';
      case 2
        current_load = 'vel';
      case 3
        error(['not implemented']);
      case 4
        current_load = 'pre';
      case 10
        % Special case: whole dumps. Check what was loaded from filename.
        if(verbose)
          disp(['[',mfilename,'] Special case: DG full wavefields dumping.']);
        end
        current_load = wvfld_path;
        for rgxp=1:numel(rgxp_wholeDumps_step)
          current_load = regexp(current_load, rgxp_wholeDumps_step{rgxp},'match','once');
        end
      otherwise
        error(['[] imagetype_wavefield_dumps not implemented']);
    end
    switch(current_load)
      case 'rho'
        nvalues_dumped = 1;
      case 'vel'
        nvalues_dumped = 2;
      case 'pre'
        nvalues_dumped = 1;
      otherwise
        error(['[] current_load not implemented']);
    end
    
    V_out.(current_load) = [V_out.(current_load);VALUES];
    
  else
    if(verbose)
      disp(['[',mfilename,'] Value file ''',wvfld_filename,''' should not be loaded. Skipping.']);
    end
    return
  end
end

