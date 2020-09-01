function [XY] = readDumps_grid(use_binary_for_wavefield_dumps, wavefield_file_path)
  if(use_binary_for_wavefield_dumps)
    % binary loading
    binfile = fopen(wavefield_file_path, 'r');
    
    machinefmt = 'single'; % This is determined by the 'write' command in 'write_wavefields_dumps.f90'.
    
    XY = fread(binfile, machinefmt);
    
    fclose(binfile);
    
    XY = reshape(XY, [2, numel(XY)/2])'; % Reshape HAS TO be like this, because binary file lists points as 'x,y,x,y,x,y' and hence first loading of inp stores them as [x;y;x;y;x;y]. N.B.: 'numel/2' would crash if wrong number of singles loaded.
    
  else
    % ascii loading
    XY = importdata(wavefield_file_path);
    
  end
end

