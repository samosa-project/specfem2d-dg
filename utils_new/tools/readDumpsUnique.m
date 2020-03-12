function [X, Y, V, imagetype_wavefield_dumps] = readDumpsUnique(OFD, IT, verbose)
  [X, Y, V, imagetype_wavefield_dumps] = readDumps(OFD, IT, verbose);
  
  % remove duplicates
  [~,iunique] = unique([X, Y], 'rows');
  X = X(iunique);
  Y = Y(iunique);
  V = V(iunique);
  
  disp(['[',mfilename,'] Made sure dumps only contain unique (x, z) couples.']);
end

