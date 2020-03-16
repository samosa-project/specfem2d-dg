function [X, Y, V] = readDumpsUnique(OFD, IT, verbose)
  [X, Y, V] = readDumps(OFD, IT, verbose);
  
  tags = {'rho', 'vel', 'pre'};
  
  % Find duplicate points.
  [~, iunique] = unique([X, Y], 'rows');
  
  % Remove duplicate points.
  X = X(iunique);
  Y = Y(iunique);
  
  % Remove duplicates, in each non-empty field.
  for t=1:numel(tags)
    if(not(isempty(V.(tags{t}))))
      V.(tags{t}) = V.(tags{t})(iunique);
    end
  end
  
  disp(['[',mfilename,'] Made sure dumps only contain unique (x, z) couples.']);
end

