% Author:        Léo Martire.
% Description:   Makes a 1D vector signal 1 peak to peak.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

% Test purposes:
% x=0:.01:3*pi;y=sin(x)+4;figure();plot(x,y); hold on;plot(x,normalise(y))

function [Yn] = normalise(Y, keep)
  if(not(exist('keep')))
    keep=1;
  end
  
%   keep
  
  if(isinteger(keep) || isfloat(keep))
    kept = Y(keep);
  else % isinteger
    if(ischar(keep))
      if(strcmp(keep,'mean'))
        kept = mean(Y);
      else % =='mean'
        error('wrong entry for variable ''keep''');
      end % =='mean'
    else % ischar
      error('wrong entry for variable ''keep''');
    end % ischar
  end % isinteger
  
  Yn = (Y-kept) / (max(Y)-min(Y)) + kept;
end