% Author:        Léo Martire.
% Description:   Produces a nice color scale for times of day.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function colours=dayColours(N)
  if(not(exist('N')))
   N=24;
  end

  midnight = [0,1,1];
  midmorning = [0,0.5,0.4];
  midday = [1,1,0];
  midafternoon = [1,0,1];
  
  RGBday=[midnight;midmorning;midday;midafternoon];
  colours = zeros(N, 3);
  for i=1:3
    colours(:, i) = interp1([0,6,12,18,24], [RGBday(:, i); RGBday(1, i)],linspace(0,23,N), 'linear');
  end
end