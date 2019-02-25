% Author:        Léo Martire.
% Description:   Computes tilt from vertical displacement at multiple
%                points in space.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function tilt = computeTilt(x_l, x_r, t_l, t_r, vz_l, vz_r)
  Nt = size(t_l,2);
  N = size(vz_l, 1);
  Nd = size(vz_r, 1);
  if(N~=Nd)
    error(['[',mfilename,'] You must provide the same number of left stations as right stations.']);
  end
  clear('Nd');
  
%   disp(['[',mfilename,'] Computing displacement by numerical integration.']);
  displ_l = zeros([N, Nt]);
  displ_r = zeros([N, Nt]);
  tilt = zeros([N, Nt]);
  for i=1:N
    displ_l(i,:) = cumtrapz(t_l(i, :), vz_l(i, :));
    displ_r(i,:) = cumtrapz(t_r(i, :), vz_r(i, :));
    tilt(i,:) = atan((displ_r(i, :) - displ_l(i, :)) / (x_r(i) - x_l(i)));
  end
end