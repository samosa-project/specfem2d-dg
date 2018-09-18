% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Computes the Richardson number for a given wind and Brunt-Väisälä frequency.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [R] = Richardson_number(WIND, D, N)
  % W (m/s)
  % N (rad/s)
  R = (N ./ (D * WIND)) .^ 2.0;
end