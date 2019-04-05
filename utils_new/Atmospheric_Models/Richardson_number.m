% Author:        Léo Martire.
% Description:   Computes the Richardson number for a given wind and
%                Brunt-Väisälä frequency.
% Notes:         The Richardson number should be >1 - or at least >0.25 -
%                to prevent fluid instabilities.
%
% Usage:
%   [R] = Richardson_number(WIND, D, N)
% with:
%   WIND a wind array,
%   D    a differentiation matrix (typically, see
%        differentiation_matrix.m),
%   N    a Brunt-Väisälä frequency array (/!\: N, not N^2).
% yields:
%   R    the Richardson number array.

function [R] = Richardson_number(WIND, D, N)
  % W (m/s)
  % N (rad/s)
  R = (N ./ (D * WIND)) .^ 2.0;
end