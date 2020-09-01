% Author:        Léo Martire.
% Description:   Encodes some useful parameters for the LNS 2D atmospheric models.
% Notes:         N. A.
%
% Usage:
%   [asc__nb_significantdigits, bin__precision] = bg_model_parameters()

function [asc__nb_significantdigits, bin__precision] = bg_model_parameters()
  asc__nb_significantdigits = 16;
  bin__precision = 'real*8'; % You should not change this, because the Fortran code depends on it.
end

