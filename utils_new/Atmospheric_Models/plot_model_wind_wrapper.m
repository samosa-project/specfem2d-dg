% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   [fh] = plot_model_wind_wrapper(spcfm_file)
% with:
%   TODO.
% yields:
%   TODO.

function [fh] = plot_model_wind_wrapper(spcfm_file)
  [Z, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, WN, WE, ~, ~, ~, ~] = extract_atmos_model(spcfm_file, 3, 0, 0);
  fh = plot_model_wind(Z, WN, WE);
end