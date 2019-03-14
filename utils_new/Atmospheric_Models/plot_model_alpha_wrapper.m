% Author:        Léo Martire.
% Description:   Wrapper for plot_model_alpha (only needing a model file).
% Notes:         Needs:
%                a) .m scripts and functions (if not alongside this
%                   script, recover via Léo):
%                  1) extract_atmos_model.m
%                  2) plot_model_alpha.m
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [fh] = plot_model_alpha_wrapper(modelfile, freq, zminmax, Nsel)
  [Z, RHO, ~, C, ~, ~, ~, ~, ~, ~, MUVOL, ~, ~, ~, ~, ~, ~, FR, SVIB] = extract_atmos_model(modelfile, 3, 0, 0);
  fh = plot_model_alpha(freq, Z, RHO, C, MUVOL, FR, SVIB, zminmax, Nsel);
end

