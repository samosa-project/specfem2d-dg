% Author:        Léo Martire.
% Description:   Loads an atmospheric model from a MCD output (Raphaël F.
%                Garcia scripts) file.
% Notes:         TODO.
%
% Usage:
%   [Z, RHO, T, C, P, H, G, NBVSQ, KAP, MU, MUvol, Wnorth, Weast, ...
%    Cp, Cv, GAM, FR, SVIB] = extract_atmos_model_MCDRa(DATAFILE)
% with:
%   TODO.
% yields:
%   TODO.

function [Z, RHO, T, C, P, H, G, NBVSQ, KAP, MU, MUvol, MUvolrot, U, V, ...
          Cp, Cv, GAM, FR, SVIB] = extract_atmos_model_MCDRa(DATAFILE)
  if(not(exist(DATAFILE,'file')))
    error(['[',mfilename,', ERROR] File ''',DATAFILE,''' does not exist. Aborting.']);
  end
  model = importdata(DATAFILE, ' ', 1);
  Z = model.data(:, 1);
  RHO = model.data(:, 2);
  T = model.data(:, 3); % Unused by SPECFEM-DG.
  C = model.data(:, 4); % Only used by SPECFEM-DG for snapshots' background color.
  P = model.data(:, 5);
  H = model.data(:, 6); % Unused by SPECFEM-DG.
  G = model.data(:, 7);
  NBVSQ = model.data(:, 8); % Unused by SPECFEM-DG.
  %=kek.data(:,9); % ncut ?
  KAP = model.data(:, 10);
  MU = model.data(:, 11);
  MUvol = model.data(:, 12); % Unused by SPECFEM-DG.
  MUvolrot = model.data(:,13); % muvolrot?
  FR = model.data(:, 14);
  SVIB = model.data(:, 15);
  U = model.data(:, 16); % meridional (positive towards east)
  V = model.data(:, 17); % zonal (positive toward north)
  Cp = model.data(:, 18); % Unused by SPECFEM-DG.
  Cv = model.data(:, 19); % Used by SPECFEM-DG-FNS for temperature computation.
  GAM = model.data(:, 20);
end

