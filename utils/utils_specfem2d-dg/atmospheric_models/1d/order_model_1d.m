% Author:        Léo Martire.
% Description:   Orders quantities as required for the 1D background model files used by the DG extension (FNS and LNS).
% Notes:         Useful only for writing files (to be read by SPECFEM2D-DG).
%                This order should always agree with:
%                - the relevant SPECFEM routine ('define_external_model_DG_only' in 'src/specfem2D/define_external_model.F90'), and
%                - the MSISE-HWM wrapper (MSISE_HWM_wrapper).
%
% Usage:
%   [orderedQuantities] = order_model_1d(Z, RHO, T, C, P, H, G, NBVSQ, KAP, MU, MUVOL, Wnorth, Weast, W, Cp, Cv, GAM, FR, SVIB)
% with:
%   TODO
% yields:
%   TODO

function [orderedQuantities] = order_model_1d(Z, RHO, T, C, P, H, G, NBVSQ, KAP, MU, MUVOL, Wnorth, Weast, W, Cp, Cv, GAM, FR, SVIB)
  mainQuantities = [Z, RHO, T, C, P, H, G, NBVSQ, KAP, MU, MUVOL, Wnorth, Weast, W, Cp, Cv, GAM];
  
  if(not(exist('FR', 'var')) && not(exist('SVIB', 'var')))
    orderedQuantities = mainQuantities;
  
  elseif(exist('FR', 'var') && exist('SVIB', 'var'))
    orderedQuantities = [mainQuantities, FR, SVIB];
  
  else
    error(['[',mfilename,', ERROR] FR and SVIB must be either both given or both not given.']);
  
  end
  
end

