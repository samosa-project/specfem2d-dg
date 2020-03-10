% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   modelconv_specfem2geoac(spcfm_file, geoac_file)
% with:
%   TODO.
% yields:
%   TODO.

function [] = modelModif_removeAttenuation(spcfm_file, outfile)

  if(nargin<2)
    error(['[',mfilename,', ERROR] Not enough input arguments. Needs ''spcfm_file, outfile''.']);
  end
  
  keep_mu=-1;
  while(not(ismember(keep_mu,[0,1])))
    keep_mu=input(['[', mfilename, '] Keep MU ? (0 for no, 1 for yes) > ']);
  end
  keep_kappa=-1;
  while(not(ismember(keep_kappa,[0,1])))
    keep_kappa=input(['[', mfilename, '] Keep KAPPA? (0 for no, 1 for yes) > ']);
  end
  keep_frsvib=-1;
  while(not(ismember(keep_frsvib,[0,1])))
    keep_frsvib=input(['[', mfilename, '] Keep CO2 attenuation (f_r, S_{vib})? (0 for no, 1 for yes) > ']);
  end
  
  [Z, RHO, T, C, P, H, G, NBVSQ, KAPPA, MU, MUVOL, WN, WE, W, CP, CV, GAMMA, FR, SVIB] = extract_atmos_model(spcfm_file, 3, 0, -1);
  
  if(not(keep_mu))
    MU = 0*MU;
    MUVOL = 0*MUVOL;
  end
  
  if(not(keep_kappa))
    KAPPA = 0*KAPPA;
  end
  
  if(keep_frsvib)
    rewrite_atmos_model(outfile, spcfm_file, Z, RHO, T, C, P, H, G, NBVSQ, KAPPA, MU, MUVOL, WN, WE, W, CP, CV, GAMMA, FR, SVIB);
  else
    rewrite_atmos_model(outfile, spcfm_file, Z, RHO, T, C, P, H, G, NBVSQ, KAPPA, MU, MUVOL, WN, WE, W, CP, CV, GAMMA);
  end
end