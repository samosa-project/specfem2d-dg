% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [alpha_vol] = relax_alphavol(freq, RHO, C, MUVOL)
  freq=reshape(freq,[1,max(size(freq))]);
  C=reshape(C,[max(size(C)),1]);
  RHO=reshape(RHO,[max(size(RHO)),1]);
  MUVOL=reshape(MUVOL,[max(size(MUVOL)),1]);
  % Check provided dimensions.
  dimensions=[size(C);size(RHO);size(MUVOL)];
  if(any(dimensions(:,1)/dimensions(1,1)~=1) || any(dimensions(:,2)/dimensions(1,2)~=1))
    disp(['  [',mfilename,', ERROR] Dimensions of inputs:']);
    dimensions
    error(['  [',mfilename,', ERROR] Dimensions of all atmospheric datasets must agree.']);
  end
  omega = 2*pi*freq;
  
  alpha_vol = (MUVOL*omega.^2) ./ (2*RHO.*C.^3);
end

