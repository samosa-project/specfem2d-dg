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

function [alpha_vib] = relax_alphavib(freq,c,tau_eps,tau_sig)
  freq=reshape(freq,[1,max(size(freq))]);
  c=reshape(c,[max(size(c)),1]);
  tau_eps=reshape(tau_eps,[max(size(tau_eps)),1]);
  tau_sig=reshape(tau_sig,[max(size(tau_sig)),1]);
  % Check provided dimensions.
  dimensions=[size(c);size(tau_eps);size(tau_sig)];
  if(any(dimensions(:,1)/dimensions(1,1)~=1) || any(dimensions(:,2)/dimensions(1,2)~=1))
    disp(['  [',mfilename,', ERROR] Dimensions of inputs:']);
    dimensions
    error(['  [',mfilename,', ERROR] Dimensions of all atmospheric datasets must agree.']);
  end
%   Nalt=numel(c)
%   Nfreq=numel(freq)
  
%   a1 = (2*pi*pi./c).*(tau_eps-tau_sig);
%   a2 = (freq.^2);
%   a3 = (1+ 2*pi*2*pi*tau_eps.*tau_sig.*freq.^2 );
%   size(a1)
%   size(a2)
%   size(a1*a2)
%   size(a3)
%   alpha_vib = (a1*a2) ./ a3;
%   size(alpha_vib)
%   alpha_vib = ( ((2*pi*pi./c).*(tau_eps-tau_sig)) * (freq.^2)) ./ (1+2*pi*2*pi*tau_eps.*tau_sig.*freq.^2 );
  alpha_vib = ( ((2*pi^2./c).*(tau_eps-tau_sig)) * (freq.^2)) ./ (1 + 4*pi^2 * (tau_eps.*tau_sig) * (freq.^2) );
end

