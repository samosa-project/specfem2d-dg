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

function sQ=smoooth_n_splinnnn(Z,Q,n)
  spline=fit(Z,smoooth(Q,n),'smoothingspline','smoothingparam',1e-10);
  sQ=spline(Z);
end

function sQ=smoooth(Q,n)
  sQ=smoothdata(Q,'gaussian',n);
%   disp(strcat("[WARNING] Smoothed quantity ", inputname(1), " by a Gaussian filter over ", num2str(n), " elements."));
end