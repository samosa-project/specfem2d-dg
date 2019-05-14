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

function [critAngRad, critAngDeg] = criticalAirToGroundTransmissionAngle(c, vp, vs)
  critAngRad = asin(c / min(vp,vs));
  critAngDeg = critAngRad * 180 / pi;
end

