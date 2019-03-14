% Author:        Léo Martire.
% Description:   Finds the point where two lines would intersect, given two
%                points on each line, in 2D.
% Notes:         TODO.
%
% Usage:
%   [intersectionPoint] = intersectLines(Line1P1, Line1P2, Line2P1, Line2P2)
% with:
%   TODO.
% yields:
%   TODO.

function [intersectionPoint] = intersectLines(Line1P1, Line1P2, Line2P1, Line2P2)
  
%   x = [0 0; 6 6];  %# Starting points in first row, ending points in second row
%   y = [0 6; 6 0];
  i=1; x = [Line1P1(i) Line2P1(i); Line1P2(i) Line2P2(i)];  %# Starting points in first row, ending points in second row
  i=2; y = [Line1P1(i) Line2P1(i); Line1P2(i) Line2P2(i)];
  dx = diff(x);  %# Take the differences down each column
  dy = diff(y);
  M=[dx;dy];
  den=det(M);
  if(den==0)
    error(['[',mfilename,', ERROR] Lines are parallel. Either there is an infinite number of points in both, or none.']);
  end
  ua = (dx(2)*(y(1)-y(3))-dy(2)*(x(1)-x(3)))/den;
  ub = (dx(1)*(y(1)-y(3))-dy(1)*(x(1)-x(3)))/den;
  xi = x(1)+ua*dx(1);
  yi = y(1)+ua*dy(1);
  intersectionPoint=[xi,yi];
end

