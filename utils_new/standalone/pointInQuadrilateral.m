
% x=[-1 0 1.2 0.5]; y=[0 1 0.2 -2];
% figure(); plot(x,y); hold on; plot(P(1), P(2), '.'); grid on;

function [pointIsInQuadrilateral] = pointInQuadrilateral(x, y, P)
  area_quadrilateral = 0.5*abs(x(1)*y(2)-x(2)*y(1) + x(2)*y(3)-x(3)*y(2) + x(3)*y(4)-x(4)*y(3) + x(4)*y(1)-x(1)*y(4));
  
  area_triangles=zeros(4,1);
  for c=1:4
    A = [x(c), y(c)];
    cp = mod(c, 4)+1;
    B = [x(cp), y(cp)];
%     [c, cp]
    area_triangles(c) = 0.5*abs(A(1)*(B(2)-P(2)) + B(1)*(P(2)-A(2)) + P(1)*(A(2)-B(2)));
  end
  
%   format longg; [sum(area_triangles), area_quadrilateral], sum(area_triangles)-area_quadrilateral
  
  tol = 3*eps;
  if( sum(area_triangles)-area_quadrilateral <= tol)
    pointIsInQuadrilateral = 1;
  elseif( sum(area_triangles)-area_quadrilateral > tol)
    pointIsInQuadrilateral = 0;
  else
    error(['[',mfilename,', ERROR] The sum of the areas of the triangles should not be in any case smaller than the area of the quadrilateral.']);
    pointIsInQuadrilateral = -1;
  end
end