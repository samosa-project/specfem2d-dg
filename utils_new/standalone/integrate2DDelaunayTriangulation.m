% cf. https://stackoverflow.com/questions/23688669/how-to-integrate-over-a-discrete-2d-surface-in-matlab

function int = integrate2DDelaunayTriangulation(delaunayTri, Z)
  P = delaunayTri.Points;
  T = delaunayTri.ConnectivityList;
  
  d21 = P(T(:,2),:) - P(T(:,1),:);
  d31 = P(T(:,3),:) - P(T(:,1),:);
  
  areas = abs( 1/2 * (d21(:,1) .* d31(:,2) - d21(:,2) .* d31(:,1)));
  
  int = areas'* mean(Z(T), 2);
end