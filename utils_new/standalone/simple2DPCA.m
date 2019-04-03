% Author:        Léo Martire.
% Description:   Very basic PCA.
% Notes:         TODO.
%
% Usage:
%   [eigenVals, eigenVect] = simple2DPCA(DATA)
% with:
%   TODO.
% yields:
%   eigenVals a vector of the generalized eigenvalues (not sorted),
%   eigenVect a full matrix V whose columns are the corresponding
%             eigenvectors.

% One-liner to produce test data:
%   noi=0.75; rota=pi/6; t=0:.01:8*pi;v=[cos(t);2*sin(t)]'+noi*rand(numel(t),2);R=[cos(-rota),sin(-rota);-sin(-rota),cos(-rota)]; for i=1:numel(t);v(i,:)=(R*v(i,:)')';end;
%   t=linspace(0,2*pi,4); v=[1,0;0,2;-1,0;0,-2];
%   figure();plot(t,v(:,1));hold on;plot(t,v(:,2));
%   figure();plot(v(:,1),v(:,2));daspect([1 1 1]);
%   plot_polarisation(t, detrend(v(:,1)), detrend(v(:,2)));

function [eigenVals, eigenVect] = simple2DPCA(DATA)
  DATA = reshape(DATA, max(size(DATA)), min(size(DATA)));
  d = size(DATA,2);
  if(d ~= 2)
    error(['[',mfilename,', ERROR] Basic PCA only checked for dimension 2 data.']);
  end
  n = size(DATA,1);
  
  MEAN = mean(DATA,1); % barycentre
  
  COV = (1/(n-1)) * (DATA-MEAN)' * (DATA-MEAN); % build covariance matrix
  
  [eigenVect,eigenVals] = eig(COV, 'vector');
%   eigenVals = diag(eigenVals);
  
  % TODO: sort correctly
end

