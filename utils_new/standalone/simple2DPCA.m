% Author:        Léo Martire.
% Description:   Very basic PCA.
% Notes:         TODO.
%
% Usage:
%   [eigenVals, eigenVect] = simple2DPCA(DATA)
% with:
%   TODO.
% yields:
%   TODO.

function [eigenVals, eigenVect] = simple2DPCA(DATA)
  DATA = reshape(DATA,max(size(DATA),min(size(DATA))));
  d=size(DATA,2);
  if(d~=2)
    error(['[',mfilename,', ERROR] Basic PCA only implemented for dimension 2 data.']);
  end
  
  b = -(sum(DATA(:, 1).^2) - sum(DATA(:, 2).^2)) / (2 * sum(DATA(:,1) .* DATA(:, 2)));
  eigenVect(:, 1) = [1, b + sqrt(b^2 + 1)];
  eigenVect(:, 2) = [1, b - sqrt(b^2 + 1)];
  
  for i=1:d
    eigenVect(:,i) = eigenVect(:,i)/norm(eigenVect(:,i));
  end
  
  eigenVals = [-1,-1];
%   disp(['[',mfilename,'] Eigenvalues not implemented yet.']);
end

