% Author:        Léo Martire.
% Description:   Plots effective sound speed of an atmospheric model file 
%                in many directions.
% Notes:         TODO.
% Usage:
%   plot_model_ceff(atmmodelfile)
% with:
%   TODO.
% yields:
%   TODO.

function [] = plot_model_ceff(DATAFILE)
  nAngles = 200;
  wanted_dz = 1000; % [m]

  [Z, ~, T, C, ~, ~, ~, ~, ~, ~, ~, WN, WE, ~, ~, ~, ~, ~, ~] = extract_atmos_model(DATAFILE, 3, 0, -1);

  sel = 1:floor(wanted_dz/mean(diff(Z))):numel(Z); % subsample for performance
  % sel=[200,300];
  Z=Z(sel);
  T=T(sel);
  C=C(sel);
  WN=WN(sel);
  WE=WE(sel);

  angg = linspace(0,360,nAngles+1)*pi/180;
  angg(end) = [];

  % colzzz = jet(nAngles);
  % colzzz = dayColours(nAngles);
  colzzz = hsv(nAngles);

  UNITV = [cos(angg);sin(angg)]';

  X = zeros(numel(Z), numel(angg), 3);
  for aa = 1:numel(angg)
    theta = angg(aa);
    for zz = 1:numel(Z)
      ceff = C(zz) + cos(theta)*WN(zz) - sin(theta)*WE(zz);
      X(zz, aa, :) = [UNITV(aa,:)*ceff, Z(zz)];
    end
  end

  figure();
  % ONE = ones(numel(Z),1);
  ONE = [1,1];
  for aa = 1:numel(angg)
    plot3(X(:,aa,1),X(:,aa,2),X(:,aa,3), 'color', colzzz(aa,:),'linewidth',1); hold on;
    plot3(X(1,aa,1)*ONE,X(1,aa,2)*ONE,[min(Z), max(Z)], 'color', 'k','linewidth',10); hold on;
  %   scatter3(X(:,aa,1),X(:,aa,2),X(:,aa,3)); hold on;
  end
  % figure();
  % for zz = 1:numel(Z)
  %   plot3(X(zz,:,1),X(zz,:,2),X(zz,:,3)); hold on;
  % end
end