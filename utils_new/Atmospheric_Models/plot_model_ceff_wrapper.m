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

function [] = plot_model_ceff_wrapper(DATAFILE, maxz)
  set(0, 'defaultTextFontSize', 20); set(0, 'defaultAxesFontSize', 20);
  if(not(exist('maxz','var')))
    maxz_provided=0;
  else
    maxz_provided=1;
  end
%   nAngles = 200;
  wanted_dz = 10; % [m]

  [Z, ~, T, C, ~, ~, ~, ~, ~, ~, ~, V, U, ~, ~, ~, ~, ~, ~] = extract_atmos_model(DATAFILE, 3, 0, -1);

  sel = 1:floor(wanted_dz/mean(diff(Z))):numel(Z); % subsample
  
  if(not(any(sel)))
    disp(['[',mfilename,', ERROR] sel returned nothing, setting to sel everything.']);
    sel = 1:numel(Z);
  end
  
  % sel=[200,300];
  Z=Z(sel);
  T=T(sel);
  C=C(sel);
  V=V(sel);
  U=U(sel);
  
%   max(Z),pause
  
  if(maxz_provided)
    sel = (Z<=maxz);
    Z=Z(sel);
    T=T(sel);
    C=C(sel);
    V=V(sel);
    U=U(sel);
  end
  
%   max(Z),pause
  
%   angg = linspace(0,360,nAngles+1)*pi/180;
%   angg(end) = [];

  % colzzz = jet(nAngles);
  % colzzz = dayColours(nAngles);
%   colzzz = hsv(nAngles);

%   UNITV = [cos(angg);sin(angg)]';

%   X = zeros(numel(Z), numel(angg), 3);
%   for aa = 1:numel(angg)
%     theta = angg(aa);
%     for zz = 1:numel(Z)
% %       ceff = C(zz) + cos(theta)*WN(zz) - sin(theta)*WE(zz);
%       ceff = C(zz) + uv2azimuth(theta, WE(zz), WN(zz));
%       X(zz, aa, :) = [UNITV(aa,:)*ceff, Z(zz)];
%     end
%   end

%   figure();
%   % ONE = ones(numel(Z),1);
%   ONE = [1,1];
%   for aa = 1:numel(angg)
%     plot3(X(:,aa,1),X(:,aa,2),X(:,aa,3), 'color', colzzz(aa,:),'linewidth',1); hold on;
%     plot3(X(1,aa,1)*ONE,X(1,aa,2)*ONE,[min(Z), max(Z)], 'color', 'k','linewidth',10); hold on;
%   %   scatter3(X(:,aa,1),X(:,aa,2),X(:,aa,3)); hold on;
%   end
%   % figure();
%   % for zz = 1:numel(Z)
%   %   plot3(X(zz,:,1),X(zz,:,2),X(zz,:,3)); hold on;
%   % end
  
%   figure();
%   polarplot(angg, Z);
  
  nAz = 100;
  A = linspace(0,360,nAz+1)*pi/180;
%   A = linspace(0,360,nAz+1);
  [ZZ, AA] = meshgrid(Z, A);
  WIND = uv2azimuth(AA*180/pi, U', V');
%   WIND = uv2azimuth(AA*180/pi, WN', WE');
%   CEFF = repmat(reshape(C,1,numel(C)),[nAz+1,1]) + sin(AA).*WE'+cos(AA).*WN'; name=['$c_\mathrm{eff}$'];
  CEFF = repmat(reshape(C,1,numel(C)),[nAz+1,1]) + WIND; name=['$c_\mathrm{eff}$'];
%   CEFF = repmat(reshape(C,1,numel(C)),[101,1])+cos(AA)*10; % TES PURPOSES
%   CEFF = WIND; name=['wind'];
%   CEFF=CEFF/CEFF(1);
%   figure();plot(max(CEFF), Z);,pause
  
  plot_radial(Z, A, CEFF, name);
  
  spl=split(DATAFILE,'/');
  set(gcf,'name',spl{end});
  customSaveFig(gcf, regexprep(DATAFILE,'\.dat','_ceff'),{'jpg'});
end