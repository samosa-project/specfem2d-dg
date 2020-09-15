% Prepare ground model.
[s, res] = system([path_to_LITHO,' -p ',sprintf('%.4f',lat),' ',sprintf('%.4f',lon),'']);
res = split(res, newline);
res(end) = [];
model = zeros(1, 6);
c = 1;
for i = numel(res):-1:1
  cur_z_r_vp_vs_qk_qm = str2num(regexp(res{i}, ['( +-?[0-9]+\.[0-9]* ){6}'], 'match', 'once'));
  if(not(isempty(cur_z_r_vp_vs_qk_qm)))
    model(c, :) = cur_z_r_vp_vs_qk_qm;
    c = c + 1;
  end
end

% model(:,1) = - (model(:,1) - min(model(:,1))); % remove negative depths by sliding whole model upward
id_positive_depths = find(model(:,1)>0);
model = model([id_positive_depths(1)-1; id_positive_depths], :); % select only positive depths, and the first layer of negative depth
model(:,1) = -model(:,1); % switch to z-coord instead of depth
model(1, 1) = 0; % set the first layer to be at the interface (topography will take over above here)
% model, pause

model(logical([0; diff(model(:,1))~=0]), :) = []; % make sure no duplicate interface
% model, pause

model(model(:,1)<depthmax, :) = []; % remove layers too deep
% model, pause

model = flipud(model);
itooclose = find(diff(model(:,1))<200);
model(itooclose, :) = []; % find layers too close to each other, remove the top one
% model, pause

NMATERIALS = size(model, 1)+1;

% Deduce interfaces and dx.
% f0 = 2; % [hz]
nptsperwavelength = 2;
% interfaces = [-41140, -31140, -16560, 0, 15e3];
interfaces = [depthmax, model(:,1)', altitudemax];
igrd = find(interfaces==0);
% dx = [2700, 1800, 800, 110, 132];
dx = [(model(:,3)'/f0)/nptsperwavelength, 110, 132];
% dx(1:NMATERIALS-1) = (dx(1:NMATERIALS-1) + 2*dx(NMATERIALS))/3; % reduce dx to allow a nice transition to air

% Print it.
disp('# Number of models.');
disp(['nbmodels                        = ',num2str(NMATERIALS),'']);
disp('#   acoustic:    model_number 1  rho  Vp   0   0   0   QKappa Qmu 0   0   0   0    0    0');
disp('#   elastic:     model_number 1  rho  Vp   Vs  0   0   QKappa Qmu 0   0   0   0    0    0');
% disp('# Remark: for elastic media, Vp and Vs must be the unrelaxed velocities (following the viscoelastic behaviour at infinite frequency). Attenuation can be taken into account by setting the parameters from the ''Attenuation'' section above as wanted.');
% disp('1 1 1.164  349.    0. 0 0   10.   10. 0 0 0 0 0 0 # air');
fmt = '%7.2f';
disp([num2str(1),' 1 ',sprintf(fmt, 1.4),'  ',sprintf(fmt, 340),' ',sprintf(fmt, 0),' 0 0 ',sprintf(fmt, 9999),' ',sprintf(fmt, 9999),' 0 0 0 0 0 0 # Air.']);
disp(['# LITHO1.0 (doi 10.1002/2013JB010626) model computed at (',sprintf('%.4f',lat),'°N, ',sprintf('%.4f',lon),'°E):']);
model = flipud(model);
for i=1:(NMATERIALS-1)
  zminmodel = model(i, 1);
  if(i==NMATERIALS-1)
    zmaxmodel = depthmax;
  else
    zmaxmodel = model(i+1, 1);
  end
  rho = model(i, 2);
  vp = model(i, 3);
  vs = model(i, 4);
  qk = model(i, 5);
  qm = model(i, 6);
  if(qk==0)
    qk = qm-100;
  end
%   rcsq=(vs./vp).^2;
%   qp = ((1-rcsq).*qk.^(-1) + rcsq.*qm.^(-1)).^(-1);
%   qs = qm;
  disp([num2str(i+1),' 1 ',sprintf(fmt, rho),'  ',sprintf(fmt, vp),' ',sprintf(fmt, vs),' 0 0 ',sprintf(fmt, qk),' ',sprintf(fmt, qm),' 0 0 0 0 0 0 # Ground layer, from z=',sprintf('%7.0f',zminmodel),' to z=',sprintf('%7.0f',zmaxmodel),'.']);
%   3 1 2830  6500. 3710. 0 0 2704. 2850. 0 0 0 0 0 0 # middle crust: 16560-31140 (Qp=2750, Qs=2850)
%   4 1 2920  6900. 3930. 0 0 2954. 3100. 0 0 0 0 0 0 # lower crust: 31140-44100 (Qp=3000, Qs=3100)
end
