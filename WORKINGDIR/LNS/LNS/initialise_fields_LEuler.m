function [drho, rho0dvx, rho0dvy, dE] = initialise_fields_LEuler(NX,NY)

  drho = zeros(NX,NY);
  rho0dvx = zeros(NX,NY);
  rho0dvy = zeros(NX,NY);
  dE     = zeros(NX,NY);

  %% Vectorized Navier-Stokes
%   if(~is_backward)
%     gamma_atmos = zeros(NX,NY)+1.4;
%     rho   = zeros(NX,NY)+1.;
%     rhovx = rho*0;
%     rhovy = rho*0;
%     c     = 140;
%     p     = (c.^2).*rho./gamma_atmos;
%     E     = zeros(NX,NY)+p./(gamma_atmos-1.);
%   else
%     gamma_atmos = zeros(NX,NY)+1.4;
%     c     = 140;
%     rho   = zeros(NX,NY);
%     rhovx = zeros(NX,NY);
%     rhovy = zeros(NX,NY);
%     E     = zeros(NX,NY);    
%   end

  %% Vectorized Euler + perturbations
  % rho   = zeros(NX,NY);
  % rhovx = zeros(NX,NY);
  % rhovy = zeros(NX,NY);
  % E     = zeros(NX,NY);

  %% Iterative
  % for i = 1:NX
  %     for j = 1:NY
  %         gamma_atmos(i,j) = 1.4
  %         rho(i,j)   = 1.;
  %         rhovx(i,j) =
  %         rhovy(i,j) =
  %         E(i,j)     =
  %     end
  % end

end