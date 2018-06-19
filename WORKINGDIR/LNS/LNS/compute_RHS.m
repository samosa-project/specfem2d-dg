% function [RHS_rho, RHS_rhovx, RHS_rhovy, RHS_E] = compute_RHS(t, DELTAX, NX, NY, drho1, dvx, dvy, dp, gamma_atmos, c, forward, isource, jsource)
function [RHS_q1, RHS_q21, RHS_q22, RHS_q3] = compute_RHS(it, DX, q1, q21, q22, q3, initial_state, gamma_atmos, c, forward, isource, jsource, VSTF)
  
  % Linearised Euler:
  % q1=rho'
  % q21=rho_0v'_x
  % q22=rho_0v'_y
  % q3=E'
  
  % Handle to derivation function (for more practical tests).
%   deriv=@(DX,NX,S,FWD) deriv_u_func(DX, NX, S, forward)';
%   deriv=@(DX,NX,S,FWD) deriv_u_func2(DX, NX, S, forward)';
  deriv=@(DX,NX,S,FWD) deriv_u_func3(DX, NX, S, forward)';
%   deriv=@(DX,NX,S,FWD) deriv_u_func4(DX, NX, S, forward)';
  % 2 & 3 are ok, 3 seems to perform best.
  
  % Remove eventual singleton dimensions.
  q1=squeeze(q1);
  q21=squeeze(q21);
  q22=squeeze(q22);
  q3=squeeze(q3);
  
  NX=size(q1,1);
  NY=size(q1,2);
  
  % Initialisation.
  flux = zeros(4,NX,NY);
  RHS_q1 = zeros(NX,NY);
  RHS_q21 = zeros(NX,NY);
  RHS_q22 = zeros(NX,NY);
  RHS_q3 = zeros(NX,NY);

%   rho    = 1.;
%   rhocsq = rho*(c)^2;
%   g_z    = 9.81;
  g_z    = 0; % Deactivate gravity.
  mu = 1;
  lambda = 1;
  kappa = 1;
  
  % Linearised Euler (q1=rho', q21=rho_0v'_x, q22=rho_0v'_y, q3=E'): define more practical quantities.
  %%%%%%%%%%%%%%%%%%%%%%%
  % For convective      %
  % terms (\Sigma^c).   %
  %%%%%%%%%%%%%%%%%%%%%%%
  RHO0=squeeze(initial_state(1,:,:)); % Squeeze to remove singleton dimension.
  dvx=q21./RHO0; % RHO0 is positive (>0), this should not cause errors.
  dvy=q22./RHO0; % RHO0 is positive (>0), this should not cause errors.
  V0X=squeeze(initial_state(2,:,:));
  V0Y=squeeze(initial_state(3,:,:));
  E0=squeeze(initial_state(4,:,:));
  P0=(c.^2).*squeeze(initial_state(1,:,:))./gamma_atmos; % Initial pressure (ideal gas hypothesis). % THIS IS MOVEABLE OUTSIDE TIME LOOP.
  P=(gamma_atmos-1)*(E0+q3 - 0.5*(RHO0+q1).*((V0X+dvx).^2+(V0Y+dvy).^2)); % Total pressure from total energy (E0+q3) and total kinetic energy.
  dp=P-P0; % Pressure perturbation, from initial pressure and total pressure.
  %%%%%%%%%%%%%%%%%%%%%%%
  % For diffusive       %
  % terms (\Sigma^d).   %
  %%%%%%%%%%%%%%%%%%%%%%%
  dT = (1/8.3144598) * (P./(RHO0+q1) - P0./RHO0); % Temperature perturbation, from ideal gas law.
  
  %%%%%%%%%%%%%%%%%%%%%%%
  % Derivatives.        %
  %%%%%%%%%%%%%%%%%%%%%%%
  % Compute $\partial_xv_{0,x}$ and $\partial_xv_{0,y}$. % V0 TERMS ARE MOVEABLE OUTSIDE TIME LOOP.
  DXV0X=zeros(size(V0X)); % FOR DIFFUSIVE ONLY.
  DXV0Y=zeros(size(V0Y)); % FOR DIFFUSIVE ONLY.
  DXdvx=zeros(size(dvx)); % FOR DIFFUSIVE ONLY.
  DXdvy=zeros(size(dvy)); % FOR DIFFUSIVE ONLY.
  DXdT=zeros(size(dT)); % FOR DIFFUSIVE ONLY.
  for i = 1:NX
    DXV0X(i,:) = deriv(DX, NY, V0X(i,1:NY), forward)'; % FOR DIFFUSIVE ONLY
    DXV0Y(i,:) = deriv(DX, NY, V0Y(i,1:NY), forward)'; % FOR DIFFUSIVE ONLY
    DXdvx(i,:) = deriv(DX, NY, dvx(i,1:NY), forward)'; % FOR DIFFUSIVE ONLY
    DXdvy(i,:) = deriv(DX, NY, dvy(i,1:NY), forward)'; % FOR DIFFUSIVE ONLY
    DXdT(i,:) = deriv(DX, NY, dT(i,1:NY), forward)'; % FOR DIFFUSIVE ONLY
  end
  % Compute $\partial_y$ terms. % V0 TERMS ARE MOVEABLE OUTSIDE TIME LOOP.
  DYV0X=zeros(size(V0X)); % CONVECTIVE
  DYV0Y=zeros(size(V0Y)); % FOR DIFFUSIVE ONLY
  DYdvx=zeros(size(dvx)); % FOR DIFFUSIVE ONLY
  DYdvy=zeros(size(dvy)); % FOR DIFFUSIVE ONLY
  DYdT=zeros(size(dT)); % FOR DIFFUSIVE ONLY.
  for j = 1:NY
    DYV0X(:,j) = deriv(DX, NX, V0X(1:NX,j), forward)'; % CONVECTIVE
    DYV0Y(:,j) = deriv(DX, NX, V0Y(1:NX,j), forward)'; % FOR DIFFUSIVE ONLY
    DYdvx(:,j) = deriv(DX, NX, dvx(1:NX,j), forward)'; % FOR DIFFUSIVE ONLY
    DYdvy(:,j) = deriv(DX, NX, dvy(1:NX,j), forward)'; % FOR DIFFUSIVE ONLY
    DYdT(:,j) = deriv(DX, NX, dvy(1:NX,j), forward)'; % FOR DIFFUSIVE ONLY
  end
  DIVV0=DXV0X+DYV0Y; % FOR DIFFUSIVE ONLY
  DIVdv=DXdvx+DYdvy; % FOR DIFFUSIVE ONLY
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Divergence terms (1st     %
  % order terms).             %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  DIFFUSIVE_TENSOR_ENERGY_X = mu*( 2*(DXV0X.*dvx+DXdvx.*V0X) + dvy.*(DXV0Y+DYV0X) + V0Y.*(DXdvy+DYdvx) ) ...
                              + lambda*(dvx.*DIVV0 + V0X.*DIVdv) + kappa*DXdT;
  DIFFUSIVE_TENSOR_ENERGY_Y = mu*( 2*(DYV0Y.*dvy+DYdvy.*V0Y) + dvx.*(DXV0Y+DYV0X) + V0X.*(DXdvy+DYdvx) ) ...
                              + lambda*(dvy.*DIVV0 + V0Y.*DIVdv) + kappa*DYdT; % Cf. Brouillons 180619.
  
  % Compute x-axis divergence component (\partial_x(\Sigma^c_x)).
  for j = 2:NY-1
%     gammal_x(1,:,1) = gamma_atmos(1:NX,j);
    size_flux = NX;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Convective terms    %
    % (\Sigma^c).         %
    %%%%%%%%%%%%%%%%%%%%%%%
    % Linearised Euler (q1=rho', q21=rho_0v'_x, q22=rho_0v'_y, q3=E'):
    flux(1,:,j) = RHO0(1:NX,j).*dvx(1:NX,j) + q1(1:NX,j).*V0X(1:NX,j); % $\rho_0v'_x + \rho'v_{0,x}$.
    flux(2,:,j) = V0X(1:NX,j).*q21(1:NX,j) + dp(1:NX,j); % $\rho_0v_{0,x}v'_x + p'$.
    flux(3,:,j) = V0Y(1:NX,j).*q21(1:NX,j); % $\rho_0v_{0,y}v'_x$.
    flux(4,:,j) = (E0(1:NX,j)+P0(1:NX,j)).*dvx(1:NX,j) + (q3(1:NX,j)+dp(1:NX,j)).*V0X(1:NX,j); % $(E_0+p_0)v'_x + (E'+p')v_{0,x}$.
%     % Navier-Stokes
%     flux(1,:,j) = q21(1:NX,j);
%     %     ! Flux for momentum equation
%     flux(2,:,j) = (gammal_x - 1.) .* ( q3(1:NX,j) - (1./(2.*q1(1:NX,j))).*( ...
%     q21(1:NX,j).^2 + q22(1:NX,j).^2 ) ) ...
%     + (q21(1:NX,j).^2)./q1(1:NX,j);
%     flux(3,:,j) = (q21(1:NX,j).*q22(1:NX,j))./q1(1:NX,j);    
%     %     ! Flux for energy equation
%     %     ! Here total energy is "p"
%     flux(4,:,j) = ( gammal_x.*q3(1:NX,j) ...
%     - ((gammal_x - 1.)./(2.*q1(1:NX,j))).*( ...
%     q21(1:NX,j).^2 + q22(1:NX,j).^2 ) ) .* q21(1:NX,j)./q1(1:NX,j);
%     %% Euler + perturbations
%     Flux(1,:,j) = rho*dvx(1:NX,j);
%     %     ! Flux for momentum equation
%     Flux(2,:,j) = dp(1:NX,j)/rho;
%     Flux(3,:,j) = 0;
%     %     ! Flux for energy equation
%     %     ! Here total energy is "p"
%     Flux(4,:,j) = rhocsq*dvx(1:NX,j) + rho*dvy(1:NX,j)*g;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Diffusive terms     %
    % (\Sigma^d).         %
    %%%%%%%%%%%%%%%%%%%%%%%
    flux(2,:,j) = flux(2,:,j) + 2*mu*DXdvx(1:NX,j) + lambda*DIVdv(1:NX,j); % First component of $\Sigma_v(v')_{1,:}$.
    flux(3,:,j) = flux(3,:,j) + mu*(DXdvy(1:NX,j) + DYdvx(1:NX,j)); % Second component of $\Sigma_v(v')_{1,:}$.
    flux(4,:,j) = flux(4,:,j) + DIFFUSIVE_TENSOR_ENERGY_X(1:NX,j); % Cf. Brouillons 180619.
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Derivation.         %
    %%%%%%%%%%%%%%%%%%%%%%%
%     Deriv_x_rho   = deriv(DELTAX, size_flux, Flux(1,:,j), forward);
%     Deriv_x_rhovx = deriv(DELTAX, size_flux, Flux(2,:,j), forward);
%     Deriv_x_rhovy = deriv(DELTAX, size_flux, Flux(3,:,j), forward);
%     Deriv_x_E     = deriv(DELTAX, size_flux, Flux(4,:,j), forward);
    Deriv_x_q1 = deriv(DX, size_flux, flux(1,:,j), forward);
    Deriv_x_q21 = deriv(DX, size_flux, flux(2,:,j), forward);
    Deriv_x_q22 = deriv(DX, size_flux, flux(3,:,j), forward);
    Deriv_x_q3 = deriv(DX, size_flux, flux(4,:,j), forward);
    % Update RHS by adding component. Do not forget the "-" here because the divergence is not classicaly on the RHS.
    RHS_q1(1:NX,j) = RHS_q1(1:NX,j) - Deriv_x_q1';
    RHS_q21(1:NX,j) = RHS_q21(1:NX,j) - Deriv_x_q21';
    RHS_q22(1:NX,j) = RHS_q22(1:NX,j) - Deriv_x_q22';
    RHS_q3(1:NX,j) = RHS_q3(1:NX,j) - Deriv_x_q3';
  end

  % Compute y-axis divergence component (\partial_y(\Sigma^c_y) and \partial_y(\Sigma^d_y)).
  for i = 2:NX-1
%     gammal_y(1,1,:) = gamma_atmos(i,1:NY);
    size_flux = NY;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Convective terms    %
    % (\Sigma^c).         %
    %%%%%%%%%%%%%%%%%%%%%%%
    % Linearised Euler (q1=rho', q21=rho_0v'_x, q22=rho_0v'_y, q3=E'):
    flux(1,i,:) = RHO0(i,1:NY).*dvy(i,1:NY) + q1(i,1:NY).*V0Y(i,1:NY); % $\rho_0v'_y + \rho'v_{0,y}$.
    flux(2,i,:) = V0X(i,1:NY).*q22(i,1:NY); % $\rho_0v_{0,x}v'_y$.
    flux(3,i,:) = V0Y(i,1:NY).*q22(i,1:NY) + dp(i,1:NY); % $\rho_0v_{0,y}v'_y + p'$.
    flux(4,i,:) = (E0(i,1:NY)+P0(i,1:NY)).*dvy(i,1:NY) + (q3(i,1:NY)+dp(i,1:NY)).*V0Y(i,1:NY); % $(E_0+p_0)v'_y + (E'+p')v_{0,y}$.
%     % Navier-Stokes
%     flux(1,i,:) = q22(i,1:NY);
%     %     ! Flux for momentum equation
%     flux(2,i,:) = (q21(i,1:NY).*q22(i,1:NY))./q1(i,1:NY);
%     flux(3,i,:) = (gammal_y - 1.) .* ( q3(i,1:NY) ...
%     - (1./(2.*q1(i,1:NY))).*( q21(i,1:NY).^2 + q22(i,1:NY).^2 ) ) ...
%     + (q22(i,1:NY).^2)./q1(i,1:NY);
%     %     ! Flux for energy equation
%     %     ! Here total energy is "p"
%     flux(4,i,:) = ( gammal_y.*q3(i,1:NY) ...
%     - ((gammal_y - 1.)./(2.*q1(i,1:NY))).*( q21(i,1:NY).^2 ...
%     + q22(i,1:NY).^2 ) ) .* q22(i,1:NY)./q1(i,1:NY);
%     %% Euler + perturbations
%     Flux(1,i,:) = rho*dvy(i,1:NY);
%     %     ! Flux for momentum equation
%     Flux(2,i,:) = 0;
%     Flux(3,i,:) = dp(i,1:NY)/rho;
%     %     ! Flux for energy equation
%     %     ! Here total energy is "p"
%     Flux(4,i,:) = rhocsq*dvy(i,1:NY) + rho*dvy(i,1:NY)*g;

    %%%%%%%%%%%%%%%%%%%%%%%
    % Diffusive terms     %
    % (\Sigma^d).         %
    %%%%%%%%%%%%%%%%%%%%%%%
    flux(2,i,:) = flux(2,i,:) + (mu*(DXdvy(i,1:NY) + DYdvx(i,1:NY))       ); % First component of $\Sigma_v(v')_{2,:}$.
    flux(3,i,:) = flux(3,i,:) + (2*mu*DYdvy(i,1:NY) + lambda*DIVdv(i,1:NY)); % Second component of $\Sigma_v(v')_{2,:}$.
    flux(4,i,:) = flux(4,i,:) + (DIFFUSIVE_TENSOR_ENERGY_Y(i,1:NY)        );
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Derivation.         %
    %%%%%%%%%%%%%%%%%%%%%%%
%     Deriv_y_rho   = deriv(DELTAX, size_flux, Flux(1,i,:), forward);
%     Deriv_y_rhovx = deriv(DELTAX, size_flux, Flux(2,i,:), forward);
%     Deriv_y_rhovy = deriv(DELTAX, size_flux, Flux(3,i,:), forward);
%     Deriv_y_E     = deriv(DELTAX, size_flux, Flux(4,i,:), forward);
    Deriv_y_q1 = deriv(DX, size_flux, flux(1,i,:), forward);
    Deriv_y_q21 = deriv(DX, size_flux, flux(2,i,:), forward);
    Deriv_y_q22 = deriv(DX, size_flux, flux(3,i,:), forward);
    Deriv_y_q3 = deriv(DX, size_flux, flux(4,i,:), forward);
%     % Boundary conditions (maybe?):
%     Deriv_y_q1([1,NY])=0; % Continuity.
%     Deriv_y_q21([1,NY])=0; % Continuity of v_x.
%     Deriv_y_q22([1,NY])=(q21(2)-q21(1))/DX; % Wall for v_y.
%     Deriv_y_q22([NY])=(q21(end)-q21(end))/DX; % Wall for v_y.
%     Deriv_y_q3([1,NY])=0; % Continuity.
    % Update RHS by adding component. Do not forget the "-" here because the divergence is not classicaly on the RHS.
    RHS_q1(i,1:NY) = RHS_q1(i,1:NY) - Deriv_y_q1;
    RHS_q21(i,1:NY) = RHS_q21(i,1:NY) - Deriv_y_q21;
    RHS_q22(i,1:NY) = RHS_q22(i,1:NY) - Deriv_y_q22;
    RHS_q3(i,1:NY) = RHS_q3(i,1:NY) - Deriv_y_q3;
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 0th order terms.          %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Remember to put the "-" if needed.
  % Linearised Euler (q1=rho', q21=rho_0v'_x, q22=rho_0v'_y, q3=E'):
%   RHS_q1 = RHS_q1; % Always zero.
  % x-Momentum:
  % 1st term: zero since $g_x=0$.
  % 2nd term: zero if $\partial_tv_0=0$.
  % 3rd term: only $\rho_0v'_y\partial_yv_{0,x}$ if model is stratified (\partial_x v_{0,x}=0) and if initial velocity is fully horizontal (v_{0,y}=0).
  RHS_q21 = RHS_q21 - q22.*DYV0X;
  % y-Momentum:
  % 1st term: $\rho'g_z$.
  % 2nd term: zero if $\partial_tv_0=0$.
  % 3rd term: zero if initial velocity is fully horizontal (v_{0,y}=0).
  RHS_q22 = RHS_q22 - g_z*q1;
%   RHS_q3 = RHS_q3 + g_z*(q22 + q1.*V0Y); % $g_z(\rho_0v'_z + \rho'v_{0,y})$.
  RHS_q3 = RHS_q3 + g_z*q22; % Only $g_z\rho_0v'_z$ if $g_x=0$ and initial velocity is fully horizontal ($v_{0,y}=0$).
  
  % %% Euler with gravity
  % RHS_rhovy = RHS_rhovy + reshape(drho1(1,:,:),[NX,NY])*g;

  % %% Boundary conditions
  % for j = [1 NY]
  % for i = 1:NX
  %     
  % gammal_y = gamma_atmos(i,j);
  % 
  % if(j==1)
  % dv_dy = (4*(dvy(1,i,j+1)/drho1(1,i,j+1))-(dvy(1,i,j+2)/drho1(1,i,j+2)))/(2*DELTAX);
  % else
  % dv_dy = (4*(dvy(1,i,j-1)/drho1(1,i,j-1))-(dvy(1,i,j-2)/drho1(1,i,j-2)))/(2*DELTAX);
  % end
  % RHS_rho(i,j)   = -drho1(1,i,j)*dv_dy;
  % RHS_rhovx(i,j) = -dvx(1,i,j)*dv_dy;
  % RHS_rhovy(i,j) = -0.;
  % RHS_E(i,j)     = -(gammal_y.*dp(1,i,j) ...
  %                 - ((gammal_y - 1.)./(2.*drho1(1,i,j))).*( ...
  %                         dvx(1,i,j).^2 + dvy(1,i,j).^2 ))*dv_dy;
  % end
  % end
  % 
  % for j = 1:NY
  % for i = [1 NX]
  %     
  % gammal_y = gamma_atmos(i,j);
  %    
  % if(i==1)
  % dv_dx = (4*(dvx(1,i+1,j)/drho1(1,i+1,j))-(dvx(1,i+2,j)/drho1(1,i+2,j)))/(2*DELTAX);
  % else
  % dv_dx = (4*(dvx(1,i-1,j)/drho1(1,i-1,j))-(dvx(1,i-2,j)/drho1(1,i-2,j)))/(2*DELTAX);
  % end
  % RHS_rho(i,j)   = -drho1(1,i,j)*dv_dx;
  % RHS_rhovx(i,j) = -0.;
  % RHS_rhovy(i,j) = -dvy(1,i,j)*dv_dx;
  % RHS_E(i,j)     = -(gammal_y.*dp(1,i,j) ...
  %                 - ((gammal_y - 1.)./(2.*drho1(1,i,j))).*( ...
  %                         dvx(1,i,j).^2 + dvy(1,i,j).^2 ))*dv_dx;
  % end
  % end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Source.                   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Remember to put the "-" if needed.
  RHS_q3(isource, jsource) = RHS_q3(isource, jsource) + VSTF(it);
  
  
  if(false)
    perio = 0.5;
    a  = pi*pi*((1/perio).^2);
    to = 0.5;
    RHS_rhovx(isource,jsource) = RHS_rhovx(isource,jsource) -2.d0*a*(t-to)*exp(-a*(t-to).^2);
    RHS_rhovy(isource,jsource) = RHS_rhovy(isource,jsource) -2.d0*a*(t-to)*exp(-a*(t-to).^2);
  end
  % RHS_E(end/2,end/2) = RHS_E(end/2,end/2) -2.d0*a*(t-to)*exp(-a*(t-to).^2);
  % RHS_E(isource,jsource) = RHS_E(isource,jsource) - exp(-a*(t-to).^2);

  % RHS_matrix(1,:,:) = RHS_rho;
  % RHS_matrix(2,:,:) = RHS_rhovx;
  % RHS_matrix(3,:,:) = RHS_rhovy;
  % RHS_matrix(4,:,:) = RHS_E;
end