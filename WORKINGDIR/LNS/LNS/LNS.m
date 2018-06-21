% Author:        Léo Martire and Quentin Brissaud.
% Mail:          TODO.
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
close all;
clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs.                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NX = 150+1;
NY = 20+1;
% nt = 500;
nt = 100;
DT = 0.01;
IT_IMAGE = 10;
% is_backward = false;
dump_fields = false;
% IT_DUMP  = 10;
DX = 5.;
size_x = NX*DX;
size_y = NY*DX;
x = linspace(0,size_x,NX);
y = linspace(0,size_y,NY);
tt = 0:DT:(nt-1)*DT;
[XX,YY]=meshgrid(x,y);
factor_source=0;
customzlim=[-1,1]*350;

global include_diffusive_terms
include_diffusive_terms = 1; % If 0, linearised Euler. If 1, linearised Navier-Stokes.

% other variables
% x0=1;
% x=linspace(0,5,NX);
% Flux = exp(-(x-x0).^2)+x;
% Flux_analytical = -2*(x-x0).*exp(-(x-x0).^2)+1;
% size_flux = NX;
% DELTAX = x(2)-x(1);
% u = deriv_u_func(DELTAX, size_flux, Flux);
% figure;
% plot(x,u,x,Flux_analytical)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants.                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global g_z mu lambda kappa gamma_atmos soundspeed
g_z = 9.81; % Gravity.
% g_z = 0; % Deactivate gravity.
mu = 1.25e-5; % Dynamic viscosity for dry air at 20 °C.
lambda = (2/3)*mu; % Second viscosity lambda=zeta-(2/3)mu, with zeta=(2/3)mu generally speaking.
kappa = 0.025; % Thermal conductivity for dry air at 20 °C.
gamma_atmos=1.4; % Classical atmospheric $\gamma$.
% soundspeed=349; % Sound speed.
soundspeed=140; % Sound speed. Artificially slow, for CFL, for tests.

global deriv % Handle to derivation function (for more practical tests). 2 & 3 are ok, 3 seems to perform best.
% deriv=@(DX,NX,S,FWD) deriv_u_func(DX, NX, S, FWD)';
% deriv=@(DX,NX,S,FWD) deriv_u_func2(DX, NX, S, FWD)';
deriv=@(DX,NX,S,FWD) deriv_u_func3(DX, NX, S, FWD)';
% deriv=@(DX,NX,S,FWD) deriv_u_func4(DX, NX, S, FWD)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation.                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constitutive variables (allocated by initialisation below).
% drho = zeros(NX,NY);
% rho0dvx = zeros(NX,NY);
% rho0dvy = zeros(NX,NY);
% dE = zeros(NX,NY);

% Auxiliary variables for time scheme.
% aux_drho = zeros(3,NX,NY);
% aux_rho0dvx = zeros(3,NX,NY);
% aux_rho0dvy = zeros(3,NX,NY);
% aux_dE = zeros(3,NX,NY);
aux_drho = zeros(NX,NY);
aux_rho0dvx = zeros(NX,NY);
aux_rho0dvy = zeros(NX,NY);
aux_dE = zeros(NX,NY);

% Order 0 and 1 terms.
RHS_drho = zeros(NX,NY);
RHS_rho0dvx = zeros(NX,NY);
RHS_rho0dvy = zeros(NX,NY);
RHS_dE = zeros(NX,NY);

% Outputs.
sis1 = zeros(nt,1); % Recording station for 1st equation.
sis21 = zeros(nt,1); % Recording station for 2nd equation, 1st term.
sis22 = zeros(nt,1); % Recording station for 2nd equation, 2nd term.
sis3 = zeros(nt,1); % Recording station for 3rd equation.
% kernel_gamma = zeros(NX,NY);
i_seismo = round(length(x)/3); % Recording station x location.
j_seismo = round(length(y)/2); % Recording station y location.

% Integration scheme coefficients.
% % TODO: Those are the 44LRK coefficients, are the 44LSRK coefficients the same?
% rk_a =           [0., 0.5, 0.5, 1.];
% rk_b = 1. / 6. * [1., 2.,  2.,  1.];
% it_scheme_range = 1:4;
% [Carpenter and Kennedy, 1994] Carpenter, M. H. and Kennedy, C. A. (1994). Fourth-order 2n-storage runge-kutta schemes.
rk_a = [0., -567301805773.0/1357537059087.0, -2404267990393.0/2016746695238.0, -3550918686646.0/2091501179385.0, -1275806237668.0/842570457699.0];
rk_b = [1432997174477.0/9575080441755.0, 5161836677717.0/13612068292357.0, 1720146321549.0/2090206949498.0, 3134564353537.0/4481467310338.0, 2277821191437.0/14882151754819.0];
rk_c = [0., 1432997174477.0/9575080441755.0, 2526269341429.0/6820363962896.0, 2006345519317.0/3224310063776.0, 2802321613138.0/2924317926251.0];
it_scheme_range = 1:5;
% % [Berland et al., 2006] Berland, J., Bogey, C., and Bailly, C. (2006). Low-dissipation and low-dispersion fourth-order Runge-Kutta algorithm. Computers & Fluids, 35(10):1459–1463.
% rk_a = [0,             -0.737101392796, -1.634740794341, -0.744739003780, -1.469897351522, -2.813971388035];
% rk_b = [0.03291860514,  0.823256998200,  0.381530948900,  0.200092213184,  1.718581042715,  0.27];
% rk_c = [0,              0.032918605146,  0.249351723343,  0.466911705055,  0.582030414044,  0.847252983783];
% it_scheme_range = 1:6;

% Initial fields.
[drho, rho0dvx, rho0dvy, dE] = initialise_fields_LEuler(NX,NY); % Basically, zeros everywhere.
global RHO0 V0X V0Y E0 P0 % Set initial fields.
RHO0 = 1.4*ones(NX,NY); % Initial density (classical atmospheric ground density).
V0X = 10.*ones(NX,NY); % Initial x-velocity.
V0Y = 0.*ones(NX,NY); % Initial y-velocity.
P0 = (soundspeed.^2).*RHO0./gamma_atmos; % Initial pressure (ideal gas hypothesis).
E0 = P0./(gamma_atmos-1.);

% Source term.
global isource jsource VSTF
VSTF = factor_source*STF(tt); % Compute source time function (STF).
% isource = round(length(x)/1.5); % Source x location.
% jsource = round(length(y)/2.5); % Source y location.
isource = round(length(x)*0.25); % Source x location.
jsource = round(length(y)*0.5); % Source y location.

% Compute derivatives of initial fields.
global DXV0X DXV0Y DYV0X DYV0Y DIVV0
DXV0X=zeros(size(V0X)); % FOR DIFFUSIVE ONLY.
DXV0Y=zeros(size(V0Y)); % FOR DIFFUSIVE ONLY.
DYV0X=zeros(size(V0X)); % CONVECTIVE
DYV0Y=zeros(size(V0Y)); % FOR DIFFUSIVE ONLY
if(include_diffusive_terms)
  for i = 1:NX
    DXV0X(i,:) = deriv(DX, NY, V0X(i,1:NY), 1)'; % FOR DIFFUSIVE ONLY
    DXV0Y(i,:) = deriv(DX, NY, V0Y(i,1:NY), 1)'; % FOR DIFFUSIVE ONLY
  end
end
for j = 1:NY
  DYV0X(:,j) = deriv(DX, NX, V0X(1:NX,j), 1)'; % CONVECTIVE
  if(include_diffusive_terms)
    DYV0Y(:,j) = deriv(DX, NX, V0Y(1:NX,j), 1)'; % FOR DIFFUSIVE ONLY
  end
end
if(include_diffusive_terms)
  DIVV0=DXV0X+DYV0Y; % FOR DIFFUSIVE ONLY
end

% Time range (forward, eventually vs. backward).
it_range     = 1:nt;
% if(is_backward)
%     it_range     = nt:-1:1;
%     it_RK4_range = 4:-1:1;
% end

forward = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time iterations.                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it = it_range
  % Initialise auxiliary variables.
%   aux_drho(2,:,:) = drho(:,:);
%   aux_drho(3,:,:) = drho(:,:);
%   aux_rho0dvx(2,:,:) = rho0dvx(:,:);
%   aux_rho0dvx(3,:,:) = rho0dvx(:,:);
%   aux_rho0dvy(2,:,:) = rho0dvy(:,:);
%   aux_rho0dvy(3,:,:) = rho0dvy(:,:);
%   aux_dE(2,:,:) = dE(:,:);
%   aux_dE(3,:,:) = dE(:,:);
  aux_drho(:,:) = zeros(size(drho));
  aux_rho0dvx(:,:) = zeros(size(drho));
  aux_rho0dvy(:,:) = zeros(size(drho));
  aux_dE(:,:) = zeros(size(drho));

  % Time scheme sub-iterations.
  for subit = it_scheme_range
%     if(it==70)
%       error('kek');
%     end
    % RK4: for each sub-step i, y^{n+1} = y^{n+1} + DT*b_i*k_i, where k_i=L(t_n+DT*c_i,y^{n}+DT*a_{i,i}).
    t = (it-1 + rk_c(subit))*DT; % Local time (t_n+DT*c_i).
    
    % Compute RHS using previous state variable as input and add source terms.
    [RHS_drho, RHS_rho0dvx, RHS_rho0dvy, RHS_dE] = compute_RHS(it, DX, drho(:,:), rho0dvx(:,:), rho0dvy(:,:), dE(:,:), forward);
    
    % Update auxiliary variable.
    % Note: "RHS_*" variables are updated at each iteration and sub-iteration (see below), and are zero at very first iteration (it=1, subit=1).
%     aux_drho(1,:,:) = aux_drho(3,:,:) + DT*rk4a(subit) * reshape(RHS_drho,[1,NX,NY]);
%     aux_rho0dvx(1,:,:) = aux_rho0dvx(3,:,:) + DT*rk4a(subit) * reshape(RHS_rho0dvx,[1,NX,NY]);
%     aux_rho0dvy(1,:,:) = aux_rho0dvy(3,:,:) + DT*rk4a(subit) * reshape(RHS_rho0dvy,[1,NX,NY]);
%     aux_dE(1,:,:) = aux_dE(3,:,:) + DT*rk4a(subit) * reshape(RHS_dE,[1,NX,NY]);
%     aux_drho(1,:,:) = aux_drho(3,:,:) + DT*rk_a(subit) * permute(RHS_drho,[3,1,2]); % Permute is a trick to inject a singleton dimension in order for Matlab to accept to add the NX*NY matrix to the 1*NX*NY matrix.
%     aux_rho0dvx(1,:,:) = aux_rho0dvx(3,:,:) + DT*rk_a(subit) * permute(RHS_rho0dvx,[3,1,2]);
%     aux_rho0dvy(1,:,:) = aux_rho0dvy(3,:,:) + DT*rk_a(subit) * permute(RHS_rho0dvy,[3,1,2]);
%     aux_dE(1,:,:) = aux_dE(3,:,:) + DT*rk_a(subit) * permute(RHS_dE,[3,1,2]);
    aux_drho(:,:) = rk_a(subit)*aux_drho(:,:) + DT*RHS_drho;
    aux_rho0dvx(:,:) = rk_a(subit)*aux_rho0dvx(:,:) + DT*RHS_rho0dvx;
    aux_rho0dvy(:,:) = rk_a(subit)*aux_rho0dvy(:,:) + DT*RHS_rho0dvy;
    aux_dE(:,:) = rk_a(subit)*aux_dE(:,:) + DT*RHS_dE;

    % Apply boundary conditions to auxiliary variables.
%     [aux_drho(1,:,:),aux_rho0dvx(1,:,:),aux_rho0dvy(1,:,:),aux_dE(1,:,:)] = apply_bc(aux_drho(1,:,:),aux_rho0dvx(1,:,:),aux_rho0dvy(1,:,:),aux_dE(1,:,:),gamma_atmos,c,t);
    [aux_drho(:,:),aux_rho0dvx(:,:),aux_rho0dvy(:,:),aux_dE(:,:)] = apply_bc(aux_drho(:,:),aux_rho0dvx(:,:),aux_rho0dvy(:,:),aux_dE(:,:),gamma_atmos,soundspeed,t);

%         if(it_RK4 > 1)
%            forward = ~forward;
%         end

    % Compute RHS using auxiliary variable as input.
%     if(~is_backward)
%       [RHS_drho, RHS_rho0dvx, RHS_rho0dvy, RHS_dE] = compute_RHS(t, DX, NX, NY, drho, rho0dvx, rho0dvy, dE, gamma_atmos, c, forward, isource, jsource);
%     else
%       [RHS_drho, RHS_rho0dvx, RHS_rho0dvy, RHS_dE] = compute_RHS_backward(t, DX, NX, NY, drho, rho0dvx, rho0dvy, dE, gamma_atmos, c, VSTF(:,it), forward,i_seismo,j_seismo);
%     end
%     [RHS_drho, RHS_rho0dvx, RHS_rho0dvy, RHS_dE] = compute_RHS(t, DX, aux_drho(1,:,:), aux_rho0dvx(1,:,:), aux_rho0dvy(1,:,:), aux_dE(1,:,:), gamma_atmos, c, forward, isource, jsource, initial_state);
%     % Add source terms.
%     RHS_dE(isource, jsource) = VSTF(it);
    
    % Update main variable and apply boundary conditions.
%     drho(:,:) = drho(:,:) + DT*rk_b(subit) * RHS_drho;
%     rho0dvx(:,:) = rho0dvx(:,:) + DT*rk_b(subit) * RHS_rho0dvx;
%     rho0dvy(:,:) = rho0dvy(:,:) + DT*rk_b(subit) * RHS_rho0dvy;
%     dE(:,:) = dE(:,:) + DT*rk_b(subit) * RHS_dE;
    drho(:,:) = drho(:,:) + rk_b(subit)*aux_drho;
    rho0dvx(:,:) = rho0dvx(:,:) + rk_b(subit)*aux_rho0dvx;
    rho0dvy(:,:) = rho0dvy(:,:) + rk_b(subit)*aux_rho0dvy;
    dE(:,:) = dE(:,:) + rk_b(subit)*aux_dE;
    [drho,rho0dvx,rho0dvy,dE] = apply_bc(drho,rho0dvx,rho0dvy,dE,gamma_atmos,soundspeed,t);
  end

  % Synthetics.
  sis1(it) = drho(i_seismo, j_seismo);
  sis21(it) = rho0dvx(i_seismo, j_seismo);
  sis22(it) = rho0dvy(i_seismo, j_seismo);
  sis3(it) = dE(i_seismo, j_seismo);
%     sisp(it)  = (gamma_atmos(i_seismo,j_seismo)-1.)*(p(i_seismo,j_seismo)-(1/(2*rho1(i_seismo,j_seismo)))*...
%         (vx(i_seismo,j_seismo)^2+vy(i_seismo,j_seismo)^2));

  % Plot.
  if(mod(it,IT_IMAGE) == 0)
%   if(mod(it,IT_IMAGE) == 0 || abs(it-140)<20)
%     qplot=drho;
%     qplot=rho0dvx;
%     qplot=rho0dvy;
    qplot=dE;
    f=figure(1000);
    if(it==IT_IMAGE)
      f.set('pos',get( 0, 'Screensize' )*0.75);
    end
    
%     pcolor(XX',YY',squeeze(qplot)); % Lighter at execution.
    
    surf(XX',YY',squeeze(qplot)); % Prettier.
%     DDD=get(gca,'Children'); DDD=DDD.ZData;
%     if(max(max(DDD))-min(min(DDD))>0)
%       set(gca, 'dataaspectratio', [1,1,3*(max(max(DDD))-min(min(DDD)))/(max(x)-min(x))]); clear('DDD');
%     end
%     zlim(customzlim);
    
    xlabel('$x$ (m)'); ylabel('$y$ (m)');
    shading flat;
    colorbar;
    axis([min(x), max(x), min(y), max(y)]);
    title(['Iteration ', num2str(it), ', $t=',num2str(t),'$ s']);
    view([-0.25,-0.6,0.25]);
%     rho = 1.;
%     pp     = (c.^2).*rho./gamma_atmos(1);
%     E2     = pp./(gamma_atmos(1)-1.);
% 
%     clear fig gca;
%     figure(1000);
%     fig=surf(XX,YY,p'-E2,'FaceColor','interp');
%     fig.EdgeColor = 'none';
%     hold on;
%     colorbar;
%     plot3(x(isource),y(jsource),0, 'go');
%     plot3(x(i_seismo),y(j_seismo),0, 'ro');
%     saveas(fig,strcat('output_p_',num2str(it),'.png'));
%     close;
  end

%   %% Dump field
%   if(~is_backward && dump_fields && mod(it,IT_DUMP) == 0)
%   dlmwrite(strcat('./field_dump_p_',num2str(it),'.dat'),p,'delimiter','\t');
%   end

%   %% Kernel computation
%   if(mod(it,IT_DUMP) == 0 && is_backward)
%   [kernel_gamma] = compute_kernel(kernel_gamma, rho1, vx, vy, p, it, DT, DX, forward, NX, NY, c, gamma_atmos);
%   end
  
  if(mod(it,10)==0)
    disp(['Iteration ', num2str(it), ', t=', num2str(t), ' s.']);
  end
end

% % Write seismograms
% fid=fopen('./sisvy.dat','w');
% for i=1:length(tt)
%    fprintf(fid,'%12.2f %12.6e %12.6e %12.6e %12.6e \n',...
%        tt(i),sis1(i),sis21(i),sis22(i),sis3(i));
% end
% fclose(fid)

toc