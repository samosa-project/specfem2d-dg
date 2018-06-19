function [RHS_rho, RHS_rhovx, RHS_rhovy, RHS_E] = compute_RHS_backward(t, DELTAX, NX, NY, drho1, dvx, dvy, dp, gamma_atmos, c, source_term, forward, i_seismo, j_seismo)

%% Initialization
Flux = zeros(4,NX,NY);
RHS_rho   = zeros(NX,NY);
RHS_rhovx = zeros(NX,NY);
RHS_rhovy = zeros(NX,NY);
RHS_E     = zeros(NX,NY);
% RHS_matrix = (4,NX,NY);

rho    = 1.;
rhocsq = rho*(c)^2;
g      = 9.81;
p      = (c^2)*rho/gamma_atmos(1);

%% Compute RHS
for i = 1:NX
    
    gammal_y(1,1,:) = gamma_atmos(i,1:NY);
    size_flux = NY;
    
    %% Navier-Stokes
    Sigma_drho_x = 0;
    Sigma_drho_y = 0;
    Flux(1,i,:) =  Sigma_drho_x*dvx(1,i,1:NY) + Sigma_drho_y*dvy(1,i,1:NY);
%     ! Flux for momentum equation
    Flux(3,i,:) = drho1(1,i,1:NY)+(((gammal_y./(gammal_y-1.)).*p)./rho).*dp(1,i,1:NY);
    Flux(2,i,:) = 0;
%     ! Flux for energy equation
%     ! Here total energy is "p"
    Flux(4,i,:) = (gammal_y-1.).*dvy(1,i,1:NY);
            
    Deriv_y_rho   = deriv_u_func3(DELTAX, size_flux, Flux(1,i,:), forward);
    Deriv_y_rhovx = deriv_u_func3(DELTAX, size_flux, Flux(2,i,:), forward);
    Deriv_y_rhovy = deriv_u_func3(DELTAX, size_flux, Flux(3,i,:), forward);
    Deriv_y_E     = deriv_u_func3(DELTAX, size_flux, Flux(4,i,:), forward);
    
    RHS_rho(i,1:NY)   = RHS_rho(i,1:NY)   - Deriv_y_rho';
    RHS_rhovx(i,1:NY) = RHS_rhovx(i,1:NY) - Deriv_y_rhovx';
    RHS_rhovy(i,1:NY) = RHS_rhovy(i,1:NY) - Deriv_y_rhovy';
    RHS_E(i,1:NY)     = RHS_E(i,1:NY)     - Deriv_y_E';
 
end

for j = 1:NY
    
    gammal_x(1,:,1) = gamma_atmos(1:NX,j);
    size_flux = NX;
    
    %% Navier-Stokes
    Sigma_drho_x = 0.;
    Sigma_drho_y = 0.;
    Flux(1,:,j) =  Sigma_drho_x*dvx(1,1:NX,j) + Sigma_drho_y*dvy(1,1:NX,j);
%     ! Flux for momentum equation
    Flux(2,:,j) = drho1(1,1:NX,j)+(((gammal_x./(gammal_x-1.)).*p)./rho).*dp(1,1:NX,j);
    Flux(3,:,j) = 0;
%     ! Flux for energy equation
%     ! Here total energy is "p"
    Flux(4,:,j) = (gammal_x-1.).*dvx(1,1:NX,j);
 
    Deriv_x_rho   = deriv_u_func3(DELTAX, size_flux, Flux(1,:,j), forward);
    Deriv_x_rhovx = deriv_u_func3(DELTAX, size_flux, Flux(2,:,j), forward);
    Deriv_x_rhovy = deriv_u_func3(DELTAX, size_flux, Flux(3,:,j), forward);
    Deriv_x_E     = deriv_u_func3(DELTAX, size_flux, Flux(4,:,j), forward);
    
    RHS_rho(1:NX,j)   = RHS_rho(1:NX,j)   - Deriv_x_rho;
    RHS_rhovx(1:NX,j) = RHS_rhovx(1:NX,j) - Deriv_x_rhovx;
    RHS_rhovy(1:NX,j) = RHS_rhovy(1:NX,j) - Deriv_x_rhovy;
    RHS_E(1:NX,j)     = RHS_E(1:NX,j)     - Deriv_x_E;
                    
end
 
RHS_rho(i_seismo,j_seismo)   = RHS_rho(i_seismo,j_seismo)   - source_term(1);
RHS_rhovx(i_seismo,j_seismo) = RHS_rhovx(i_seismo,j_seismo) - source_term(2);
RHS_rhovy(i_seismo,j_seismo) = RHS_rhovy(i_seismo,j_seismo) - source_term(3);
RHS_E(i_seismo,j_seismo)     = RHS_E(i_seismo,j_seismo)     - source_term(4);

end