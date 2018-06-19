function [kernel_gamma] = compute_kernel(kernel_gamma, rho1, vx, vy, p, it, DELTAT, DELTAX, forward, NX, NY, c, gamma_atmos)

rho = 1.;
pp     = (c.^2).*rho./gamma_atmos(1);
E2     = pp./(gamma_atmos(1)-1.);

data_forward = load(strcat('./field_dump_p_',num2str(it),'.dat')) - E2;
data_forward = data_forward';

Flux = zeros(1,NY);
for i = 1:NX
    
    size_flux = NY;
    
    Flux(1,:) =  vy(i,1:NY);
            
    Deriv_x_vx           = deriv_u_func3(DELTAX, size_flux, Flux(1,:), forward);
    kernel_gamma(i,1:NY) = kernel_gamma(i,1:NY) - DELTAT*it*(data_forward(i,1:NY).*Deriv_x_vx');
 
end

Flux = zeros(NX,1);
for j = 1:NY
    
    size_flux = NX;
    
    Flux(:,1) =  vx(1:NX,j);
            
    Deriv_y_vy           = deriv_u_func3(DELTAX, size_flux, Flux(:,1), forward);
    kernel_gamma(1:NX,j) = kernel_gamma(1:NX,j) - DELTAT*it*(data_forward(1:NX,j).*Deriv_y_vy);
 
end

end