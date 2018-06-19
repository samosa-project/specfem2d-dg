function deriv_u = deriv_u_func4(DELTAX, size_flux, Flux, forward)

%% Initialization
RHS_Flux    = zeros(size_flux,1);
Flux_diag   = zeros(size_flux,1);
Flux_diagp1 = zeros(size_flux,1);
Flux_diagm1 = zeros(size_flux,1);
deriv = zeros(size_flux);

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     ! Compute boundary terms 
%     ! From eq. (1.7) Gabrro, 2009 part 2 p.13 
    i = 1;

    % Second order                                   
    RHS_Flux(i) = (1/(1*DELTAX)) * ( (-5/2) * Flux(i) ...
                                           +2 * Flux(i+1) ...
                                           +(1/2) * Flux(i+2));   %left boundary
                                       
                                       %    RHS_Flux(i) = Flux_analytical(i);
   
%    ! From eq. (1.6) Gabrro, 2009 part 2 p.13    
   i = 2;
   RHS_Flux(i)    = (1/(1*DELTAX)) * ( -(3/4) * Flux(i-1) ...
                                                - 0 * Flux(i) ...
                                                + (3/4) * Flux(i+1) );
%    ! From eq. (1.7) Gabrro, 2009 part 2 p.13 
   i = size_flux;                                          
%    RHS_Flux(i) = Flux_analytical(i);
   
   RHS_Flux(i) = (1/(1*DELTAX)) * ( (5/2) * Flux(i) ...
                                           -2 * Flux(i-1) ...
                                           -(1/2) * Flux(i-2));

%    ! From eq. (1.6) Gabrro, 2009 part 2 p.13                                              
   i = size_flux-1;
   RHS_Flux(i)    = (1/(1*DELTAX)) * ( (3/4) * Flux(i-1) ...
                                                - 0 * Flux(i) ...
                                                - (3/4) * Flux(i+1) );
                                            
alf=7/9/DELTAX;
bet=1/36/DELTAX;
%    ! From eq. (1.5) Gabrro, 2009 part 2 p.13 
   for i = 3:size_flux-2
%         ! Forward coefficients
        if(forward)
        RHS_Flux(i)    = (1/(1)) * ( -bet * Flux(i-2) ...
                                                - alf * Flux(i-1) ...
                                                + 0 * Flux(i) ...
                                                + alf * Flux(i+1) ...
                                                + bet * Flux(i+2) );
        else
        RHS_Flux(i)    = (1/(48*DELTAX)) * ( -13 * Flux(i-2) ...
                                                - 76 * Flux(i-1) ...
                                                - 54 * Flux(i) ...
                                                + 148 * Flux(i+1) ...
                                                - 5 * Flux(i+2) );    
        end
   end
   
%     ! Coefficient from Gabarro, phD thesis, 2009 part 2 p.12/13 
    Flux_diag(1)        = 1;
    Flux_diag(2)        = 1;
    Flux_diag(size_flux-1)     = 1/4;
    Flux_diag(size_flux)       = 1;
    Flux_diag(3:size_flux-2)   = 1;
    
    Flux_diagm1(1)      = 0;
    Flux_diagm1(2)      = 1/4;
    Flux_diagm1(size_flux-1)   = 1/4;
%     Flux_diagm1(size_flux)     = 0;
    Flux_diagm1(size_flux)     = 2;
    Flux_diagm1(3:size_flux-2) = 1/3;
    
%     Flux_diagp1(1)      = 0;
    Flux_diagp1(1)      = 2;
    Flux_diagp1(2)      = 1/4;
    Flux_diagp1(size_flux-1)   = 1/4;
    Flux_diagp1(size_flux)     = 0;
    Flux_diagp1(3:size_flux-2) = 1/3;
    
 %% Build derivation matrix
 for i = 1:size_flux
    deriv(i,i) = Flux_diag(i);
    if(i>1)
       deriv(i,i-1) = Flux_diagm1(i); 
    end
    if(i<size_flux)
       deriv(i,i+1) = Flux_diagp1(i); 
    end
 end
 
 %% Matrix inversion
%  deriv_u = tridag(Flux_diagm1,Flux_diag,Flux_diagp1,RHS_Flux,size_flux);
 deriv_u = deriv\RHS_Flux;

end