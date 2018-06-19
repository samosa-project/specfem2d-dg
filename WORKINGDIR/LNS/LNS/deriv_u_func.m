function deriv_u = deriv_u_func(DELTAX, size_flux, Flux, forward)

%% Initialization
RHS_Flux    = zeros(size_flux,1);
Flux_diag   = zeros(size_flux,1);
Flux_diagp1 = zeros(size_flux,1);
Flux_diagm1 = zeros(size_flux,1);
deriv = sparse(size_flux,size_flux);

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     ! Compute boundary terms 
%     ! From eq. (1.7) Gabrro, 2009 part 2 p.13 
    i = 1;
    RHS_Flux(i) = (1/(60*DELTAX)) * ( -137 * Flux(i) ...
                                           +300 * Flux(i+1) ...
                                           -300 * Flux(i+2) ...
                                           +200 * Flux(i+3) ...
                                           -75 * Flux(i+4) ...
                                           +12 * Flux(i+5) );
                                       
                                       %    RHS_Flux(i) = Flux_analytical(i);
   
%    ! From eq. (1.6) Gabrro, 2009 part 2 p.13    
   i = 2;
   RHS_Flux(i)    = (1/(12*DELTAX)) * ( -43 * Flux(i-1) ...
                                                - 80 * Flux(i) ...
                                                + 108 * Flux(i+1) ...
                                                + 16 * Flux(i+2) ...
                                                - Flux(i+3) );
%    ! From eq. (1.7) Gabrro, 2009 part 2 p.13 
   i = size_flux;
   RHS_Flux(i) = (1/(60*DELTAX)) * ( 137 * Flux(i) ...
                                           -300 * Flux(i-1) ...
                                           +300 * Flux(i-2) ...
                                           -200 * Flux(i-3) ...
                                           +75 * Flux(i-4) ...
                                           -12 * Flux(i-5) );                                            
%    RHS_Flux(i) = Flux_analytical(i);

%    ! From eq. (1.6) Gabrro, 2009 part 2 p.13                                              
   i = size_flux-1;
   RHS_Flux(i) = (1/(12*DELTAX)) * ( 43 * Flux(i+1) ...
                                                + 80 * Flux(i) ...
                                                - 108 * Flux(i-1) ...
                                                - 16 * Flux(i-2) ...
                                                + Flux(i-3) );
   
%    ! From eq. (1.5) Gabrro, 2009 part 2 p.13 
   for i = 3:size_flux-2
%         ! Forward coefficients
        if(forward)
        RHS_Flux(i)    = (1/(48*DELTAX)) * ( 5 * Flux(i-2) ...
                                                - 148 * Flux(i-1) ...
                                                + 54 * Flux(i) ...
                                                + 76 * Flux(i+1) ...
                                                + 13 * Flux(i+2) );
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
    Flux_diag(2)        = 8;
    Flux_diag(size_flux-1)     = 8;
    Flux_diag(size_flux)       = 1;
    Flux_diag(3:size_flux-2)   = 3;
    
    Flux_diagm1(1)      = 0;
    Flux_diagm1(2)      = 6;
    Flux_diagm1(size_flux-1)   = 1;
%     Flux_diagm1(size_flux)     = 0;
    Flux_diagm1(size_flux)     = 2;
    Flux_diagm1(3:size_flux-2) = 1;
    
%     Flux_diagp1(1)      = 0;
    Flux_diagp1(1)      = 2;
    Flux_diagp1(2)      = 1;
    Flux_diagp1(size_flux-1)   = 6;
    Flux_diagp1(size_flux)     = 0;
    Flux_diagp1(3:size_flux-2) = 1;
    
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