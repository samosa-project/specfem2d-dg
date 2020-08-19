
function [R, T] = RTCoefs(fts0_stf1, alpha__1, rho__1, alpha__2, beta__2, rho__2, i_inc)
  % See MotleyToolbox/geophysics/fluid_solid_scattering_matrix.mw
  
  normalise = 0;
  
  if(fts0_stf1)
    % STF
%     i__2 = i_inc; % i_inc is i2 in this case
%     i__1 = snells(alpha__2, alpha__1, i__2); % deduce i1
%     j__2 = snells(alpha__2, beta__2, i__2); % deduce j2
    p = sin(i_inc)/alpha__2; % Compute slowness from solid side.
    D = 8*beta__2^4*p^4*rho__2-4*beta__2^2*p^2*rho__2+rho__1+rho__2; % Intermediate variable, same as in other case.
    
    T = -4*alpha__2*rho__2*(p^2*beta__2^2-1/2)/(alpha__1*D); % transmission
    R(1) = (4*beta__2^2*p^2*rho__2+rho__1-rho__2)/D; % reflection to P-wave
    R(2) = -8*p^2*alpha__2*beta__2*rho__2*(p^2*beta__2^2-1/2)/D; % reflection to S-wave
    
    if(normalise)
      % normalise by sqrt( (outgoing_impedance*cos(outgoing)) / (incoming_impedance*cos(incoming)) )
      T = T * sqrt((rho__1*alpha__1*alpha__1)/(rho__2*alpha__2*alpha__2));
      R(1); % unchanged
      R(2) = R(2) * sqrt((beta__2*beta__2)/(alpha__2*alpha__2));
    end
    
  else
    % FTS
%     i__1 = i_inc; % i_inc is i1 in this case
%     i__2 = snells(alpha__1, alpha__2, i__1); % deduce i2
%     j__2 = snells(alpha__1, beta__2, i__1); % deduce j2
    p = sin(i_inc)/alpha__1; % Compute slowness from fluid side.
    D = 8*beta__2^4*p^4*rho__2-4*beta__2^2*p^2*rho__2+rho__1+rho__2; % Intermediate variable, same as in other case.
    
    T(1) = -4*rho__1*alpha__1*(p^2*beta__2^2-1/2)/(alpha__2*D); % transmission to P-wave
    T(2) = -4*beta__2*p^2*rho__1*alpha__1/D; % transmission to S-wave
    R = (8*beta__2^4*p^4*rho__2-4*beta__2^2*p^2*rho__2-rho__1+rho__2)/D; % reflection
    
    if(normalise)
      error(['[',mfilename,', ERROR] Not implemented.']);
    end
    
  end
end