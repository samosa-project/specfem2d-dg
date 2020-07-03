function [R, T] = ReflexionTransmissionCoefs(Z_1, Z_2, i_1, i_2)
  % Coefficients from Fresnel theory.
  % E.g., Born, M., & Wolf, E. (1981). Principles of Optics: Electromagnetic Theory of Propagation Interference and Diffraction of Light (6th ed.). Cambridge University Press.
  
  R = (Z_2*cos(i_1) - Z_1*cos(i_2)) / (Z_2*cos(i_1) + Z_1*cos(i_2));
  T = 2*Z_1*cos(i_1) / (Z_2*cos(i_1) + Z_1*cos(i_2));
end