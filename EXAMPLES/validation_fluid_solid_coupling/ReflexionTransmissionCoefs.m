function [R, T] = ReflexionTransmissionCoefs(Z_1, Z_2, theta_1, theta_2)
  R = (Z_2*cos(theta_1) - Z_1*cos(theta_2)) / (Z_2*cos(theta_1) + Z_1*cos(theta_2));
  T = 2*Z_1*cos(theta_1) / (Z_2*cos(theta_1) + Z_1*cos(theta_2));
end