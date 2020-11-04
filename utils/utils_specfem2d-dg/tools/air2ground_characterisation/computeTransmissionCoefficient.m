function [Ta2g_x, Ta2g_z, Tg2a_x, Tg2a_z] = computeTransmissionCoefficient(atmfile,parfile)
  [ZgS,ZgP] = computeImpedanceGround(parfile); % get impedances
  Za = computeImpedanceAir(atmfile); % get impedance
  
  Ta2g_x = 2*ZgS / (Za + ZgS); % transmission from air to ground S-waves
  Ta2g_z = 2*ZgP / (Za + ZgP); % transmission from air to ground P-waves
  
  Tg2a_x = 2*Za / (Za + ZgS); % transmission from ground S-waves to air
  Tg2a_z = 2*Za / (Za + ZgS); % transmission from ground P-waves to air
end