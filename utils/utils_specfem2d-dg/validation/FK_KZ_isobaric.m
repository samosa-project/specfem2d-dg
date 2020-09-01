function [KZ, corrFact] = FK_KZ_isobaric(Omega_intrinsic, Omega, KX, SOUNDSPEED)
  % Test new formula (WORKS QUITE GOOD ON ISOBARIC, WITH Mz = ifftn(filt.*TFMo);)
%     KZ = sqrt( (wind_x.^2 - SOUNDSPEED.^2).*KX.^2 - 2*KX.*Omega + Omega.^2) ./ SOUNDSPEED;
  % Working formula. Nothing else needed.
  KZ = sqrt( (Omega_intrinsic ./ SOUNDSPEED).^2 - KX.^2); corrFact=0;
%     KZ = sqrt(-(SOUNDSPEED * KX - wind_x * KX + Omega) .* (SOUNDSPEED * KX + wind_x * KX - Omega)) / SOUNDSPEED; % Maple raw input works as well.
end