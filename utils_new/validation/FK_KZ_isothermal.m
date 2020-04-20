function [KZ, corrFact] = FK_KZ_isothermal(Omega_intrinsic, Omega, KX, SOUNDSPEED, H, GRA, GAM, wind_x)
  disp(['[',mfilename,', WARNING] ISOTHERMAL VERIFICATION IS NOT FUNCTIONNAL YET.']);
  disp(['[',mfilename,', WARNING] ISOTHERMAL VERIFICATION IS NOT FUNCTIONNAL YET.']);
  disp(['[',mfilename,', WARNING] ISOTHERMAL VERIFICATION IS NOT FUNCTIONNAL YET.']);
%     % isothdecay = z/H. First try somewhat working. Needs correcting factor 1 (see below).
%     KZ = sqrt( (Omega_intrinsic ./ SOUNDSPEED).^2 - KX.^2 - (1./(2*H.*GAM)).^2 ) + (2*GAM-1)*1j/(2*GAM*H);
%     corrFact=1;

%     % isothDecay = 0. Needs correcting factor 2 (see below).
%     KZ = 0.5 * ( -1j/(H*GAM) + sqrt( 4*( KX.^2 + (Omega_intrinsic/SOUNDSPEED).^2 ) - (H*GAM)^(-2) ) );
%     KZ = 1 / GRA / H / GAM * (-1j * GRA + sqrt(-GRA * (4 * H ^ 2 * GRA * GAM ^ 2 * KX .^ 2 - 4 * H * GAM * wind_x ^ 2 * KX .^ 2 + 8 * H * GAM .* Omega .* wind_x .* KX - 4 * H * GAM * Omega .^ 2 + GRA))) / 2;
%     corrFact=2;

  % isothdecay = z/(2H). Needs correcting factor 3 (see below).
  KZ = (1j * GAM * GRA + -1j * GRA + sqrt(-GRA * (4 * H ^ 2 * GRA * GAM ^ 2 * KX .^ 2 - 4 * H * GAM * wind_x ^ 2 * KX .^ 2 + 8 * H * GAM * Omega .* wind_x .* KX - 4 * H * GAM * Omega .^ 2 + GRA))) / GRA / GAM / H / 2;
  corrFact=3;
end