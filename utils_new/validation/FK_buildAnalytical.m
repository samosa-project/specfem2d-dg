function [t, x, y, TFM0, KX, KZ] = FK_buildAnalytical(t_0, dt_anal, ...
                                       xmin, xmax, ymin, ymax, dx, dy, ...
                                       zstattab, ...
                                       externalDGAtmosModel, USE_ISOTHERMAL_MODEL, SOUNDSPEED, H, GRA, NSQ, wind_x, wind_y, GAM, ...
                                       TYPE_FORCING, T_0, mult_tSpan, mult_xSpan, mult_ySpan)
  tmax_anal = t_0 + max(zstattab)/SOUNDSPEED(1) + 0.5*T_0 + 1*T_0; % go as far as 2 times the time necessary (roughly) for the signal to pass last station (t_0 + travel time + half period), + full period for safety
  % dt_anal   = T_0 / (2*100); % ensure fmax=1/dt >> 2*f0=2/T_0 by setting 1/dt = (2*2.5)/T_0
  dx_anal = dx;
  dy_anal = dy;
  % syn = zeros(nstat,length(Ztime(nstat,:))) ;
  % synX = zeros(nstat,length(Ztime(nstat,:))) ;
  % synZ = zeros(nstat,length(Ztime(nstat,:))) ;
  % NFFT2 = 2^nextpow2(length(Ztime(nstat,:))); 
  % syn = zeros(nstat, nSamp); %unused?
  % synX = zeros(nstat, nSamp); %unused?
  % synZ = zeros(nstat, nSamp); %unused?
  % NFFT2 = 2^nextpow2(nSamp)*extent_time;
  % NFFT2 = 2^nextpow2(nSamp)*mult_tSpan;
  NFFT2 = 2^nextpow2(ceil(tmax_anal/dt_anal))*mult_tSpan;
  % dt=256*dt; % ?????
  % Fourier space
  NFFT1 = 2^nextpow2((xmax-xmin)/dx_anal)*mult_xSpan;
  NFFT3 = 2^nextpow2((ymax-ymin)/dy_anal)*mult_ySpan;
  disp(['[',mfilename,'] FFT3D matrix size: ',sprintf('%.1e', NFFT1*NFFT2*NFFT3),' reals.']);
  disp(['[',mfilename,'] Usually, 9e7 is ok, 1e8 chugs a bit, 2e8 chugs HARD.']);
  k = zeros(NFFT1, NFFT2, NFFT3);
  % t = zeros(1,NFFT2) ;
  % t = dt * (0:1:NFFT2-1);
  % t = subsampledDt * (0:1:NFFT2-1);
  t = dt_anal * (0:1:NFFT2-1);
  % maxwind = max(t);
  x = dx_anal * (0:1:NFFT1-1) + xmin;
  y = dy_anal * (0:1:NFFT3-1) + ymin;
  % IN THIS SENSE? -> NO!
  %omega = 2.0*pi()*(1.0/(dt*NFFT2))*[ [0.0-[NFFT2/2-1:-1:1]] [0:1:NFFT2/2]];
  %kx = 2.0*pi()*(1.0/(dx*NFFT1))*[[0.0-[NFFT1/2-1:-1:1]] [0:1:NFFT1/2] ];
  % Define the starting Matrix
  %     tendsig=1500.0;
  [X, Y, T] = meshgrid(x, y, t);
  % [X, T] = meshgrid(x, t);
  %     Xp = X+2500;
  %     Yp = Y+2500;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % [Brissaud et al., 2016, Section 5.1]: "Calculation of the forcing signal for the whole time domain along the forcing boundary or at the point source."
  switch(TYPE_FORCING)
    % See expressions in 'boundary_terms_DG.f90'.
    case(1)
  %     Mo = 2 * (   exp(-((T-(t_0-T_0/4)-t0)/(T_0/4)).^2) ...
  %                - exp(-((T-(t_0+T_0/4)-t0)/(T_0/4)).^2) ); % not really pretty
      Mo = -0.5*(2.*(pi/T_0)^2) * (T-t_0) .* exp(-((T-t_0)*pi/T_0).^2); % neary equivalent, and an actual Gaussian derivative
      displayMo = squeeze(Mo(1, 1, :));
    case(2)
      Mo = 2 * (   exp(-((X-(X_0-P_0/4))/(P_0/4)).^2) ...
                 - exp(-((X-(X_0+P_0/4))/(P_0/4)).^2) ) .* ...
               (   exp(-((T-(t_0-T_0/4)-t0)/(T_0/4)).^2) ...
                 - exp(-((T-(t_0+T_0/4)-t0)/(T_0/4)).^2) ); %...
      %              .* exp(-(sqrt((X-xo).^2+(Y-yo).^2)/(lambdo/4)).^2);
    otherwise
      error(['[', mfilename, ',ERROR] Bad TYPE_FORCING.']);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % [Brissaud et al., 2016, Section 5.1]: "Calculation of the 3D (or 2D) Fourier transform (spatial and time transformations) of that function."
  TFM0 = fftn(Mo);
  %%%%
  %%%% IMPORTANT CORRECTION
  %%%%
  % fft2 of matlab perform the projection on the function basis of the type:
  % exp[i(Kx*X+Ky*Y)] whereas we use exp[i(Kx*X-omega*t)] !!!
  % so:
  % * Ky = -omega
  % * and positive frequencies are in the second part of the fft table
  % this is corrected by inserting minus sign in the expression below and by
  % using positive frequencies of the second part of the fft table to ensure
  % symmetry of the fft (see below)
  % the formulas themselves (for dispersion relation) are not modified
  %%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculation of wavenumbers kx = 2π/λx and ky = 2π/λy for all spatial wavelengths λ_{x,y} [Brissaud et al., 2016, Section 5.1].
  % omega = 2.0*pi()*(1.0/(dt*NFFT2))*[[0:1:NFFT2/2] [0.0-[NFFT2/2-1:-1:1]]];
  omega = 2.0*pi()*(1.0/(dt_anal*NFFT2))*[[0:1:NFFT2/2] [0.0-[NFFT2/2-1:-1:1]]];
  kx =    2.0*pi()*(1.0/(dx     *NFFT1))*[[0:1:NFFT1/2] [0.0-[NFFT1/2-1:-1:1]]];
  ky =    2.0*pi()*(1.0/(dy     *NFFT3))*[[0:1:NFFT3/2] [0.0-[NFFT3/2-1:-1:1]]];
  [KX, KY, Omega] = meshgrid(kx, ky, omega);
  Nsqtab = NSQ(1, 1) + 0.0*Omega; % Assuming constant Nsq
  % KX = Omega/velocity(1);
  % see occhipinti 2008 for analytical solution (appendix A1)
  onestab = 0.0*Nsqtab + 1.0;
  Omega_intrinsic = Omega - wind_x*KX - wind_y*KY;

  % Calculation of kz from dispersion relations for all wavenumbers kx, ky and time frequencies [Brissaud et al., 2016, Appendix B].
  % [Brissaud et al., 2016, Eq. (B4)] ?? See brouillons 190412
  % [Fritts and Alexander, 2003, Eq. (22)] ?? Not  ( 1 - Nsqtab./(Omega.^2)) ?? See brouillons 190412
  % Base formula:
  % KZ = sqrt(Nsqtab.*(KX.*KX)./((Omega-(wi)*KX).*(Omega-(wi)*KX))-(KX.*KX) + (1 - Nsqtab./(Omega.*Omega) ).*((Omega-(wi)*KX).*(Omega-(wi)*KX))/(SOUNDSPEED(1)^2) );
  % KZ = sqrt(Nsqtab.*(KX.*KX)./((Omega-(wind_x)*KX).*(Omega-(wind_x)*KX))-(KX.*KX) + (1 - Nsqtab./(Omega.*Omega) ).*((Omega-(wind_x)*KX).*(Omega-(wind_x)*KX))/(SOUNDSPEED(1)^2) );
  % Factorised base formula.
  % KZ = sqrt(   (Nsqtab.*KX.^2) ./ ((Omega-wind_x*KX).^2) ...
  %            - KX.^2 ...
  %            + ( 1 - Nsqtab./(Omega.^2)) ...
  %              .* ((Omega-wind_x*KX)./SOUNDSPEED(1)).^2      );
  % Factorised base formula + generalisation to stick to Fritts and Alexander.
  % KZ = sqrt(   Nsqtab .* (KX.^2+KY.^2) ./ (Omega_intrinsic.^2) ...
  %            - (KX.^2+KY.^2) ...
  %            - onestab/(4*H^2) ...
  %            + ( 1 - Nsqtab./(Omega.^2)) ... % THIS TERM INCREASES LOCAL PROPAGATION SPEED
  %              .* (Omega_intrinsic ./ SOUNDSPEED(1)).^2            );
  % Full Fritts and Alexander (22) with f=0.
  % KZ = sqrt(   (KX.^2+KY.^2) .* (Nsqtab./(Omega_intrinsic.^2) - 1) ...
  %            - onestab/(4*H^2) + (Omega_intrinsic ./ SOUNDSPEED(1)).^2 );
  % Full Fritts and Alexander (22) with f=0 and no H term.
  % KZ = sqrt(   (KX.^2+KY.^2) .* (Nsqtab./(Omega_intrinsic.^2) - 1) ...
  %            + (Omega_intrinsic ./ SOUNDSPEED(1)).^2 );
  % Full Fritts and Alexander (24) with f=0.
  % KZ = sqrt(   (KX.^2+KY.^2) .* (Nsqtab./(Omega_intrinsic.^2) - 1) ...
  %            - onestab/(4*H^2)                                         );
  % Old formula.
  % KZ = sqrt( Nsqtab.*(KX.*KX + KY.*KY)./(Omega_intrinsic.*Omega_intrinsic)-(KX.*KX + KY.*KY) - onestab/(4*H*H) + (Omega_intrinsic.*Omega_intrinsic)/(SOUNDSPEED(1)^2) );
  % Test new formula
  if(externalDGAtmosModel)
    error('no formula for external atmos model');
  else
    if(USE_ISOTHERMAL_MODEL)
      [KZ, corrFact] = FK_KZ_isothermal(Omega_intrinsic, Omega, KX, SOUNDSPEED, H, GRA, GAM, wind_x);
    else
      [KZ, corrFact] = FK_KZ_isobaric(Omega_intrinsic, Omega, KX, SOUNDSPEED);
    end
  end
  indNanKZ = find(isnan(KZ));
  indInfKZ = find(isinf(KZ));
  indImagKZ = find(imag(KZ)<0);
  disp(['[',mfilename,'] Number of NaNs in KZ: ',num2str(numel(indNanKZ)),' (',sprintf('.0f',(numel(indNanKZ)/numel(KZ))*100),'%).']);
  if(numel(indNanKZ))
    disp(['[',mfilename,']   Setting those NaNs to zeros.']);
    KZ(indNanKZ) = 0.0;
  end
  disp(['[',mfilename,'] Number of Infs in KZ: ',num2str(numel(indInfKZ)),' (',sprintf('.0f',(numel(indInfKZ)/numel(KZ))*100),'%).']);
  if(numel(indInfKZ))
    disp(['[',mfilename,']   Setting those Infs to zeros.']);
    KZ(indInfKZ) = 0.0;
  end
  disp(['[',mfilename,'] Number of KZ such that Im(KZ)<0: ',num2str(numel(indImagKZ)),' (',sprintf('.0f',(numel(indImagKZ)/numel(KZ))*100),'%).']);
  if(numel(indImagKZ))
    disp(['[',mfilename,']   Imaginary part of KZ should be positive in order to attenuate the signal. Setting those to their conjugate (only flips the sign of the imaginary part).']);
    KZ(indImagKZ) = conj(KZ(indImagKZ));
  end
  % real(KZ) should be positive for positive frequencies and negative for
  % negative frequencies in order to shift signal in positive times
  % restore the sign of KZ depending on Omega-wi*KX
  %     KZnew=real(KZ).*sign((Omega-wind_x*KX)).*sign(KX)+1i*imag(KZ);
  % !!! Why KZ should have a sign opposite to Omega for GW NOT UNDERSTOOD !!!
  % => because vg perpendicular to Vphi ?
%   KZnew = 0.0 - real(KZ).*sign(Omega_intrinsic) + 1i*imag(KZ); KZ = KZnew;
  KZnew = 0.0 - real(KZ).*sign(Omega_intrinsic); KZ = KZnew;

  %%%%%%%%%%%%%%%%%%%%%%%
  % CORRECTING FACTORS. %
  %%%%%%%%%%%%%%%%%%%%%%%
  switch corrFact
    case 1
      % Correcting factor 1 for isothermal case (use it with the Maple script configured with the isothermal decay = z/H).
      KZ = KZ -1j*(1/H + 2/(SOUNDSPEED^2));
      disp(['[',mfilename,'] Used a correcting factor: kz = kz - i*( 1/H+2/(c^2) )']);
      % 1/H comes from the fact that isothDecay=1/H (developing K \cdot X makes this 1/H go outside the i and finds back its place). TBH I do not remember how I found the second factor, sheer luck is not excluded.
    case 2
      % Correcting factor 2 for isothermal case (use it with the Maple script configured without the isothermal decay).
      KZ = KZ -1j*(8.589e-05);
      disp(['[',mfilename,'] Used a correcting factor: kz = kz - i*8.589e-5']);
      % Cannot explain why this value.
    case 3
      % Correcting factor 3 for isothermal case (use it with the Maple script configured with isothermal decay = z/(2H)).
  %     KZ = KZ -1j*(6.449e-05);
      KZ = KZ -1j*(1/(2*H) + 2/(SOUNDSPEED^2));
      disp(['[',mfilename,'] Used a correcting factor: kz = kz - i*( 1/(2H)+2/(c^2) )']);
      % Same rationale as case 1.
  end
  % KZ = KZ -1j*(2/(SOUNDSPEED^2)); % Test correcting factor for isothermal case.
  % KZ = KZ -1j*(1.1643e+05); % Test correcting factor for isothermal case.
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%

  % restore negative frequencies from postive ones:
  %      KZ(:,NFFT2/2+2:NFFT2)=0.0-real(KZ(:,NFFT2/2:-1:2))+1i*imag(KZ(:,NFFT2/2:-1:2));
  % restore negative frequencies from postive ones:
  % Corrected for "%%%%%%    IMPORTANT CORRECTION   %%%%"
  %      KZ(:,NFFT2/2:-1:2)=0.0-real(KZ(:,NFFT2/2+2:NFFT2))+1i*imag(KZ(:,NFFT2/2+2:NFFT2));
  % Above not taken into account for "filt", because only positive computations 
end