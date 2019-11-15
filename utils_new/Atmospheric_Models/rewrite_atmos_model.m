% Author:        Léo Martire.
% Description:   Outputs an atmospheric model file to the SPECFEM-DG
%                'external_DG' format.
% Notes:         N/A.
%
% Usage 1:
%   rewrite_atmos_model(output_file, header_file, Z, RHO, T, C, P, H, ...
%                       G, NBVSQ, KAP, MU, MUVOL, Wnorth, Weast, W, ...
%                       Cp, Cv, GAM)
% with:
%   output_file an output file path,
%   header_file a file containing a three-line header,
%   Z           an altitude array [m],
%   RHO         a density array [kg/m^3],
%   T           a temperature array [K] (not used by SPECFEM-DG),
%   C           a speed of sound array [m/s] (only used by SPECFEM-DG for
%               snapshots' background color),
%   P           a pressure array [Pa],
%   H           a scale height array [m] (not used by SPECFEM-DG),
%   G           a gravity array [m/s^2],
%   NBVSQ       a Brunt-Vaisälä frequency array [Hz] (not used by
%               SPECFEM-DG),
%   KAP         a thermal conductivity array [J/(s.m.K)],
%   MU          a viscosity array [kg/(s.m)],
%   MUVOL       a volumic viscosity array [kg/(s.m)] (not used by
%               SPECFEM-DG),
%   Wnorth      a meridional wind array [m/s],
%   Weast       a zonal wind array [m/s],
%   W           a wind array [m/s],
%   Cp          an isobaric specfic heat capacity array [J/(mol.K)] (not
%               used by SPECFEM-DG),
%   Cv          an isochoric specific heat capacity array [J/(mol.K)]
%               (used by SPECFEM-DG-FNS for temperature computation),
%   GAM         an adiabatic ratio array [1],
% outputs the model to the given file.
%
% Usage 2:
%   rewrite_atmos_model(output_file, header_file, Z, RHO, T, C, P, H, ...
%                       G, NBVSQ, KAP, MU, MUVOL, Wnorth, Weast, W, ...
%                       Cp, Cv, GAM, FR, SVIB)
% with:
%   FR          a relaxation frequency array,
%   SVIB        a relaxation strength array,
% does the same thing as Usage 1, but adds the relaxation parameters to
% the file.

function [] = rewrite_atmos_model(output_file, HF, Z, RHO, T, C, P, H, G, NBVSQ, KAP, MU, MUVOL, Wnorth, Weast, W, Cp, Cv, GAM, FR, SVIB, ask0_overwrite1)
  % parameters
  defaultColumnHead = 'z[m] rho[kg/(m^3)] T[K] c[m/s] p[Pa] H[m] g[m/(s^2)] N^2[Hz] kappa[J/(s.m.K)] mu[kg(s.m)] mu_vol[kg/(s.m)] w_M[m/s] w_Z[m/s] w_P[m/s] c_p[J/(mol.K)] c_v[J/(mol.K)] gamma[1]';
  headersizeneededifnumericalarray = 10; % [year, dayofyear, secondsofday, lat, lon, LST, F107A, F107, AP, SWITCHPLANET]
  SWITCHPLANET_acceptable = [3, 4]; % earth (3) and mars (4)
  
  % Check provided dimensions.
  dimensions=[size(Z);size(RHO);size(T);size(C);size(P);size(H);size(G);size(NBVSQ);size(KAP);size(MU);size(MUVOL);size(Wnorth);size(Weast);size(W);size(Cp);size(Cv);size(GAM)];
  if(any(dimensions(:,1)/dimensions(1,1)~=1) || any(dimensions(:,2)/dimensions(1,2)~=1))
    disp(['[',mfilename,', ERROR] Dimensions of inputs:']);
    dimensions
    error(['[',mfilename,', ERROR] Dimensions of all datasets must agree.']);
  end
  % treat header input
  if(isempty(HF))
    % test if empty
    HF_empty0_file1_array2 = 0;
  else
    % test if char array
    if(ischar(HF))
      % test if it is a path to a file
      if(exist(HF,'file'))
        HF_empty0_file1_array2 = 1;
      else
        error(['header_file is a character array, but does not point at an existing file (''',HF,''' does not exist)']);
      end
    else
      % test if numerical array
%       if( (isfloat(header_file) || isinteger(header_file)) && numel(header_file)==headersizeneededifnumericalarray && ismember(header_file(10),SWITCHPLANET_acceptable))
      if(isstruct(HF))
        HF_empty0_file1_array2 = 2;
      else
        error(['header_file must be either [], a character array, or a structure (with at least fields PLANET, lat, and lon, and others for some planets']);
      end
    end
  end
  % Treat eventual FR/SVIB presence in arguments' list.
  if(exist('FR') && exist('SVIB'))
    disp(['[',mfilename,'] FR and SVIB given. They will be added to atmospheric model file.']);
    frsvibgiven=1;
  else
    frsvibgiven=0;
  end
  if(exist('FR') && any(size(Z)~=size(FR)))
    error(['[',mfilename,', ERROR] Dimensions of FR must agree with other dimensions.']);
  end
  if(exist('SVIB') && any(size(Z)~=size(SVIB)))
    error(['[',mfilename,', ERROR] Dimensions of SVIB must agree with other dimensions.']);
  end
  if(xor(~exist('FR'), ~exist('SVIB')))
    error(['[',mfilename,', ERROR] FR and SVIB must be given together.']);
  end
  
  if(not(exist('ask0_overwrite1')))
    ask0_overwrite1=0;
  end
  
  % Informative print.
  disp(['[',mfilename,'] Maximum acceptable DX for CFD computation with this atmospheric model is ',num2str(min(C).*(1-max(W./C))),' / (min(c)*(1-max(Mach))).']);
  
  % Check if output file exists.
  if(exist(output_file,'file'))
    if(ask0_overwrite1)
      % overwrite
      decision = 1;
      disp(['[', mfilename, ', INFO] Overwriting existing file.']);
    else
      % else, ask
      decision=-1;
      while(not(ismember(decision,[0,1])))
        decision=input([char(strcat("[", mfilename, "] File (", output_file,') exists.\n[',mfilename,'] Overwrite? (0 for no, 1 for yes) >')),' ']);
      end
    end
    
    if(decision==0)
      disp(['[', mfilename, '] Overwrite cancelled, stopping script.']); return;
    end
  end
  
  switch(HF_empty0_file1_array2)
    case 0
      % do nothing (will be prepared later)
    case 1
      % Open "old" model to read its header.
      f_old = fopen(HF, 'r');
      if(f_old==-1)
        disp(strcat("[",mfilename,", ERROR] Cannot open old data file ", HF,').'))
      end
    case 2
      % extract from array
      %   IYD=year*1000.0+day from msiswrapper
      %[year, dayofyear, secondsofday, lat, lon, LST, F107A, F107, AP]
%       HF_year = HF(1);
%       HF_day = HF(2);
%       HF_sec = HF(3); 
%       HF_lat = HF(4); HF_format_latlon = '%.6f';
%       HF_lon = HF(5);
      
%       HF_glat = HF_lat;
%       HF_glon = HF_lon;
%       HF_LST = HF(6); HF_format_LST = '%.7f';
%       HF_F107A = HF(7); HF_format_solarparams = '%.6f';
%       HF_F107 = HF(8);
%       HF_AP = HF(9);
%       HF_planetstr=num2str(HF(10));
      
      HF_format_float = '%.6f';
      
      headerline{1} = ['PLANET=',num2str(HF.planet),' lat ',sprintf(HF_format_float,HF.lat),' lon ',sprintf(HF_format_float,HF.lon),' '];
      
      switch(HF.planet)
        case 3
          % earth
          f = 1./298.257223563; % from msisehwm_wrapper, flattening of Earth deduced from the WGS84 value of 1/f. Used to compute the geodetic latitude (GLAT).
          HF_glat = atan( tan(HF.lat*pi/180.) / ((1.-f)^2) ) * 180./pi; % from msisehwm_wrapper
          HF_glon = HF_lon;
          HF_IYD = HF.year*1000.0+HF.doy; % from msisehwm_wrapper
          HF_LST=HF.sod/3600.+HF_glon/15.; % from msisehwm_wrapper
          if(HF_LST>24.0)
            HF_LST=HF_LST-24.0;
          end
          headerline{1} = [headerline{1}, 'year ',sprintf('%d',HF.year),' day ',sprintf('%d',HF.doy),' seconds ',sprintf(HF_format_float,HF.sod)];
          headerline{2} = [' IYD ',sprintf('%d',HF_IYD),' SEC ',sprintf(HF_format_float,HF.sod),...
                           ' GLAT ',sprintf(HF_format_float,HF_glat),' GLON ',sprintf(HF_format_float,HF_glon),...
                           ' STL ',sprintf(HF_format_LST,HF_LST),' F107A ',sprintf(HF_format_solarparams,HF.f107a),...
                           ' F107 ',sprintf(HF_format_solarparams,HF.f107),' AP ',sprintf(HF_format_solarparams,HF.ap)];
        case 4
          % mars
          headerline{1} = [headerline{1}, 'LS ',sprintf(HF_format_float,HF.ls),' LT ',sprintf(HF_format_float,HF.lt)];
          headerline{2} = [''];
        otherwise
          error();
      end
      if(frsvibgiven)
        headerline{3} = [defaultColumnHead, ' fr[Hz] Svib[1]'];
      else
        headerline{3} = defaultColumnHead;
      end
    otherwise
      error('error, and script should not be able to go down this far');
  end
  
  % Open "new" model to write in it.
  f_new = fopen(output_file, 'w');
  if(f_new==-1)
    error(strcat("[",mfilename,", ERROR] Cannot open new data file ", output_file,').'))
  end
  
  % Write header.
  count=1;
  while(count<=3)
    switch(HF_empty0_file1_array2)
      % emtpy header
      case 0
        if(count<=2) % if first two lines, print this
          fprintf(f_new,'NO HEADER PROVIDED');
        else % if last line of header, print comlum head
          fprintf(f_new,defaultColumnHead);
          if(frsvibgiven)
            fprintf(f_new,' fr[Hz] Svib[1]');
          end
        end
      case 1
        % file exists, copy it
        fprintf(f_new, fgetl(f_old));
      case 2
        % print header from array
        fprintf(f_new, headerline{count});
      otherwise
        error('error, and script should not be able to go down this far');
    end
    fprintf(f_new, "\n");
    count=count+1;
  end
  
%z[m] rho[kg/(m^3)] T[K] c[m/s] p[Pa] H[m] g[m/(s^2)] N^2[1/s] kappa[J/(s.m.K)] mu[kg(s.m)] mu_vol[kg/(s.m)] w_M[m/s] w_Z[m/s] w_P[m/s] c_p[J/(mol.K)] c_v[J/(mol.K)] gamma[1] fr[Hz] Svib[1]
  for i=1:numel(Z)
    line_numbers=[Z(i), RHO(i), T(i), C(i), P(i), H(i), G(i), NBVSQ(i), KAP(i), MU(i), MUVOL(i), Wnorth(i), Weast(i), W(i), Cp(i), Cv(i), GAM(i)];
    fprintf(f_new,'%15.8e ',line_numbers);
    if(frsvibgiven)
      frsvib_line_numbers=[FR(i),SVIB(i)];
      fprintf(f_new,'%15.8e ',frsvib_line_numbers);
    end
    fprintf(f_new, "\n");
  end
  
  fclose('all');
  disp(strcat("[",mfilename,"] Model output to file: '", output_file, "'."));
end

function checknonsense(output_file, header_file, Z, RHO, T, C, P, H, G, NBVSQ, KAP, MU, MUVOL, Wnorth, Weast, W, Cp, Cv, GAM, FR, SVIB)
  if(any(RHO<=0))
    error();
  end
end