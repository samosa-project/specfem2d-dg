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

function [] = rewrite_atmos_model(output_file, header_file, Z, RHO, T, C, P, H, G, NBVSQ, KAP, MU, MUVOL, Wnorth, Weast, W, Cp, Cv, GAM, FR, SVIB)
  % Check provided dimensions.
  dimensions=[size(Z);size(RHO);size(T);size(C);size(P);size(H);size(G);size(NBVSQ);size(KAP);size(MU);size(MUVOL);size(Wnorth);size(Weast);size(W);size(Cp);size(Cv);size(GAM)];
  if(any(dimensions(:,1)/dimensions(1,1)~=1) || any(dimensions(:,2)/dimensions(1,2)~=1))
    disp(['[',mfilename,', ERROR] Dimensions of inputs:']);
    dimensions
    error(['[',mfilename,', ERROR] Dimensions of all datasets must agree.']);
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
  
  % Informative print.
  disp(['[',mfilename,'] Maximum acceptable DX for CFD computation with this atmospheric model is ',num2str(min(C).*(1-max(W./C))),' / (max(f)*NPointsPerWavelength).']);
  
  % Check if output file exists.
  if(exist(output_file,'file'))
    decision=-1;
    while(not(ismember(decision,[0,1])))
      decision=input([char(strcat("[",mfilename,"] File (", output_file,') exists.\n[',mfilename,'] Overwrite? (0 for no, 1 for yes) >')),' ']);
    end
    if(decision==0)
      disp(['[',mfilename,'] Overwrite cancelled, stopping script.']); return;
    end
  end
  
  % Copy header if it is provided.
  if(isempty(header_file))
    olddatafileprovided=0;
  else
    olddatafileprovided=1;
  end
  if(olddatafileprovided)
    % Open "old" model to read its header.
    f_old = fopen(header_file, 'r');
    if(f_old==-1)
      disp(strcat("[",mfilename,", ERROR] Cannot open old data file ", header_file,').'))
    end
  end
  
  % Open "new" model to write in it.
  f_new = fopen(output_file, 'w');
  if(f_new==-1)
    error(strcat("[",mfilename,", ERROR] Cannot open new data file ", output_file,').'))
  end
  
  % Write header.
  count=1;
  while(count<=3)
    if(olddatafileprovided)
      fprintf(f_new, fgetl(f_old));
    else
      if(count<=2)
        fprintf(f_new,'NO HEADER PROVIDED');
      else
        fprintf(f_new,'z[m] rho[kg/(m^3)] T[K] c[m/s] p[Pa] H[m] g[m/(s^2)] N^2[Hz] kappa[J/(s.m.K)] mu[kg(s.m)] mu_vol[kg/(s.m)] w_M[m/s] w_Z[m/s] w_P[m/s] c_p[J/(mol.K)] c_v[J/(mol.K)] gamma[1]');
        if(frsvibgiven)
          fprintf(f_new,' fr[Hz] Svib[1]');
        end
      end
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