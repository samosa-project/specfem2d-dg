% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [dx_c_mach] = model_getdxdt(atmospheric_model_file)
  format compact;
  set(0, 'DefaultLineLineWidth', 3); set(0, 'DefaultLineMarkerSize', 8);
  set(0, 'defaultTextFontSize', 14); set(0, 'defaultAxesFontSize', 14);
  set(0, 'DefaultTextInterpreter', 'latex');
  set(0, 'DefaultLegendInterpreter', 'latex');
  
  % Parameters.
  np=2;
  CFL=0.49;
  percentGLL=0.17267316464601141;
  f0 = [];
  while (not(length(f0) == 1))
    f0 = input(['[',mfilename,'] Main temporal frequency [Hz]? > ']);
  end
  paramtxt=['$f_0=',sprintf('%.3e',f0),'$, ',num2str(np),' elts./wavelength'];
  
  % Extract model.
  [Z, ~, ~, C, ~, ~, ~, ~, ~, ~, ~, ~, ~, W, ~, ~, ~] = extract_atmos_model(atmospheric_model_file, 3, 0, 0);
  
  % Compute Mach.
  Mach = W./C;
  
  % Deduce DX, DT.
  dx_c = C./(f0*np);
  dx_c_mach = C.*(1-Mach)/(f0*np);
  dt_c_cfl = CFL*dx_c_mach*percentGLL./C;
  dt_c_mach = (1-Mach)./(np*f0*(1+Mach));
  
  % Figures.
  figure();
  plot(dx_c_mach,Z);
  xlabel('$\Delta x$ [m]');
  ylabel('altitude $z$ [m]');
  title({['maximum acceptable $\Delta x$, accounting for Mach'],['(',paramtxt,')']});
  set(gca,'ticklabelinterpreter','latex');
  grid on;
  
%   figure();
%   plot(dtaccountingformach,Z);
%   xlabel('maximum acceptable $\Delta t$, accounting for Mach number [s]');
%   ylabel('altitude $z$ [m]');
%   title(['$f_0=',sprintf('%.3e',f0),'$']);
%   set(gca,'ticklabelinterpreter','latex');
  
  % Displays.
  disp(['[',mfilename,'] Maximum acceptable DX with this atmospheric model is ',sprintf('%.3e',min(dx_c)),' [m] (c).']);
  disp([blanks(length(mfilename)+2),'    -        -      -   -    -        -       -    -  ',sprintf('%.3e',min(dx_c_mach)),' [m] (c, & Mach).']);
  disp(['[',mfilename,'] Maximum acceptable DT with this atmospheric model is ',sprintf('%.3e',min(dt_c_cfl)),' [s] (Mach-friendly DX, CFL ',sprintf('%.3e',CFL),', & ',sprintf('%.2g',percentGLL*100),'%GLL).']);
  disp([blanks(length(mfilename)+2),'    -        -      -   -    -        -       -    -  ',sprintf('%.3e',min(dt_c_mach)),' [s] (c, & Mach).']);
  
  generatemeshfemparams = [];
  while (not(length(generatemeshfemparams) == 1 && ismember(generatemeshfemparams,[0,1])))
    generatemeshfemparams = input(['[',mfilename,'] Generate meshfem parameters (0 for no, 1 for yes)? > ']);
  end
  
  if(generatemeshfemparams)
    dx=dx_c_mach; % choice

    zmax = [];
    while (not(length(zmax) == 1 && zmax<=max(Z) && zmax>min(Z)))
      zmax = input(['[',mfilename,'] Maximum altitude of interest (',num2str(min(Z)),' < zmax < ',num2str(max(Z)),') [m]? > ']);
    end
    maxlayerz=ceil(length(Z(Z<=zmax))/5);
    if(maxlayerz==1)
      Nlayerz = 1;
    else
      Nlayerz = [];
      while (not(length(Nlayerz) == 1 && Nlayerz<=maxlayerz && Nlayerz>=1))
        Nlayerz = input(['[',mfilename,'] Number of layers (1 <= N <= ',num2str(maxlayerz),') [1]? > ']);
      end
    end
    nx = [];
    while (not(length(nx) == 1))
      nx = input(['[',mfilename,'] nx (as in parfile) [1]? > ']);
    end
    
    IDregionDG=1; % ID of the DG material in parfile.

    interfaces=linspace(min(Z),zmax,Nlayerz+1);

    for l=1:Nlayerz
      bestdx(l)=min(dx(Z>interfaces(l) & Z<interfaces(l+1)));
      nelts(l)=ceil((interfaces(l+1)-interfaces(l))/bestdx(l));
    end
    
    format_layerz=['%',num2str(floor(log10(Nlayerz))+1),'i'];
    format_elts=['%',num2str(floor(log10(sum(nelts)))+1),'i'];
    
    disp(' ');
    disp('# Number of Interfaces:');
    disp(num2str(Nlayerz+1));
    for i=1:Nlayerz+1
      disp(['# Interface n°',num2str(i),':']);
      disp('2');
      disp(['-1d10 ',sprintf('%.15e',interfaces(i))]);
      disp([' 1d10 ',sprintf('%.15e',interfaces(i))]);
    end
    txtlayer={};
    disp('# Layers:');
    for l=1:Nlayerz
      txtlayer{l}=['Air layer n°',sprintf(format_layerz,l),' (from z = ',sprintf('%6.5g',interfaces(l)),' to z = ',sprintf('%6.5g',interfaces(l+1)),', dx = ',sprintf('%4.3g',(interfaces(l+1)-interfaces(l))/nelts(l)),').'];
      disp([num2str(nelts(l)), ' # ',txtlayer{l}]);
    end

    disp(' ');
    disp(['nbregions = ',num2str(Nlayerz)])
    previouslyregisteredelements=0;
    for l=1:Nlayerz
      disp(['1 ',num2str(nx), ' ', sprintf(format_elts,previouslyregisteredelements+1), ' ', sprintf(format_elts,previouslyregisteredelements+nelts(l)), ' ',num2str(IDregionDG), ' # ',txtlayer{l}]);
      previouslyregisteredelements=previouslyregisteredelements+nelts(l);
    end
    
    disp(' ');
    disp(['[',mfilename,'] Copy-paste the two previous blocks, respectively in the interfacefile and the parfile.']);
  end
end