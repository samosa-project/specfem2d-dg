% Author:        Léo Martire.
% Description:   Parses a parfile, deduce optimal DZ from user-provided
%                thicknesses, a main frequency, and a number of points per
%                wavelength.
% Notes:         Needs:
%                a) .m scripts and functions (if not alongside this
%                   script, recover via Léo):
%                  1) extract_atmos_model.m
%                  2) model_viscoelastic_getinterfaces.m (only if solid
%                     models are to be read)
%
% Usage:
%   model_getdxdt(atmospheric_model_file)
% with:
%   atmospheric_model_file the path to an atmospheric model file.

function model_getdxdt()
  addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new'); % for model_viscoelastic_getinterfaces.m
  
  % Parameters.
  np=2;
  CFL=0.49;
  percentGLL=0.17267316464601141;
  f0 = [];
  IDregionDG=1; % ID of the DG material in parfile.
  plot_fluidmodel_dz=1;
  
  fluid_model_type='"external-DG" fluid';
  solid_model_type='"parfile" viscoelastic';
  
  % Ask user.
  fluidmodel = [];
  while (not(length(fluidmodel) == 1 && ismember(fluidmodel,[0,1])))
    fluidmodel = input(['[',mfilename,'] ', fluid_model_type, '    model to be loaded (0 for no, 1 for yes)? > ']);
  end
  solidmodel = [];
  while (not(length(solidmodel) == 1 && ismember(solidmodel,[0,1])))
    solidmodel = input(['[',mfilename,'] ', solid_model_type, ' model to be loaded (0 for no, 1 for yes)? > ']);
  end
  if(not(any([fluidmodel, solidmodel])))
    error(['[',mfilename,', ERROR] Neither external-DG atmospheric model nor parfile viscoelastic model. Nothing can be done here.']);
  end
  
  % Ask user for dominant frequency.
  while (not(length(f0) == 1))
    f0 = input(['[',mfilename,'] Main temporal frequency [Hz]? > ']);
  end
  paramtxt=['$f_0=',sprintf('%.3e',f0),'$, ',num2str(np),' elts./wavelength'];
  
  % Load
  if(fluidmodel)
    % Ask for file.
    atmospheric_model_file=input(['[',mfilename,'] Path to ',fluid_model_type,' model file (atmospheric_model.dat)? > '],'s');
    % Extract atmospheric model.
    [Z, ~, ~, C, ~, ~, ~, ~, ~, ~, ~, ~, ~, W, ~, ~, ~] = extract_atmos_model(atmospheric_model_file, 3, 0, 0);
    % Compute atmospheric model Mach.
    Mach = W./C;
    % Deduce some atmospheric (DX, DT).
    dx_c = C./(f0*np);
    dx_c_mach = C.*(1-Mach)/(f0*np);
    dt_c_cfl = CFL*dx_c_mach*percentGLL./C;
    dt_c_mach = (1-Mach)./(np*f0*(1+Mach));
  end
  if(solidmodel)
    [interfaces_solid, nelts_solid, IDregionSolid] = model_viscoelastic_getinterfaces(f0, np);
    Nlayerz_solid=numel(nelts_solid);
  end
  
  % Figures.
  if(fluidmodel)
    if(plot_fluidmodel_dz)
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
    end

    % Displays.
    disp(['[',mfilename,'] Maximum acceptable DX with this atmospheric model is ',sprintf('%.3e',min(dx_c)),' [m] (c).']);
    disp([blanks(length(mfilename)+2),'    -        -      -   -    -        -       -    -  ',sprintf('%.3e',min(dx_c_mach)),' [m] (c, & Mach).']);
    disp(['[',mfilename,'] Maximum acceptable DT with this atmospheric model is ',sprintf('%.3e',min(dt_c_cfl)),' [s] (Mach-friendly DX, CFL ',sprintf('%.3e',CFL),', & ',sprintf('%.2g',percentGLL*100),'%GLL).']);
    disp([blanks(length(mfilename)+2),'    -        -      -   -    -        -       -    -  ',sprintf('%.3e',min(dt_c_mach)),' [s] (c, & Mach).']);
  end
  
  % Ask user.
  disp(' ');
  generatemeshfemparams = [];
  while (not(length(generatemeshfemparams) == 1 && ismember(generatemeshfemparams,[0,1])))
    generatemeshfemparams = input(['[',mfilename,'] Interactively generate meshfem parameters (0 for no, 1 for yes)? > ']);
  end
  
  if(generatemeshfemparams)
    if(fluidmodel)
      % Compute fluid layers.
      dz=dx_c_mach; % choice for atmospheric model
      zmax = [];
      while (not(length(zmax) == 1 && zmax<=max(Z) && zmax>min(Z)))
        zmax = input(['[',mfilename,'] Maximum atmospheric altitude of interest (',num2str(min(Z)),' < zmax < ',num2str(max(Z)),') [m]? > ']);
      end
      maxlayerz=ceil(length(Z(Z<=zmax))/5);
      if(maxlayerz==1)
        Nlayerz_fluid = 1;
      else
        Nlayerz_fluid = [];
        while (not(length(Nlayerz_fluid) == 1 && Nlayerz_fluid<=maxlayerz && Nlayerz_fluid>=1))
          Nlayerz_fluid = input(['[',mfilename,'] Number of atmospheric layers (1 <= N <= ',num2str(maxlayerz),') [1]? > ']);
        end
      end
      %interfaces_fluid=linspace(min(Z),zmax,Nlayerz_fluid+1);
      interfaces_fluid=linspace(min(Z),zmax+min(Z),Nlayerz_fluid+1) - min(Z); % Shift wanted model, in order to glue bottom to z=0 (convention).
    end
    
    % Safeguard checking: first air interface MUST BE exactly the same as last solid interface (if any).
    if(fluidmodel && solidmodel)
      if(not(interfaces_fluid(1)==interfaces_solid(end)))
        error(['[',mfilename,', ERROR] First air interface (z=',num2str(interfaces_fluid(1)),') is not exactly the same as last solid interface (z=',num2str(interfaces_solid(end)),').']);
      end
    end
    
    if(fluidmodel)
      % Compute elements per layer.
      for l=1:Nlayerz_fluid
        bestdx_fluid(l)=min(dz(Z>interfaces_fluid(l) & Z<interfaces_fluid(l+1)));
        nelts_fluid(l,1)=ceil((interfaces_fluid(l+1)-interfaces_fluid(l))/bestdx_fluid(l));
      end
    end
    if(solidmodel)
      dz_solid=diff(interfaces_solid)./nelts_solid;
    end
    
    if(fluidmodel && solidmodel)
      overall_minmaxdz = [min(min(bestdx_fluid),min(dz_solid)), max(max(bestdx_fluid),max(dz_solid))];
    elseif(not(fluidmodel) && solidmodel)
      overall_minmaxdz = [min(dz_solid), max(dz_solid)];
    elseif(fluidmodel && not(solidmodel))
      overall_minmaxdz = [min(bestdx_fluid), max(bestdx_fluid)];
    else
      error('ERROR THAT SHOULD NOT BE REACHED');
    end
    
    
    % Ask user for parfile nx (print only, has no impact on dz determination).
    nx = [];
    while (not(length(nx) == 1 && nx>0))
      disp(['[',mfilename,'] With these layers, ',num2str(overall_minmaxdz(1)),' <= dz <= ',num2str(overall_minmaxdz(2)),'.']);
      nx = input(['[',mfilename,'] nx (as in parfile) [1]? This has nothing to do with parameters computation, it''s only for printing. > ']);
    end
    
    % Prepare formats.
    if(fluidmodel && solidmodel)
      Nlayerz_total=Nlayerz_fluid+Nlayerz_solid;
      nelts_total=[nelts_solid;nelts_fluid];
    elseif(not(fluidmodel) && solidmodel)
      Nlayerz_total=Nlayerz_solid;
      nelts_total=nelts_solid;
    elseif(fluidmodel && not(solidmodel))
      Nlayerz_total=Nlayerz_fluid;
      nelts_total=nelts_fluid;
    else
      error('ERROR THAT SHOULD NOT BE REACHED');
    end
    format_layerz=['%',num2str(floor(log10(Nlayerz_total))+1),'i'];
    format_elts=['%',num2str(floor(log10(sum(nelts_total)))+1),'i'];
    format_dx='%10.5g';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interfaces for              %
    % interfacefile.              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(' ');
    disp('# Number of Interfaces:');
    disp(num2str(Nlayerz_total+1));
    if(solidmodel)
      % Print solid interfaces.
      for i=1:Nlayerz_solid % Do not print last interface, as it is the same as the first air one.
        if i==1; disp(['# Solid interface n°',num2str(i),' (bottom of mesh):']);
        elseif i==(Nlayerz_solid+1); disp(['# Solid interface n°',num2str(i),' (ground):']);
        else; disp(['# Solid interface n°',num2str(i),':']);
        end
        disp('2');
        disp(['-1d10 ',sprintf('%.15e',interfaces_solid(i))]);
        disp([' 1d10 ',sprintf('%.15e',interfaces_solid(i))]);
      end
    end
    if(fluidmodel)
      % Print fluid interfaces.
      for i=1:Nlayerz_fluid+1
        if i==1; disp(['# Air interface n°',num2str(i),' (ground):']);
        elseif i==(Nlayerz_fluid+1); disp(['# Air interface n°',num2str(i),' (top of the mesh):']);
        else; disp(['# Air interface n°',num2str(i),':']);
        end
        disp('2');
        disp(['-1d10 ',sprintf('%.15e',interfaces_fluid(i))]);
        disp([' 1d10 ',sprintf('%.15e',interfaces_fluid(i))]);
      end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Layers for interfacefile.   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    txtlayer_fluid={};
    disp('# Number of elements per layer:');
    if(solidmodel)
      % Print solid layers.
      for l=1:Nlayerz_solid
        txtlayer_solid{l}=['Solid layer n°',sprintf(format_layerz,l),' (from z = ',sprintf('%6.5g',interfaces_solid(l)),' to z = ',sprintf('%6.5g',interfaces_solid(l+1)),', dz = ',sprintf(format_dx,(interfaces_solid(l+1)-interfaces_solid(l))/nelts_solid(l)),').'];
        disp([sprintf(format_elts,nelts_solid(l)), ' # ',txtlayer_solid{l}]);
      end
    end
    if(fluidmodel)
      % Print fluid layers.
      for l=1:Nlayerz_fluid
        txtlayer_fluid{l}=['Air   layer n°',sprintf(format_layerz,l),' (from z = ',sprintf('%6.5g',interfaces_fluid(l)),' to z = ',sprintf('%6.5g',interfaces_fluid(l+1)),', dz = ',sprintf(format_dx,(interfaces_fluid(l+1)-interfaces_fluid(l))/nelts_fluid(l)),').'];
        disp([sprintf(format_elts,nelts_fluid(l)), ' # ',txtlayer_fluid{l}]);
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Regions for parfile.        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(' ');
    previouslyregisteredelements=0;
    disp(['nbregions = ',num2str(Nlayerz_total)])
    if(solidmodel)
      % Print solid layers.
      for l=1:Nlayerz_solid
        disp(['1 ',num2str(nx), ' ', sprintf(format_elts,previouslyregisteredelements+1), ' ', sprintf(format_elts,previouslyregisteredelements+nelts_solid(l)), ' ',num2str(IDregionSolid(l)), ' # ',txtlayer_solid{l}]);
        previouslyregisteredelements=previouslyregisteredelements+nelts_solid(l);
      end
    end
    if(fluidmodel)
      % Print fluid layers.
      for l=1:Nlayerz_fluid
        disp(['1 ',num2str(nx), ' ', sprintf(format_elts,previouslyregisteredelements+1), ' ', sprintf(format_elts,previouslyregisteredelements+nelts_fluid(l)), ' ',num2str(IDregionDG), ' # ',txtlayer_fluid{l}]);
        previouslyregisteredelements=previouslyregisteredelements+nelts_fluid(l);
      end
    end
    
    disp(' ');
    disp(['[',mfilename,'] Copy-paste the two previous blocks, respectively in the interfacefile and the parfile.']);
    
    if(fluidmodel)
      disp(' ');
      disp(['[',mfilename,'] For fluid, dt_max = CFL * min(dx) * pGLL / max(c). Here, safest is to choose dt = ',sprintf('%.3e',0.49*min(bestdx_fluid)*percentGLL/max(C)),'.']);
    end
  end
end