function [] = model_getdxdt(atmospheric_model_file)

  format compact;
  set(0, 'DefaultLineLineWidth', 3); set(0, 'DefaultLineMarkerSize', 8);
  set(0, 'defaultTextFontSize', 14); set(0, 'defaultAxesFontSize', 14);
  set(0, 'DefaultTextInterpreter', 'latex');
  set(0, 'DefaultLegendInterpreter', 'latex');
  
  [Z, ~, ~, C, ~, ~, ...
   ~, ~, ~, ~, ~, ~, ~, W, ~, ~, ~] = ...
   extract_atmos_model(atmospheric_model_file, 3, 0, 0);
 
  f0 = [];
  while (not(length(f0) == 1))
    f0 = input(['[',mfilename,'] Main temporal frequency [Hz]? > ']);
  end
  f0txt=['$f_0=',sprintf('%.3e',f0),'$'];
  
  np=2;
  CFL=0.49;
  percentGLL=0.18;
  
  Mach=W./C;
  
  dx_c=C./(f0*np);
  dx_c_mach=C.*(1-Mach)/(f0*np);
  
  figure();
  plot(dx_c_mach,Z);
  xlabel('$\Delta x$ [m]');
  ylabel('altitude $z$ [m]');
  title({['maximum acceptable $\Delta x$'],['accounting for Mach'],['(',f0txt,')']});
  set(gca,'ticklabelinterpreter','latex');
  grid on;
  
  dt_c_cfl=CFL*dx_c_mach*percentGLL./C;
  dt_c_mach=(1-Mach)./(np*f0*(1+Mach));
  
%   figure();
%   plot(dtaccountingformach,Z);
%   xlabel('maximum acceptable $\Delta t$, accounting for Mach number [s]');
%   ylabel('altitude $z$ [m]');
%   title(['$f_0=',sprintf('%.3e',f0),'$']);
%   set(gca,'ticklabelinterpreter','latex');
  
  disp(['[',mfilename,'] Maximum acceptable DX with this atmospheric model is ',sprintf('%.3e',min(dx_c)),' [m] (c).']);
  disp([blanks(length(mfilename)+2),'    -        -      -   -    -        -       -    -  ',sprintf('%.3e',min(dx_c_mach)),' [m] (c, & Mach).']);
  disp(['[',mfilename,'] Maximum acceptable DT with this atmospheric model is ',sprintf('%.3e',min(dt_c_cfl)),' [s] (Mach-friendly DX, CFL ',sprintf('%.3e',CFL),', & ',sprintf('%.3e',percentGLL*100),'%GLL).']);
  disp([blanks(length(mfilename)+2),'    -        -      -   -    -        -       -    -  ',sprintf('%.3e',min(dt_c_mach)),' [s] (c, & Mach).']);
  
end

