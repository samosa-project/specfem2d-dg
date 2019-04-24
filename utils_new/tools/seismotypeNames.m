% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   [variable_name, signal_letter] = seismotypeNames(seismotype, DG)
% with:
%   seismotype an integer (found in parfiles)
%   DG         a boolean encoding if the discontinuous method (SPECFEM-DG)
%              is being used.
% yields:
%   TODO.

function [variable_name, signal_letter] = seismotypeNames(seismotype, DG)
  if(not(exist('DG')))
    DG = 0;
  end
  switch(seismotype)
    case(1)
      if(DG)
        variable_name = 'Velocity';
      else
        variable_name = 'Displacement';
      end
      signal_letter = 'd' ; % -> 'a' for acceleration
      % -> 'd' for displacement or sqrt(density) *
      %    displacement
      % -> 'v' for velocity or sqrt(density) * velocity
      % Watch the right letter at the end of the output
      % ascii file from specfem2D
    case(2)
      if(DG)
        variable_name = 'Pressure';
      else
        variable_name = 'Velocity';
      end
      signal_letter = 'v' ;
    case(3)
      if(DG)
        variable_name = '$\sqrt{\rho}v_z$';
      else
        variable_name = 'Acceleration';
      end
      signal_letter = 'a' ;
    case(4)
      if(DG)
        disp(['[',mfilename,', INFO] seismotype ',num2str(seismotype),' not implemented for DG. Setting variable to the one from classic SPECFEM.']);
        variable_name = 'Pressure';
      else
        variable_name = 'Pressure';
      end
      signal_letter = 'p' ;
    case(5)
      if(DG)
        disp(['[',mfilename,', INFO] seismotype ',num2str(seismotype),' not implemented for DG. Setting variable to the one from classic SPECFEM.']);
        variable_name = 'Curl of Displacement';
      else
        variable_name = 'Curl of Displacement';
      end
      signal_letter = 'c' ;
    case(6)
      if(DG)
        disp(['[',mfilename,', INFO] seismotype ',num2str(seismotype),' not implemented for DG. Setting variable to the one from classic SPECFEM.']);
        variable_name = 'Fluid Potential';
      else
        variable_name = 'Fluid Potential';
      end
      signal_letter = 'p' ;                  
    %case(7)
    %  variable = '$\rho^{1/2}\cdot {\bf u}';
    %  signal_type = 'd' ;
    otherwise
      error(['[',mfilename,', ERROR] seismotype ',num2str(seismotype),' not implemented.']);
  end
end

