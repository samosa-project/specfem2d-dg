function [VISCOUS] = MMS_check(testCase, MU, KAPPA)
  [RHO_cst, VX_cst, VZ_cst, E_cst, dRHO_x, dRHO_z, dVX_x, dVX_z, dVZ_x, dVZ_z, dE_x, dE_z] = MMS_constants(testCase);
  MU=0.35; %test
  if(max(MU,KAPPA)>0)
    VISCOUS = 1;
%     disp(['[] VISCOUS']);
    if(strcmp(testCase,'inviscid'))
      error(['[, ERROR] Cannot test inviscid with viscosity activated in this OUTPUT_FILES.'])
    end
    if(RHO_cst~=0)
      error(['[ERROR] for viscous tests, RHO_cst must be zero, and now is not']);
    end
    if(strcmp(testCase,'kappa') && not(VX_cst==0 & VZ_cst==0))
      error(['[ERROR] for viscous kappa test, VX_cst and VZ_cst must be zero, and now are not']);
    end
    if(strcmp(testCase,'kappa') && E_cst==0)
      error(['[ERROR] for viscous kappa test, E_cst must be non zero, and now is zero']);
    end
    if(strcmp(testCase,'mu') && not(E_cst==0))
      error(['[ERROR] for viscous mu test, E_cst must be zero, and now is not']);
    end
    if(strcmp(testCase,'mu') && VX_cst==0)
      error(['[ERROR] for viscous mu test, VX_cst must be non zero, and now is zero']);
    end
  else
    VISCOUS = 0;
%     disp(['[] INVISCID']);
    if(strcmp(testCase,'kappa') || strcmp(testCase,'mu'))
      error(['[, ERROR] Cannot test viscous with viscosity deactivated in this OUTPUT_FILES.'])
    end
    if(not(VX_cst==0 & VZ_cst==0))
%       error(['[ERROR] for inviscid test, VX_cst and VZ_cst must be zero, and now are not']);
    end
    if(RHO_cst==0)
      error(['[ERROR] for inviscid test, RHO_cst must be non zero, and now is zero']);
    end
  end
end