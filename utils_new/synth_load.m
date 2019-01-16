% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clc;
% clear all;
clear('Zamp','Ztime'); disp(['[',mfilename,'] Cleared Zamp and Ztime variables.']);
% close all;
format compact;
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

addpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/Atmospheric_Models');
SPCFMloc = '/home/l.martire/Documents/SPECFEM/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters' values. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Those can be re-set inline for each OUTPUT_FILES directory (see "OUTPUT_FILES location" section below).
rescale_factor = 1; % Rescaling: by default, do no rescale.
coord_units = 'km'; % Self-explanatory.
convert_to_relative_coords = 0; pos_interface = 0; % Convert to relative coordinates: x=0 above source, z=0 on surface (defined by pos_interface).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot amplitude?
%   0 for no.
%   1 for yes, sorted by x.
%   2 for yes, sorted by z.
%   3 for yes, sorted by d.
plot_amplitude = 0;

subsample = 0; % Sub-sample? Useful for lengthy seismograms. If set to 1, sub-sample so that synthetics are nsublength (below) long.
nsublength = 1000; % Length of sub-sampled synthetics.

% Quantity to display (should be the same as the seismotype variable in parfile):
%   1 = displacement for non-DG and velocity for DG,
%   2 = velocity for non-DG and pressure perturbation (Pa) for DG.
% type_display = 1;
type_display = 2;

% Unknown:
% unknown = 'BXX'; % _x.
unknown = 'BXZ'; % _z.
% For type_display==2 and stations in DG zones, pressure perturbation (Pa) is saved both in BXX and BXZ files.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT_FILES location.       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TNTGlanes
fig_title = strcat('Validation LNS');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0_FNS/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_FNS/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_0.5cfl/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_FNS_0.5cfl/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.15/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.15_FNS/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0_refined/')
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_refined/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0_gmsh/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_gmsh/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0_gmsh_refined_121954/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_gmsh_refined_1217753/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_corrected/'); % FIRST RUN AFTER PATCH
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_corrected/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_gmsh/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_gmsh/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_gmshrefined_1218660/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_gmshrefined_1218665/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0_FNS_testlambda0.2/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_FNS_testlambda0.2/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0_FNS_testlambda1.2/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_FNS_testlambda1.2/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0_testlambda0.2/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_testlambda0.2/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0_testlambda1.2/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_testlambda1.2/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0_chglambda/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_chglambda/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0_chglambda2woperturb/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M=0.3_chglambda2woperturb/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_1222576_M0_refined/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_1222565_M.3_refined/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M0_dx1_cfl.245/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.3_dx1_cfl.245/');
rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/validation_lns/'); OFd = strcat(rootd, 'OUTPUT_FILES_M.6_dx1_cfl.245/');

% TNTGlanes
% fig_title = strcat('Tirs de Mine Glanes');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/tirdemine_75040_redone/'); OFd = strcat(rootd, 'OUTPUT_FILES_1240682_nospurious_butbadgeom/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/tir_de_mine/'); OFd = strcat(rootd, 'OUTPUT_FILES_1216334/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/tir_de_mine/'); OFd = strcat(rootd, 'OUTPUT_FILES_1204916_spuriousreflexions/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/tir_de_mine/'); OFd = strcat(rootd, 'OUTPUT_FILES_75040/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/tir_de_mine/'); OFd = strcat(rootd, 'OUTPUT_FILES_74752/');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/tntglanes_10/'); OFd = strcat(rootd, 'OUTPUT_FILES_full/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/tntglanes_10/'); OFd = strcat(rootd, 'OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/tntglanes_10/'); OFd = strcat(rootd, 'OUTPUT_FILES_long600hz/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/tntglanes_10/'); OFd = strcat(rootd, 'OUTPUT_FILES_long300hz/');

% Microbaroms ULDB.
% fig_title = strcat('Microbaroms, (49N, 178W), 6:00 UT');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1206217/'); % Basically same as 1205575.
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1205575/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1204148_LNS/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_1203633_FNS/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_74710/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/mb_gmsh/'); OFd = strcat(rootd, 'OUTPUT_FILES_74565/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/mb_huge/'); OFd = strcat(rootd, 'OUTPUT_FILES_672048/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/mb_huge/'); OFd = strcat(rootd, 'OUTPUT_FILES_642746/');

% StratoBaro, 66, June, 12:00
% fig_title = strcat('Microbaroms, lat66, June, 12:00');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/microbaroms_patch'); OFd = strcat(rootd, '/OUTPUT_FILES_668482_disp7_isp6_full/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/microbaroms_patch'); OFd = strcat(rootd, '/OUTPUT_FILES_668482_disp7/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/microbaroms_patch'); OFd = strcat(rootd, '/OUTPUT_FILES_668446_disp7_wrongstations/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_668354_testlarger_str_1e-1mps_isp6/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_656744_straight_1mps_isp6/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_656505_straight_1e-2mps_isp6/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_656465_straight_1mps_test/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_655513_analytic/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_655494_alright/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_655487/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/microbaroms_periodic'); OFd = strcat(rootd, '/OUTPUT_FILES_655369_unstable/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_650851/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_644923_EBF_ispread3/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_641616_EBF_ispread1.5/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_641395_testEBF_stopped/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_639014_long_betterEBF_crash/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/stratobaro_test_EBF/'); OFd = strcat(rootd, '/OUTPUT_FILES/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/stratobaro_test_EBF/'); OFd = strcat(rootd, '/OUTPUT_FILES_test1/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_637450_long_EBF_crash/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_634307_testexternalforcing/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_624650_long/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_624515_rpw_spatially_fixed_s0.2/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_624478_apo+rpw0.2_10.5p/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_624436_apo+rpw0.15/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_623945_apodised/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_622147_final_test/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_622125/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_622037/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_621952/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_621860/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_621802/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratobaro_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_599638_testmicrobarom/');

% StratoExplo, 66, June, 12:00
% fig_title = strcat('Atmospheric Explosions, lat66, June, 12:00');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_597316/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_597250/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_597099/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_595500_crash40k1it/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRdfsdfATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_594736_crash27kit/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_594361_dt1e-3_cancelled/');

% Seismic Hammer, soft soil.
% fig_title = strcat('Seismic Hammer Simulation (Soft Soil)'); coord_units = 'm'; convert_to_relative_coords = 1; pos_interface = 308;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final'); OFd = strcat(rootd, '/OUTPUT_FILES_669168_fullretweaked/'); rescale_factor = 236; % Same as 593959 but test with first layers changed.
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final'); OFd = strcat(rootd, '/OUTPUT_FILES_668888_stopped_12kit/'); rescale_factor = 236; % Same as 593959 but test with first layers changed.
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_axisym'); OFd = strcat(rootd, '/OUTPUT_FILES_660223_full_dec1m/'); % Same as 593959 but axisymmetric.
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final'); OFd = strcat(rootd, '/OUTPUT_FILES_627577_qk4sls_truefreesurf/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final'); OFd = strcat(rootd, '/OUTPUT_FILES_623195_qk_4sls_freesurf/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final/SH_soft_final_redone'); OFd = strcat(rootd, '/OUTPUT_FILES_610770/'); % With additionnal stations for comparison with data.
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final'); OFd = strcat(rootd, '/SH_soft_final_redone_Qkappa_616368/'); % Same as 593959 only with Qp converted to Qk and additionnal stations.
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final'); OFd = strcat(rootd, '/SH_soft_final_redone_Qkappa+f=0_618645/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final'); OFd = strcat(rootd, '/SH_soft_final_redone_Qkappa+f=f0_618882/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final'); OFd = strcat(rootd, '/SH_soft_final_redone_qk_noatt_619264');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final'); OFd = strcat(rootd, '/SH_soft_final_redone_qk_att4sls_620294');
% rootd = strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final/SH_soft_final'); OFd = strcat(rootd, '/OUTPUT_FILES_593959/'); % Original (used in paper).
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_new_model'); OFd = strcat(rootd, '/OUTPUT_FILES_551980_seismic_potential_with_memvars_solid/');

% Seismic Hammer, hard soil.
% fig_title = strcat('Seismic Hammer Simulation (Hard Soil)'); coord_units='m'; convert_to_relative_coords = 1; pos_interface=308;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final/SH_hard_final'); OFd = strcat(rootd, '/OUTPUT_FILES_593960/'); % Original (used in paper).
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_hard_axisym'); OFd = strcat(rootd, '/OUTPUT_FILES_661601_full_dec1m/'); % Same as 593960 but axisymmetric.
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_hard_axisym'); OFd = strcat(rootd, '/OUTPUT_FILES_661609_full_onlypress/'); type_display=4; unknown = 'PRE'; % Same as 661601 but only recording above ground.

% Quake, 45.
% fig_title = strcat('Quake Simulation (45$^\circ$ dip)');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_1811221612_local');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_1811221544_local');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_impedance'); OFd = strcat(rootd, '/OUTPUT_FILES');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_75118_isoth_d6_savedvdg');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_71984_isoth_d6');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_71980_isothermal');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_pot');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_71931_d9');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_71920_force_instead_of_moment_d9');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ_test_imp'); OFd = strcat(rootd, '/OUTPUT_FILES_71913');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ'); OFd = strcat(rootd, '/OUTPUT_FILES_668844_OKQ45_redone'); rescale_factor = 1e-3;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ/ON_EOS_quake_ok_45'); OFd = strcat(rootd, '/OUTPUT_FILES_583041_long');

% Quake, 0.
% fig_title = strcat('Quake Simulation (0$^\circ$ dip)');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ'); OFd = strcat(rootd, '/OUTPUT_FILES_668833_OKQ0_redone'); rescale_factor = 1e-3;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/OKQ/ON_EOS_quake_ok_0'); OFd = strcat(rootd, '/OUTPUT_FILES_586984_full');

% Tests.
% fig_title = 'test';
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/tir_mars'); OFd = strcat(rootd, '/OUTPUT_FILES/');

% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/demo_pot'); OFd = strcat(rootd, '/OUTPUT_FILES_826234/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/demo_fns'); OFd = strcat(rootd, '/OUTPUT_FILES_826226/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/demo_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_826213/');

% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/demo_fns'); OFd = strcat(rootd, '/OUTPUT_FILES_fnsf2s_local/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/demo_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_lnsf2s_local/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/demo_fns'); OFd = strcat(rootd, '/OUTPUT_FILES_fnsf2s_local_butd7/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/demo_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_lnsf2s_local_butd7/');

% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/demo_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_lns_t=105s/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/demo_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_fns_t=195s/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/mb_gmsh'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_impedance'); OFd = strcat(rootd, '/OUTPUT_FILES_extatm+oksoil+lowdt_BUTONCALMIP/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_impedance'); OFd = strcat(rootd, '/OUTPUT_FILES_extatm+oksoil+lowdt/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_impedance'); OFd = strcat(rootd, '/OUTPUT_FILES_lns/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_fts_lns_19s/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_fts_fns_26s/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_f_lns_44s/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_f_fns_62s/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_stf_lns_27s/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_lns'); OFd = strcat(rootd, '/OUTPUT_FILES_stf_fns_45s/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/axisym_test'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_external_forcing'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_stretching_wind'); OFd = strcat(rootd, '/OUTPUT_FILES/'); coord_units='m'; convert_to_relative_coords = 0; pos_interface=0;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_stretching_wind'); OFd = strcat(rootd, '/OUTPUT_FILES_observesignalinbuffer_cstrhdr/'); coord_units='m'; convert_to_relative_coords = 0; pos_interface=0;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf4/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf10_1hz/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf2_1hz/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf5_jpeguz/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf31hz_homo_otherstation/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf31hz_homo/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf4_homogenous/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_sharpstf3/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf5/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf3/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf2/'); type_display = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stations' loading.           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test if OUTPUT_FILES directory exists.
if (not(strcmp(OFd(end), '/')))
  OFd = [OFd, '/'];
end
if (not(exist(OFd, 'dir')))
  error(['[',mfilename,', ERROR] OUTPUT_FILES directory does not exist (', OFd, ').']);
end

% Load sources' positions.
pos_sources = [inf, inf]; % Allocate a row for the first source's position.
% fid = fopen([rootd, '/DATA/SOURCE']);
fid = fopen([OFd, 'SOURCE']);
if (fid == - 1)
  fid = fopen([OFd, 'input_source']);
  if (fid == - 1)
    error(['[',mfilename,', ERROR] Cannot open SOURCE file (', OFd, 'SOURCE', ').']);
  end
end
line = 0; xfound = 0; zfound = 0; stop = 0;
while (stop == 0)
  % TODO: Loop on source number.
  line = fgetl(fid);
  if length(line) > 0
    if (line == - 1)
      stop = 1;
    end
    line = regexprep(regexprep(line, ' +', ' '), '^ ', ''); % Remove multiple spaces, and then eventually remove space if it there is one as first character.
    if strcmp(line(1:2), 'xs')
      xfound = 1; pos_sources(1, 1) = str2num(regexprep(regexprep(line(3:end), ' *#.*', ''), ' *=* *', '')); % Remove comments (everything after a '#'), remove the equals sign and spaces around it, and cast it as source position.
    end
    if strcmp(line(1:2), 'zs')
      zfound = 1; pos_sources(1, 2) = str2num(regexprep(regexprep(line(3:end), ' *#.*', ''), ' *=* *', ''));
    end
  end
  if (xfound && zfound)
    break
  end
end
fclose('all');

% Load stations data (first try OUTPUT folder, then if not found, try parent DATA folder).
try
  A = importdata(strcat(OFd, 'STATIONS'));
catch
  disp(['[',mfilename,'] STATIONS file not found in OUTPUT_FILES directory.']);
  try
    A = importdata(strcat(rootd, '/DATA/STATIONS'));
    disp(['[',mfilename,'] STATIONS file found in root directory (OUTPUT_FILES*/../DATA/ folder).']);
  catch
    error(['[',mfilename,', ERROR] Cannot find STATIONS file.']);
  end
end
pos_stations = [A.data(:, 1) A.data(:, 2)]; xstattab = pos_stations(:, 1); ystattab = pos_stations(:, 2);
% Compute distance to sources.
dist_to_sources = zeros(size(pos_stations, 1), size(pos_sources, 1));
for n_source = 1:size(pos_sources, 1)
  dist_to_sources(:, n_source) = sqrt((pos_stations(:, 1) - pos_sources(n_source, 1)) .^ 2 + (pos_stations(:, 2) - pos_sources(n_source, 2)) .^ 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ask for user input.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display stations' information.
format shortG;
if (convert_to_relative_coords == 1)
  disp(['[',mfilename,'] [station_id; x; z; d] for all stations (x relative to source, z relative to ground, d relative to source):']); disp([(1:size(pos_stations, 1)); (pos_stations - [pos_sources(1, 1), pos_interface])';dist_to_sources']);
else
  disp(['[',mfilename,'] [station_id; x; z; d] for all stations (absolute x and z, d relative to source):']); disp([(1:size(pos_stations, 1)); pos_stations';dist_to_sources']);
end
format compact;

% Ask for behaviour.
display_or_load = - 1;
while (not(length(display_or_load) == 1 && ismember(display_or_load, [0, 1, 2])))
  display_or_load = input(['[',mfilename,'] Load and display (0), load only (1), or load and combine (2)? > ']);
end
% Ask for stations.
istattab = input(['[',mfilename,'] Stations (Matlab format, e.g. [1, 4, 7] or 1:20)? > ']);
istattab = reshape(istattab,[1,numel(istattab)]);
disp(['[',mfilename,'] Loading [station_id, x, z, d] (absolute x and z, d relative to source):']);
disp([istattab', pos_stations(istattab,:),dist_to_sources(istattab)]);
nstat = size(pos_stations(istattab, 1), 1);
% Ask for geometric attenuation (relies on distance to source).
geometric_attenuation = - 1;
while (not(ismember(geometric_attenuation, [0, 1, 2, 3])))
  geometric_attenuation = input(['[',mfilename,'] Apply geometric attenuation factor to data? (0 for no, 1 for d, 2 for |x|, 3 for |z|) > ']);
end
% Ask if plot y-axis should be normalised to same scale.
normalise_ylims = 0; % Default value.
if (display_or_load == 0 && nstat > 1)
  normalise_ylims = - 1;
  while (not(normalise_ylims == 0 || normalise_ylims == 1))
    normalise_ylims = input(['[',mfilename,'] Normalise y-scale? (0 for no, 1 for yes) > ']);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and eventually plot.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (display_or_load == 0)
  figure(); hold on;
end

% Loop on synthetics.
if (normalise_ylims)
  % Prepare updates of y-axis scale.
  max_ylim_plus = - Inf;
  min_ylim_minus = + Inf;
end
if (convert_to_relative_coords == 1)
  % Eventually remove source components for display.
  xstattab = xstattab - pos_sources(1, 1);
  ystattab = ystattab - pos_interface;
end
for istat = 1:nstat
  istat_glob = istattab(istat); % Recover global number of station.

  % Switch on type of display.
  if (type_display == 1)
    % Original SPECFEM2D's synthetic is displacement.
    extension = "semd"; % Because original SPECFEM2D's synthetic is displacement.
    % For stations in solid zones it's displacement. For stations in DG zones it's velocity.
    if (strcmp(unknown, 'BXZ'))
      unknown_name = 'vertical {$u_z$ (m), $v_z$ (m/s)}';
    elseif (strcmp(unknown, 'BXX'))
      unknown_name = 'horizontal {$u_x$ (m), $v_x$ (m/s)}';
    else
      error(['[',mfilename,', ERROR] The variable ''unknown'' has a non-standard value.']);
    end
  elseif (type_display == 2)
    % Original SPECFEM2D's synthetic is velocity.
    extension = "semv"; % Because original SPECFEM2D's synthetic is velocity.
    if (strcmp(unknown, 'BXZ'))
      unknown_name = '{$v_z$ (m/s), $\delta P$ (Pa)}';
    elseif (strcmp(unknown, 'BXX'))
      unknown_name = '{$v_x$ (m/s), $\delta P$ (Pa)}';
    else
      error(['[',mfilename,', ERROR] The variable ''unknown'' has a non-standard value.']);
    end
  elseif (type_display == 4)
    % Original SPECFEM2D's synthetic is pressure.
    extension = "semp"; % Because original SPECFEM2D's synthetic is pressure.
    if (strcmp(unknown, 'PRE'))
      unknown_name = '$\delta P$ (Pa)';
    else
      error(['[',mfilename,', ERROR] The variable ''unknown'' has a non-standard value.']);
    end
  end

  % Read the synthetic.
  file = strcat(OFd, 'AA.', A.textdata(istat_glob, 1), '.', unknown, '.', extension);
  data = load(file{1});
  nt = max(size(data));
  if (subsample == 1)
    % Sub-sample of records.
    nsub = ceil(nt / nsublength);
    nd = max(size(data(1:nsub:nt, 1)));
  else
    nsub = 1;
    nd = nt;
  end

  % Recover time/amplitude data.
  Ztime(istat, 1:nd) = data(1:nsub:nt, 1)';
  % Ztime(istat, :) = Ztime(istat, :) - Ztime(istat, 1); % Make time values start at zero.
  Zamp(istat, 1:nd) = data(1:nsub:nt, 2)';

  % Renormalisation (global, and geometric).
  factor = 1; % Reset to default value for each station.
  if (geometric_attenuation ~= 0)
    % If geometric_attenuation was asked by user.
    switch geometric_attenuation
      case 1
        geom_att_fact = dist_to_sources(istat_glob) ^ 0.5; % Geometric attenuation respective to raw distance to source.
      case 2
        geom_att_fact = abs(xstattab(istat_glob)) ^ 0.5; % Geometric attenuation respective to horizontal distance to source.
      case 3
        geom_att_fact = abs(ystattab(istat_glob)) ^ 0.5; % Geometric attenuation respective to vertical distance to source.
      otherwise
        error(['[',mfilename,', ERROR] Geometric attenuation parameter not implemented.']);
    end
    if (geom_att_fact == 0)
      % If exactly at zero distance, consider no geometric rescaling.
      geom_att_fact = 1;
    end
    factor = factor / geom_att_fact;
  end
  renorm_statbystat = - 1;
  if (rescale_factor ~= 1)
    % Rescaling was asked. Check again with user.
    renorm_statbystat = - 1;
    disp(strcat("    Specified rescale factor is ", num2str(rescale_factor), "."));
    inputtxt = char(strcat("    Rescale data for station ", num2str(istat_glob), "? (0 for no, 1 for yes) > "));
    while (not(ismember(renorm_statbystat, [0, 1])))
      renorm_statbystat = input(inputtxt);
      if (isempty(renorm_statbystat))
        renorm_statbystat = 1;
      end
    end
  end
  if (renorm_statbystat == 1)
    % If rescaling is actually wanted by user, do it.
    factor = factor * rescale_factor;
  end
  Zamp(istat, :) = factor * Zamp(istat, :); % Rescale.

  % Eventually, display.
  if (display_or_load == 0)
    ax(istat) = subplot(nstat, 1, istat);
    if (strcmp(coord_units, 'km'))
      legtext{istat} = strcat('S', num2str(istat_glob), ', $(x,z)=(', num2str(xstattab(istat_glob) / 1000), ',', num2str(ystattab(istat_glob) / 1000), '$) ', coord_units);
    elseif (strcmp(coord_units, 'm'))
      legtext{istat} = strcat('S', num2str(istat_glob), ', $(x,z)=(', num2str(xstattab(istat_glob)), ',', num2str(ystattab(istat_glob)), ')$ ', coord_units);
    else
      error(['coord_units = ', coord_units, 'not implemented.']);
    end
    plot(Ztime(istat, :), Zamp(istat, :));
    set(gca, 'TickLabelInterpreter', 'latex');
    % Cosmetics.
    if (istat == 1)
      title(fig_title)
    end
    if (istat == nstat)
      xlabel('time (s)')
    end
    if (istat ~= nstat)
      set(gca, 'xticklabel', []);
    end
    if (istat == round(nstat / 2))
      ylabel(unknown_name);
    end
    xlim([Ztime(1, 1), Ztime(1, end)]);
    legend(legtext{istat}, 'Location', 'northeast');
    hold on;
    if (normalise_ylims)
      % Update y-axis scale.
      if (ax(istat).YLim(1) < min_ylim_minus)
        min_ylim_minus = ax(istat).YLim(1);
      end
      if (ax(istat).YLim(2) > max_ylim_plus)
        max_ylim_plus = ax(istat).YLim(2);
      end
    end
    if (normalise_ylims)
      linkaxes(ax);
    else
      linkaxes(ax, 'x');
    end
  end
end

if (display_or_load == 0 && nstat > 1 && normalise_ylims == 1)
  % Normalise plot y-axis to same scale.
  f = gcf;
  for i = 1:length(f.Children)
    if (strcmp(f.Children(i).Type, 'axes'))
      f.Children(i).YLim = [min_ylim_minus, max_ylim_plus];
    end
  end
end

% Display information.
disp([' ']);
disp(['[',mfilename,'] Data loaded. [matlab_id, station_id, x, z, d]:']);
disp([(1:length(istattab))',istattab', pos_stations(istattab,:),dist_to_sources(istattab)]);
disp(strcat("  Example: Data of station ", num2str(istattab(1)), " are in         Zamp(", num2str(1), ", :)."));
disp(strcat("           Corresponding time values are in Ztime(", num2str(1), ", :)."));
disp([' ']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear('A', 'ans', 'ax', 'data', 'extension', 'f', 'fid', 'file', 'i', 'inputtxt', 'istat', 'istat_glob', 'line', 'max_ylim_plus', 'min_ylim_minus', 'nd', 'normalise_ylims', 'nt', 'nsub', 'nsublength', 'pos_stations', 'rootd', 'SPCFMloc', 'stop', 'subsample', 'unknown', 'xfound', 'zfound');
clear('A', 'ans', 'data', 'extension', 'f', 'fid', 'file', 'i', 'inputtxt', 'istat', 'istat_glob', 'line', 'max_ylim_plus', 'min_ylim_minus', 'nd', 'nt', 'nsub', 'nsublength', 'pos_stations', 'rootd', 'SPCFMloc', 'stop', 'subsample', 'unknown', 'xfound', 'zfound');
if (factor == 1)
  clear('factor');
end
if (geometric_attenuation == 0)
  clear('geometric_attenuation');
end
if (renorm_statbystat == 0)
  clear('renorm_statbystat', 'rescale_factor');
end

global synth_load_was_ran
synth_load_was_ran = 1;

if (display_or_load == 2)
  
  distancechoice = - 1;
  while (~ ismember(distancechoice, [1, 2, 3, 4]))
    distancechoice = input(['[', mfilename, '] Distance choice? (1 for x, 2 for |x|, 3 for z, 4 for d) > ']);
  end
  switch distancechoice
    case 1
      distance = xstattab; dist_symbol = "x"; dist_name = "horizontal distance";
    case 2
      distance = abs(xstattab); dist_symbol = "|x|"; dist_name = "horizontal distance";
    case 3
      distance = ystattab; dist_symbol = "z"; dist_name = "altitude";
    case 4
      distance = dist_to_sources; dist_symbol = "d"; dist_name = "distance";
  end
  addpath('/home/l.martire/Documents/Ongoing_Work/1811_glanes/treatment_leo');
  plot_time_v_dist(Ztime,Zamp,distance(istattab));
  
%   synth_plot;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old things.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (ismember(plot_amplitude, [1, 2, 3]))
  % Plot amplitude.
  for ii = 1:length(Zamp(:, 1))
    amp(ii) = max(Zamp(ii, :)) - min(Zamp(ii, :));
  end
  switch(plot_amplitude)
    case 1
      [sorted_dist, isort] = sort(xstattab(istattab));
    case 2
      [sorted_dist, isort] = sort(ystattab(istattab));
    case 3
      [sorted_dist, isort] = sort(dist_to_sources(istattab));
  end
  sorted_amp = amp(isort);
  figure();
  loglog(sorted_dist, sorted_amp);
end

% Seismic Hammer, hard soil.
% fig_title = strcat('Seismic Hammer Simulation (Hard Soil)'); coord_units='m'; convert_to_relative_coords = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final/SH_hard_final'); OFd = strcat(rootd, '/OUTPUT_FILES_593960/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580457_full/'); rescale_factor=8.840811261618920e-04;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580113/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580185/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580228/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580333/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS__seismic_hammer_hard_soil'); OFd = strcat(rootd, '/OUTPUT_FILES_580712/');

% Seismic Hammer, soft soil.
% fig_title = strcat('Seismic Hammer Simulation (Soft Soil)'); coord_units='m'; convert_to_relative_coords = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/SH_final/SH_soft_final'); OFd = strcat(rootd, '/OUTPUT_FILES_593959/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9113508_seismic_DG_with_memvars_solid/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9048100_seismic_DG/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9081476_seismic_potential/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9091088_seismic_DG_new_coupling/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_9102702_seismic_potential_rem_forcing/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Balloons/simulations'); OFd = strcat(rootd, '/OUTPUT_FILES_551980_seismic_potential_with_memvars_solid/');

% Tests.
% fig_title = 'test';
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf4/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf10_1hz/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf2_1hz/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf5_jpeguz/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf31hz_homo_otherstation/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf31hz_homo/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf4_homogenous/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_sharpstf3/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf5/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf3/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_stf2/'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_external_modelDG_only'); OFd = strcat(rootd, '/OUTPUT_FILES_1e0mu/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_external_modelDG_only'); OFd = strcat(rootd, '/OUTPUT_FILES_1e2mu/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_external_modelDG_only'); OFd = strcat(rootd, '/OUTPUT_FILES_1e4mu/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_external_modelDG_only'); OFd = strcat(rootd, '/OUTPUT_FILES_1e5mu/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/full_DG_square'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_stretching'); OFd = strcat(rootd, '/OUTPUT_FILES_long/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_FTS'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_coupling'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_stretching'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_stretching_wind'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_stretching_FFcounterpart'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_test_atmo'); OFd = strcat(rootd, '/OUTPUT_FILES_TEST');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_test_densitysource'); OFd = strcat(rootd, '/OUTPUT_FILES_TEST');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/test_quake'); OFd = strcat(rootd, '/OUTPUT_FILES');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES__long_working');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma0'); OFd = strcat(rootd, '/OUTPUT_FILES_test_vz');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_test_vz');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_oklahoma45'); OFd = strcat(rootd, '/OUTPUT_FILES_test_dz'); type_display = 1;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/quake_ok_45'); OFd = strcat(rootd, '/OUTPUT_FILES_narrow_okdx');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_583123_100km_3sources_nospread');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_583128_100km_3sources_nospread');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_583138_100km_3sources_spread');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_test_scale_sources');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200'); OFd = strcat(rootd, '/OUTPUT_FILES_591778_66j1200_regmukap0_softground_nocrash');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/seismic_hammer_zooms/soft/'); OFd = strcat(rootd, 'OUTPUT_FILES_583180');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/seismic_hammer_zooms/hard/'); OFd = strcat(rootd, 'OUTPUT_FILES_583194');
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/seismic_hammer_zooms/hard/'); OFd = strcat(rootd, 'OUTPUT_FILES_586795_d4_source_flipped'); rescale_factor=8.840811261618920e-04;
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/seismic_hammer_zooms/soft/'); OFd = strcat(rootd, 'OUTPUT_FILES_591776_additional_stations');fig_title = 'test_soft';
% rootd=strcat(SPCFMloc, 'specfem-dg-master/EXAMPLES/seismic_hammer_zooms/hard/'); OFd = strcat(rootd, 'OUTPUT_FILES_591777_additional_stations');fig_title = 'test_hard';

% Mars Gravity Wave.
% fig_title = strcat('Mars Gravity Wave Simulation');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_533937_vNEW_full/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_534758_long_instab/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_535011_with_FTS/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_535489_removed_discontinuity_long/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave/test_RAPHAEL'); OFd = strcat(rootd, '/OUTPUT_FILES/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_540064_FTS_no_disc_long/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9078210_spread_source/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9081352_spread_cut_source/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9091089_new_coupling/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_9103256_same_as_previous_but_factor_1/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_552471_atmo_only/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_552455_factor0p1/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_557219_long/');
% rootd=strcat(SPCFMloc, 'Ongoing_Work/Mars_Gravity_Wave'); OFd = strcat(rootd, '/OUTPUT_FILES_558183_Gderiv/');

% Mars AGW.
% rootd=strcat(SPCFMloc, 'Ongoing_Work/SPECFEM-DG_Mars_AGW_runs/explo_mars_sub'); OFd = strcat(rootd, '/OUTPUT_FILES_KappaON/');
