% Author:        Léo Martire.
% Description:   Computes ratio of pressure over vertical velocity, in
%                order to estimate impedance.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         /utils_new/synth_load.m has to have been ran
%                before.

% clear all;
% close all;
% clc;
format compact;

if (not(exist('synth_load_was_ran') && synth_load_was_ran == 1))
  error(['[', mfilename, ', ERROR] synth_load was not ran before.']);
end

idvz = - 1;
while (idvz == - 1)
  idvz = input(['[', mfilename, '] Matlab ID for v_z signal? > ']);
end
iddp = - 1;
while (iddp == - 1)
  iddp = input(['[', mfilename, '] Matlab ID for dp signal? > ']);
end

t=Ztime(1,:);
amp_vz=Zamp(idvz,:);
amp_dp=Zamp(iddp,:);
recordpointdp=['@(x,z) = (',sprintf('%.3g',xstattab(istattab(iddp))),', ',sprintf('%.3g',ystattab(istattab(iddp))),') [m]'];
recordpointvz=['@(x,z) = (',sprintf('%.3g',xstattab(istattab(idvz))),', ',sprintf('%.3g',ystattab(istattab(idvz))),') [m]'];

%   dt=mean(unique(diff(Ztime(1,:))));
dt=mean(diff(Ztime(1,:)));

smolest_period = - 1;
while (smolest_period == - 1)
  smolest_period = input(['[', mfilename, '] Smallest period to get [s]? > ']);
end

indices_to_get_smolest_period=floor(smolest_period/dt);

[vz_ue,vz_le]=envelope(amp_vz,indices_to_get_smolest_period,'analytic');
[dp_ue,dp_le]=envelope(amp_dp,indices_to_get_smolest_period,'analytic');

rel_threshold = 1e-2;
sel_id_vz = abs(amp_vz) > max(abs(amp_vz))*rel_threshold;
sel_id_dp = abs(amp_dp) > max(abs(amp_dp))*rel_threshold;
sel_id=logical(sel_id_vz.*sel_id_dp);
clear('sel_id_vz', 'sel_id_dp');

% ratio_u=dp_ue./vz_ue;
% ratio_l=dp_le./vz_le;
ratio_u=vz_ue./dp_ue;
ratio_l=vz_le./dp_le;


% figures
figure('units','normalized','outerposition',[0 0 0.5 1]);
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'DefaultTextInterpreter', 'latex'); set(0, 'DefaultLegendInterpreter', 'latex');
ax(1)=subplot(311);
plot(t, amp_dp, t, dp_le, ':', t, dp_ue, ':');
ylabel({['$\delta P$'],[recordpointdp]});
xlim([min(t), max(t)]); set(gca, 'TickLabelInterpreter', 'latex'); grid on;

ax(2)=subplot(312);
plot(t, amp_vz, t, vz_le, ':', t, vz_ue, ':');
ylabel({['$v_z$'],[recordpointvz]});
xlim([min(t), max(t)]); set(gca, 'TickLabelInterpreter', 'latex'); grid on;

ax(3)=subplot(313);
plot(t(sel_id), ratio_l(sel_id), t(sel_id), ratio_u(sel_id));
ylab=['$\left|v_z/\delta P\right|$'];
title({['Impedance ratio ',ylab,''],['(rel. signal threshold ',num2str(rel_threshold),', min. period ',num2str(smolest_period),' s)']});
ylabel(ylab);
xlabel('$t$ (s)');
xlim([min(t), max(t)]);% set(gca, 'TickLabelInterpreter', 'latex'); grid on;

linkaxes(ax,'x');

for i=1:3
  setaxformats(ax(i));
end
% xlim([12,18]);

% other displays
disp(['[', mfilename, '] Recall:     dp recorded ',recordpointdp,'.']);
disp([blanks(length(mfilename)+2),'         and vz          ',recordpointvz,'.']);

function setaxformats(ax)
  set(ax,'fontsize',18);
  set(ax, 'TickLabelInterpreter', 'latex');
  grid on;
  box on;
end