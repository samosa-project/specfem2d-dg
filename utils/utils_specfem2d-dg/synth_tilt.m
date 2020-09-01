% Author:        Léo Martire.
% Description:   TODO.
% Notes:         synth_load.m should have been ran before.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

% clear all;
% close all;
% clc;

[SPCFMEXloc] = setup_overall();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Artificial data for testing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nt=100;
% t_true=linspace(0,10,Nt);
% f_true=3/10;
% tilt_true= 0.2*(sin(2*pi*f_true*t_true) + 0.5*sin(2*pi*f_true*5*t_true));
% figure();
% subplot(321); plot(t_true,tilt_true*180/pi); subplot(322); plot(t_true,tilt_true*180/pi);
% delta=5;
% d_g_true=-delta*tan(tilt_true);
% d_d_true=+delta*tan(tilt_true);
% % figure();
% subplot(323);plot(t_true,d_g_true);subplot(324);plot(t_true,d_d_true);
% v_g_true=gradient(d_g_true,t_true);v_d_true=gradient(d_d_true,t_true);
% % figure();
% subplot(325);plot(t_true,v_g_true);subplot(326);plot(t_true,v_d_true);
% Ztime=[t_true;t_true;t_true;t_true];
% noisefac=0.25;
% noise_g=noisefac*(max(v_g_true)-min(v_g_true))*rand(size(v_g_true));
% noise_d=noisefac*(max(v_d_true)-min(v_d_true))*rand(size(v_d_true));
% Zamp=[v_g_true;v_d_true;v_g_true+noise_g;v_d_true+noise_d];
% xstattab=[-delta,delta,10-delta,10+delta];
% synth_load_was_ran=1;
% ID_l=1:2:3;
% ID_r=2:2:4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (not(exist('synth_load_was_ran', 'var') && synth_load_was_ran == 1))
  error(['[', mfilename, ', ERROR] synth_load was not ran before.']);
end

times=Ztime;
values=Zamp;
abscissas=xstattab(istattab);

% Cut times/values array based on shortest relevant array.
difftimes = times(:,2:end)-times(:,1:end-1);
relevantdifftimes= (difftimes>0);
minimum_last_relevant=+Inf;
for i=1:size(times,1)
  locminim=find(relevantdifftimes(i,:)>0,1,'last');
  if(locminim<minimum_last_relevant)
    minimum_last_relevant=locminim;
  end
end
minimum_last_relevant
times=times(:,1:minimum_last_relevant);
values=values(:,1:minimum_last_relevant);

ID_l = - 1;
while (not(min(ID_l)>=1 & max(ID_l)<=size(times,1)))
  ID_l = input(['[', mfilename, '] Matlab ID for left stations? > ']);
end
ID_r = - 1;
while (not(min(ID_r)>=1 & max(ID_r)<=size(times,1)))
  ID_r = input(['[', mfilename, '] Matlab ID for right stations? > ']);
end

% times_l = times(ID_l,:);
% times_r = times(ID_r,:);
% Nt = size(times_l,2);
% vals_l = values(ID_l,:);
% vals_r = values(ID_r,:);
% x_l = abscissas(ID_l);
% x_r = abscissas(ID_r);
% tilt = computeTilt(x_l, x_r, times_l, times_r, vals_l, vals_r);
tilt = computeTilt(abscissas(ID_l), abscissas(ID_r), times(ID_l,:), times(ID_r,:), values(ID_l,:), values(ID_r,:));

disp(['[',mfilename,'] Tilt computed for all stations.']);

disp(['[',mfilename,'] Plot one by one with:']);
disp([blanks(length(mfilename)+2),'   i = 1; figure(); plot(times_l(i, :), tilt(i, :));']);
disp([blanks(length(mfilename)+2),' With i in [1, ',num2str(N),'].']);

dist = 0.5*(x_l+x_r);
tilt2plot=tilt*180/pi;
timssss=times_l;

plot_time_v_dist(timssss, tilt2plot, dist, 'Tilt [$^\circ$]', 'horizontal position $x$');