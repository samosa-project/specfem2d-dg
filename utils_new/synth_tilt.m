% Author:        Léo Martire.
% Description:   TODO.
% Notes:         /utils_new/synth_load.m has to have been ran
%                before.
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
format compact;

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
% idleft=1:2:3;
% idright=2:2:4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (not(exist('synth_load_was_ran') && synth_load_was_ran == 1))
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

idleft = - 1;
while (not(min(idleft)>=1 & max(idleft)<=size(times,1)))
  idleft = input(['[', mfilename, '] Matlab ID for left stations? > ']);
end
idright = - 1;
while (not(min(idright)>=1 & max(idright)<=size(times,1)))
  idright = input(['[', mfilename, '] Matlab ID for right stations? > ']);
end

times_l = times(idleft,:);
times_r = times(idright,:);
Nt = size(times_l,2);
vals_l = values(idleft,:);
vals_r = values(idright,:);
x_l = abscissas(idleft);
x_r = abscissas(idright);
N=size(vals_l,1);
Nd=size(vals_r,1);
if(N~=Nd)
  error('must provide same number of left stations as right stations');
end
clear('Nd');

displ_l = zeros([N, Nt]);
displ_r = zeros([N, Nt]);
tilt = zeros([N, Nt]);
for i=1:N
  displ_l(i,:) = cumtrapz(times_l(i, :), vals_l(i, :));
  displ_r(i,:) = cumtrapz(times_r(i, :), vals_r(i, :));
  tilt(i,:) = atan((displ_r(i, :)-displ_l(i, :))/(x_r(i)-x_l(i)));
  
%   trailZ_l=find(diff(times_l(i, :))==0); % find possible trailing zeros in the times array (happens if seismograms are not all exactly the same length
%   if(numel(trailZ_l)>0); trailZ_l=[trailZ_l(1)-1,trailZ_l]; end; % add first one (line above only finds using diff and so the very first one is not detected).
%   trailZ_r=find(diff(times_r(i, :))==0); if(numel(trailZ_r)>0); trailZ_r=[trailZ_r(1)-1,trailZ_r]; end;% idem with _r
%   trailZs=sort(unique([trailZ_l,trailZ_r])); % combine both arrays
%   tilt(i,trailZs)=NaN; % set tilt as nan for those
end

disp(['[',mfilename,'] Tilt computed for all stations.']);

disp(['[',mfilename,'] Plot one by one with:']);
disp([blanks(length(mfilename)+2),'   i = 1; figure(); plot(times_l(i, :), tilt(i, :));']);
disp([blanks(length(mfilename)+2),' With i in [1, ',num2str(N),'].']);

dist = 0.5*(x_l+x_r);
tilt2plot=tilt*180/pi;
timssss=times_l;

plot_time_v_dist(timssss, tilt2plot, dist, 'Tilt [$^\circ$]', 'horizontal position $x$');