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
if (not(exist('synth_load_was_ran') && synth_load_was_ran == 1))
  error(['[', mfilename, ', ERROR] synth_load was not ran before.']);
end

Zamp_save=Zamp;
Ztime_save=Ztime;

v=Zamp;
t=Ztime;

% Select intervals to force at zero.
i=1;
tltu{i} = [];
i=2;
tltu{i} = [300,Inf];
i=3;
tltu{i} = [310, 340; ...
           430, Inf];

% Actually zero the thing.
for i=1:3
  loc_t=t(i,:);
  loc_v=v(i,:);
  loc_tltu=tltu{i};
  for j=1:size(loc_tltu,1)
    sel = find(loc_t>=loc_tltu(j,1) & loc_t<=loc_tltu(j,2)); % select ids wrt current t
    
    t_sta=loc_t(sel(1));
    t_end=loc_t(sel(end));
    v_sta=loc_v(sel(1));
    v_end=loc_v(sel(end));
    
    loc_v(sel) = interp1([t_sta, t_end], [v_sta, v_end], loc_t(sel)); % linear interpolation bewteen the two external points of the interval.
  end
  t(i,:)=loc_t;
  v(i,:)=loc_v;
end

disp(['[', mfilename, ', WARNING] Artifically set some intervals of signal to zero. Be careful of spurious behaviour.']);
disp(['[', mfilename, ', WARNING] Run the following command to make sure you have understood and replace loaded Ztime and Zamp with new arrays:']);
disp([blanks(length(mfilename)+2),'            Ztime=t; Zamp=v;']);