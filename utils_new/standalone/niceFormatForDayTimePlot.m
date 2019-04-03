% Author:        Léo Martire.
% Description:   Converts XTickLabels of some axis to a nice readable
%                format.
% Notes:         The axis' XTicks (or custom XTicks) must be in decimal
%                day format (e.g. 72.5 for midday of day 72, etc.)
%
% Usage:
%   niceFormatForDayTimePlot(axx)
%   niceFormatForDayTimePlot(axx, XTLRotation)
%   niceFormatForDayTimePlot(axx, XTLRotation, customXTicks)
% with:
%   axx          an active axis,
%   XTLRotation  an angle to rotate XTickLabels to (in [deg] clockwise),
%   customXTicks custom XTicks.
% yields:
%   N./A..
function niceFormatForDayTimePlot(axx, XTLRotation, customXTicks)
  if(nargin < 1)
    error(['[',mfilename,', ERROR] Need at least an axis.']);
  end
  if(not(exist('XTLRotation')))
    XTLRotation = 0;
  end
  if(not(exist('customXTicks')))
    customXTicks=-1;
  end
  
  if(not(numel(customXTicks==1) & customXTicks==-1))
    set(axx, 'xtick', customXTicks); % eventually set custom XTicks
  end

  XT=get(axx,'xtick'); % get current xticks

  DAYS=floor(XT); % extract days
  
  SECONDS = floor((XT-DAYS)*24*3600); % extract seconds of day

  SECONDS_STR = datestr(seconds(SECONDS), 'HH:MM'); % produce nice string for hours:minutes of day
  
  % produce nice xticklabels
  zfilllength = floor(log10(max(DAYS)))+1; % length of zfill adjusted to maximum value of days
  FULL_XTL = [];
  for d = 1:numel(DAYS)
    FULL_XTL(d,:) = [sprintf(['D%0',num2str(zfilllength),'.f'],DAYS(d)), '@', SECONDS_STR(d,:)];
  end
  FULL_XTL = char(FULL_XTL);

  set(axx,'xticklabels', FULL_XTL);
  
  set(axx,'xticklabelrotation',XTLRotation)
end