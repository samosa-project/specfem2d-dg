% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   [outputdata, nsamp] = readAndSubsampleSynth(OFd, stationGlobalNumber, unknown, extension, subsample, wanted_dt, outputSubsampleInfo)
% with:
%   outputSubsampleInfo
% yields:
%   TODO.

function [outputdata, nsamp] = readAndSubsampleSynth(OFd, stationGlobalNumber, channel, extension, subsample, wanted_dt, outputSubsampleInfo)
  if(not(exist('subsample'))); subsample=0; end;
  if(not(exist('wanted_dt'))); wanted_dt=-1; end;
  if(not(exist('outputSubsampleInfo'))); outputSubsampleInfo=1; end;
  % Read the synthetic.
%   stattag = stations_data.textdata(istat_glob, 1);
  stattag = sprintf('S%04d', stationGlobalNumber);
  OFd = char(OFd); channel = char(channel); extension = char(extension); % safety.
  file = [OFd, 'AA.', stattag, '.', channel, '.', extension];
%   data = load(file{1});
  data = load(file);
  nt = max(size(data));
  meanactualdt = mean(diff(data(:,1)));
  if (subsample == 1 && wanted_dt>meanactualdt)
    % Sub-sample of records.
    nsub = floor(wanted_dt/meanactualdt);
    nd = max(size(data(1:nsub:nt, 1)));
    if(outputSubsampleInfo == 1)
      disp(['[',mfilename,'] Subsampled synthetics by a factor ',num2str(nsub),'.']);
    end
  elseif (subsample == 1 && wanted_dt<meanactualdt)
    if(outputSubsampleInfo == 1)
      disp(['[',mfilename,'] Subsampled dt is smaller than actual dt, discarding subsampling.']);
    end
    nsub = 1;
    nd = nt;
  else
    nsub = 1;
    nd = nt;
  end
  outputdata(:, :) = data(1:nsub:nt, :);
  nsamp = size(outputdata, 1);
end