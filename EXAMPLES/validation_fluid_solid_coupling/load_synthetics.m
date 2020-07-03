function [time, amp] = load_synthetics(OFD, parfile, station_ids)
  subsample = 0; wanted_dt = -1; % do not subsample synthetics
  type_display = readExampleFiles_extractParam(parfile,'seismotype','int');
  [extension, ~] = getUnknowns(type_display, 'BXZ');
  % Load synthetics.
  nstat = numel(station_ids);
  for istat = 1:nstat
    istatglob = station_ids(istat);
    [data, nsamples] = readAndSubsampleSynth(OFD, istatglob, 'BXZ', extension, subsample, wanted_dt, istatglob);
    Ztime(istat, 1:nsamples) = data(:, 1)';
    Zamp(istat, 1:nsamples) = data(:, 2)';
    [data, nsamples] = readAndSubsampleSynth(OFD, istatglob, 'BXX', extension, subsample, wanted_dt, istatglob);
    Xtime(istat, 1:nsamples) = data(:, 1)';
    Xamp(istat, 1:nsamples) = data(:, 2)';
  %   max(abs(Xamp(istat,:)-Zamp(istat,:)))

    if(all(abs(Xamp(istat,1:nsamples)-Zamp(istat,1:nsamples))==0))
      % channel X is exactly channel Z, this means we have loaded a DG station
      % thus save either one, that is pressure
      time(istat, 1:nsamples) = Xtime(istat, 1:nsamples);
      amp(istat, 1:nsamples) = Xamp(istat, 1:nsamples);
    else
      % else, we have loaded a v station
      % thus, save velocity norm
      time(istat, 1:nsamples) = Xtime(istat, 1:nsamples);
      amp(istat, 1:nsamples) = (Xamp(istat, 1:nsamples).^2 + Zamp(istat, 1:nsamples).^2).^0.5;
    end

    clear('Xtime', 'Xamp', 'Ztime', 'Zamp');
  end
end