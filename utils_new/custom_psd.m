% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [PSD, f, PIC] = custom_psd(t, s)
  s = detrend(s); % Remove eventual linear trend.
  dt = t(2) - t(1); % Sample rate.
  fmax = 1 / dt; % Maximum frequency.
  nfft = min(2 ^ floor(log2(length(s))), length(s)); % Maximum window.
%   nfft = nfft/4; % Window of interest.
%   nfft = nfft/16; % Window of interest.
  window = hann(nfft);
  noverlap = nfft / 2;
  [PSD, f, PIC] = pwelch(s, window, noverlap, nfft, fmax);
end
