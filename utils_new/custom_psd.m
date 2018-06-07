% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [PSD, f] = custom_psd(t, s)
  s = detrend(s); % Remove eventual linear trend.
  dt = t(2) - t(1);
  fmax = 1 / dt;
  nfft = max(2 ^ floor(log2(length(s))), length(s));
  window = hann(nfft);
  noverlap = nfft / 2;
  [PSD, f] = pwelch(s, window, noverlap, nfft, fmax);
end