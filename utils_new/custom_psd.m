% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [Powerf, Freqf] = custom_psd(time, signal)
  signal = detrend(signal); % Remove eventual linear trend.
  dt = time(2) - time(1);
  Fls = 1 / dt;
  nfft = max(2 ^ floor(log2(length(signal))), length(signal));
  window = hann(nfft);
  noverlap = nfft / 2;
  [Powerf, Freqf] = pwelch(signal, window, noverlap, nfft, Fls);
end