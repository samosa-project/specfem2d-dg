% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [signal_LP, signal_HP] = custom_filter(time, signal, fcut)
  signal = detrend(signal);
  dt = time(2) - time(1);
  Fls = 1 / dt;
  flmax = Fls / 2;
  Wln = fcut / flmax;
  order = 3;
  [bl, al] = butter(order, Wln, 'low');
  signal_LP = filtfilt(bl, al, signal);
  signal_HP = signal - signal_LP;
end