% Author:        Léo Martire.
% Description:   Produces a variety of filters.
%
% Usage:
%   [freq, H_lp, H_hp, H_bp] = custom_Filter(f_sample, n_samples, ...
%                                            c_f_lp, c_f_hp, ...
%                                            d_lp, d_hp, order)
% with:
%   f_sample the sampling frequency [1/s],
%   n_samples the number of samples [1],
%   c_f_lp the low-pass cutoff frequency,
%   c_f_hp the high-pass cutoff frequency,
%   d_lp,
%   d_hp,
%   order the filter order (1, 2, or 4),
% yields:
%   freq a (n_samples) frequency array,
%   H_lp a low-pass filter,
%   H_hp a high-pass filter,
%   H_bp a band-pass filter.

function [freq, H_lp, H_hp, H_bp] = custom_Filter(f_sample, n_samples, c_f_lp, c_f_hp, d_lp, d_hp, order)
  freq = linspace(-f_sample/2, f_sample/2, n_samples);
  switch(order)
    case(1)
      H_hp = 1i*freq/c_f_hp./(1 +  1i*freq/c_f_hp);
      H_lp = 1./(1 +  1i*freq/c_f_lp);
    case(2)
      H_hp = -(freq/c_f_hp).^2 ./ (1 +  2*d_hp*1i*freq/c_f_hp - (freq/c_f_hp).^2);
      H_lp = 1./(1 +  2*d_lp*1i*freq/c_f_lp - (freq/c_f_lp).^2);
    case(4)
      H_hp = ((freq/c_f_hp).^2 ./ (1 +  2*d_hp*1i*freq/c_f_hp - (freq/c_f_hp).^2)).^2;
      H_lp = 1./(1 +  2*d_lp*1i*freq/c_f_lp - (freq/c_f_lp).^2).^2;
    otherwise
      error(['[',mfilename,', ERROR] Filter order unknown.']);
  end
  H_bp=H_lp.*H_hp;
end

