% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Asks user if a packed data should be filtered. If so,
%                ask user parameters of filter, and filter all of them.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         Needs:
%                a) .m scripts and functions (if not alongside this script, recover via Léo):
%                  1) custom_Filter.m

function [data_v_filtered] = bulkfilter(data_t,data_v)
  possible_orders=[1,2,4]; % Possible (implemented) filter orders.
  
  % Make sure data has right shape.
  % The following lines assume we have more time steps than stations to plot.
  data_t=reshape(data_t,[min(size(data_t)), max(size(data_t))]);
  data_v=reshape(data_v,[min(size(data_v)), max(size(data_v))]);
  if(not(all(size(data_t)==size(data_v))))
    error(['[',mfilename,', ERROR] time data and amplitude should have the same size, but right now do not.']);
  end
  
  nbstat=size(data_t,1);
  
  filter_data = - 1;
  % while(not(ismember(filter_data,[0,1,2,3])))
  %   filter_data=input(['  Filter (0 for no, 1 for high-pass, 2 for low-pass, 3 for band-pass)? > ']);
  while (not(ismember(filter_data, [0, 1, 2, 3])))
    filter_data = input(['[', mfilename, '] Filter (0 for no, 1 for high-pass, 2 for low-pass, 3 for band-pass)? > ']);
  end
  if (filter_data ~= 0)
    filter_fchp = -1;
    filter_fclp = -1;
    filter_dhp = 0.65;
    filter_dlp = 0.65;
    filter_ord = -1;
    while (not(ismember(filter_ord, possible_orders)))
      filter_ord = input(['[', mfilename, '] Filter order (possible values [',num2str(possible_orders),'])? > ']);
    end
    if(ismember(filter_data, [1, 3]))
      while (filter_fchp <= 0)
        filter_fchp = input(['[', mfilename, '] High-pass cutoff frequency? > ']);
      end
    end
    if(ismember(filter_data, [2, 3]))
      while (filter_fclp <= 0)
        filter_fclp = input(['[', mfilename, '] Low-pass cutoff frequency? > ']);
      end
    end
    
    if(not(ismember(filter_ord, possible_orders)))
      error(['[', mfilename, ', ERROR] Filtering type not implemented.']);
    end
    
    for i = 1:nbstat
      switch filter_data
        case 1
          [filt_f,~,filt_hhp,~]=custom_Filter(mean(1./diff(data_t(i, :))),length(data_t(i, :)),filter_fclp,filter_fchp,filter_dlp,filter_dhp,filter_ord);
          data_v_filtered(i, :) = real(ifft(fft(detrend(data_v(i, :))) .* fftshift(filt_hhp)));
          disp(['[',mfilename,'] Highpass filtered data over ',num2str([filter_fchp]),' Hz with homemade order ',num2str(filter_ord),' filter.']);
        case 2
          [filt_f,filt_hlp,~,~]=custom_Filter(mean(1./diff(data_t(i, :))),length(data_t(i, :)),filter_fclp,filter_fchp,filter_dlp,filter_dhp,filter_ord);
          data_v_filtered(i, :) = real(ifft(fft(detrend(data_v(i, :))) .* fftshift(filt_hlp)));
          disp(['[',mfilename,'] Lowpass filtered data under ',num2str([filter_fclp]),' Hz with homemade order ',num2str(filter_ord),' filter.']);
        case 3
%           [~, data_v_HP] = custom_filter(data_t(i, :), data_v(i, :), filter_fcutoff);
          [filt_f,~,~,filt_hbp]=custom_Filter(mean(1./diff(data_t(i, :))),length(data_t(i, :)),filter_fclp,filter_fchp,filter_dlp,filter_dhp,filter_ord);
          data_v_filtered(i, :) = real(ifft(fft(detrend(data_v(i, :))) .* fftshift(filt_hbp)));
%           data_v(i, :) = data_v_HP;
          disp(['[',mfilename,'] Bandpass filtered data between (',num2str([filter_fchp,filter_fclp]),') Hz with homemade order ',num2str(filter_ord),' filter.']);
%           disp(['[', mfilename, ', WARNING] Data was high-pass filtered, with cutoff frequency ', num2str(filter_fchp), '.']);
%           clear('data_v_HP');
        otherwise
          error(['[', mfilename, ', ERROR] Filtering type not implemented.']);
      end
    end
  else
    data_v_filtered=data_v;
    disp(['[', mfilename, '] Data was not filtered.']);
  end
end