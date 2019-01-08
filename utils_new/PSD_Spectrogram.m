% Author:        Alexandre Cadu.
% Mail:          ??
% Description:   Computes the power spectral density s_PSD of s and
%                associated frequencies f_PSD over N_window with
%                sampling frequency of fs. The overlapping is of 50 %.
%                The window is a square sine function
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [s_PSD, f_PSD, s_SPG, t_SPG, s_PSD_plottable, f_PSD_plottable] = PSD_Spectrogram(s, fs, N_window)
    N_t = length(s);        % Number of samples
    T_s = (N_t - 1)/fs;     % Acquisition time
    N_window_size = floor(N_t/N_window);    % Number of samples per window
    
    % PSD frequency table (Hz)
    f_PSD = linspace(-fs/2, fs/2, N_window_size);   
    
    % Signal PSD table
    s_PSD = zeros(1, N_window_size);
    
    % Spectrogram table (s)
    t_SPG = zeros(2*N_window-1, 1);
    
    % Signal spectrogram table
    s_SPG = zeros(2*N_window-1, N_window_size);   
    
    % Square sine window function
    window_table = sin(pi*(0:N_window_size-1)/(N_window_size-1));
    
    % Loop over windows
    for k_window = 1:2*N_window-1
        % Window index
        window_index_start = floor(0.5*(k_window - 1)*N_window_size + 1);
        window_index_stop = floor((0.5*k_window + 0.5)*N_window_size);
        
        % Extracted signal times window function
        s_windowed = s(window_index_start:window_index_stop) .* window_table;
        
        s_SPG(k_window, :) = abs(fft((s_windowed))).^2;
%         s_SPG(k_window, :) = abs(fft(detrend(s_windowed))).^2;
        
        t_SPG(k_window) = window_index_start/fs;
        
        % Cumulated PSD
        s_PSD = s_PSD + s_SPG(k_window, :);
    end
    
    % PSD Normalization
    s_PSD = fftshift(s_PSD)/(fs*N_window_size*(2*N_window - 1)*0.5);
    
    % SPG Normalization
    s_SPG = fftshift(s_SPG, 2)/(fs*N_window_size*0.5);
    
    % Prepare PSD for direct plotting.
    N_PSD = floor(length(s_PSD)/2 + 2);
    f_PSD_plottable=f_PSD(N_PSD:end);
    s_PSD_plottable=sqrt(s_PSD(N_PSD:end));
end