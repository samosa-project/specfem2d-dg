%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Attenuation File Analysis                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time, dataMIC, distMic_Spkr] = Load_Data_ATN(timestamp, nrun, pre, freq)
	do_plot = 0;

  % Selection of a file in the pathData directory
%   pathData = 'C:\Users\ba.chide\Desktop\Aarhus Data\Aarhus_Data_0410\';
  pathData = '/home/l.martire/Documents/data/aarhus/';
%   FileATN = ['MIC_ATN_20190410_165702_Run00322_Pre10,0 _Freq1324.tdms'];
  FileATN = ['MIC_ATN_20190410_',timestamp,'_Run',sprintf('%05d',nrun),'_Pre',regexprep(sprintf('%.1f',pre),'\.', ','),' _Freq',sprintf('%04d',freq),'.tdms'];
  
  % Parameters for each microphones.
  fs = 200e3; % sampling frequency [Hz]
  dt = 1/fs; % period [s]
  sensitivity = [1.34 1.37 1.36 1.39 1.49]; % conversion to pascal [V/Pa]
%   distMic_Spkr = [0.3, 0.5, 1, 2, 3]; % distance from the speaker [m]
  distMic_Spkr = 0.44 + [0, 0.26, 0.761, 1.76, 2.76]; % mail chide 200604@1600

  % Convert data into a .mat file and loading data
  matFile = [FileATN(1:end-5) '.mat'];
  if(exist([pathData, filesep, matFile], 'file') == 2)
    disp('This file has already been converted.');
  else
    simpleConvertTDMS([pathData, filesep, FileATN]);
  end
  matFile = matfile([pathData, filesep, matFile]);
%   matFile = save([pathData, filesep, matFile], '-v7.3');

  % Data extraction and loading into table dataMIC.
  for iMic = 1:5
%       disp(num2str(iMic));
      if iMic ==1
          dataTMP = matFile.d_Micro1(1,1);
      elseif iMic == 2
          dataTMP = matFile.d_Micro2(1,1);
      elseif iMic == 3
          dataTMP = matFile.d_Micro3(1,1);
      elseif iMic == 4
          dataTMP = matFile.d_Micro4(1,1);
      elseif iMic ==5
          dataTMP = matFile.d_Micro5(1,1);
      end
      dataMIC(:,iMic) = detrend(dataTMP.Data)/sensitivity(iMic);
  end % end of loop on microphones.

  % Produce time vector.
  time = [0:dt:(length(dataMIC(:,iMic))-1)*dt];
  
  if(do_plot)
    % Plot of time series.
    figure();
    plot(time, dataMIC(:,1));
    hold on;
    plot(time, dataMIC(:,2));
    plot(time, dataMIC(:,3));
    plot(time, dataMIC(:,4));
    plot(time, dataMIC(:,5));
    xlabel('Time [s]');
    ylabel('Sound Pressure [Pa]');
  end
end