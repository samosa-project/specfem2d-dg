function [time, Zamp, NMs, COLs, LSs, ord, sqrd_L2_Norm_of_Err_FAF, sqrd_L2_Norm_of_Err_BUF] = ABC_load_TS_and_compute_error(OFs, nbstats, subsample, subsample_dt, errorFactor, colours_runs, DSPLNM_strct)
  
  % extract OFs and colours
  largeRunOF = OFs{1};
  farFieldRunOF = OFs{2};
  bufferRunOF = OFs{3};
  colour_LAR = colours_runs{1};
  colour_FAF = colours_runs{2};
  colour_BUF = colours_runs{3};
  
  % colours and linestyles
  LS_truth = '-';
  LS_tests = '--';
  LS_errors = ':';
  
  % Ordering for plot. Lowest is below (towards background), higher is above (towards foreground).
  ord = {};
  ord.EFf=0;
  ord.EBu=1;
  ord.La=2;
  ord.Ff=3;
  ord.Bu=4;
  
  % retrieve time series
  Zamp=[]; NMs={};
  for i=1:nbstats
    % load LARGE run (assumed to be truth)
    j=ord.La*nbstats+i; iLa=j;
    [data, ~] = readAndSubsampleSynth(largeRunOF, i, 'BXZ', 'semv', subsample, subsample_dt, 1); Zamp(j,:)=data(:,2)'; if(i==1); time=data(:,1)'; end;
    NMs{j} = ['S',num2str(i),' ',DSPLNM_strct.LAR,'']; COLs{j} = colour_LAR; LSs{j} = LS_truth;
    
    % load FARFIELD run (or any other test ABC method)
    j=ord.Ff*nbstats+i; iFAF=j;
    [data, ~] = readAndSubsampleSynth(farFieldRunOF, i, 'BXZ', 'semv', subsample, subsample_dt, 1); Zamp(j,:)=data(:,2)';
    NMs{j} = ['S',num2str(i),' ',DSPLNM_strct.FAF,'']; COLs{j} = colour_FAF; LSs{j} = LS_tests;
    j=ord.EFf*nbstats+i; % save error between the previous and LARGE
    Zamp(j,:) =  errorFactor* abs(Zamp(iLa,:) - Zamp(iFAF,:));
    NMs{j} = ['S',num2str(i),' $',num2str(errorFactor),'\times|$',DSPLNM_strct.LAR,'$-$',DSPLNM_strct.FAF,'$|$']; COLs{j} = colour_FAF; LSs{j} = LS_errors;
    
    % load BUFFER run (or any other test ABC method)
    j=ord.Bu*nbstats+i; iBUF=j;
    [data, ~] = readAndSubsampleSynth(bufferRunOF, i, 'BXZ', 'semv', subsample, subsample_dt, 1); Zamp(j,:)=data(:,2)';
    NMs{j} = ['S',num2str(i),' ',DSPLNM_strct.BUF,'']; COLs{j} = colour_BUF; LSs{j} = LS_tests;
    j=ord.EBu*nbstats+i; % save error between the previous and LARGE
    Zamp(j,:) =  errorFactor* abs(Zamp(iLa,:) - Zamp(iBUF,:));
    NMs{j} = ['S',num2str(i),' $',num2str(errorFactor),'\times|$',DSPLNM_strct.LAR,'$-$',DSPLNM_strct.BUF,'$|$']; COLs{j} = colour_BUF; LSs{j} = LS_errors;
  end
  
  sqrd_L2_Norm_of_Err_FAF = [];
  sqrd_L2_Norm_of_Err_BUF = [];
  for i = 1:nbstats
    sqrd_L2_Norm_of_Err_FAF(i) = trapz(time, (Zamp(ord.EFf*nbstats + i,:)/errorFactor).^2); % compute squared L2 norm
    sqrd_L2_Norm_of_Err_BUF(i) = trapz(time, (Zamp(ord.EBu*nbstats + i,:)/errorFactor).^2); % compute squared L2 norm
  end
end

