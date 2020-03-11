function [time, Zamp, NMs, COLs, LSs] = ABC_load_TS_and_compute_error(OFs, nbstats, subsample, subsample_dt, errorFactor, DSPLNM_strct)
  
  % extract OFs
  largeRunOF = OFs{1};
  farFieldRunOF = OFs{2};
  bufferRunOF = OFs{3};
  
  % colours and linestyles
  LS_truth = '-';
  LS_tests = '--';
  LS_errors = ':';
  colour_LAR = 'k';
  colour_FAF = 'c';
  colour_BUF = 'm';
  
  % Ordering for plot. Lowest is below (towards background), higher is above (towards foreground).
  orEFf=0;
  orEBu=1;
  orLa=2;
  orFf=3;
  orBu=4;
  
  % retrieve time series
  Zamp=[]; NMs={};
  for i=1:nbstats
    % load LARGE run (assumed to be truth)
    j=orLa*nbstats+i; iLa=j;
    [data, ~] = readAndSubsampleSynth(largeRunOF, i, 'BXZ', 'semv', subsample, subsample_dt, 1); Zamp(j,:)=data(:,2)'; if(i==1); time=data(:,1)'; end;
    NMs{j} = ['S',num2str(i),' ',DSPLNM_strct.LAR,'']; COLs{j} = colour_LAR; LSs{j} = LS_truth;
    
    % load FARFIELD run (or any other test ABC method)
    j=orFf*nbstats+i; iFAF=j;
    [data, ~] = readAndSubsampleSynth(farFieldRunOF, i, 'BXZ', 'semv', subsample, subsample_dt, 1); Zamp(j,:)=data(:,2)';
    NMs{j} = ['S',num2str(i),' ',DSPLNM_strct.FAF,'']; COLs{j} = colour_FAF; LSs{j} = LS_tests;
    j=orEFf*nbstats+i; % save error between the previous and LARGE
    Zamp(j,:) =  errorFactor* abs(Zamp(iLa,:) - Zamp(iFAF,:));
    NMs{j} = ['S',num2str(i),' $',num2str(errorFactor),'\times|$',DSPLNM_strct.LAR,'$-$',DSPLNM_strct.FAF,'$|$']; COLs{j} = colour_FAF; LSs{j} = LS_errors;
    
    % load BUFFER run (or any other test ABC method)
    j=orBu*nbstats+i; iBUF=j;
    [data, ~] = readAndSubsampleSynth(bufferRunOF, i, 'BXZ', 'semv', subsample, subsample_dt, 1); Zamp(j,:)=data(:,2)';
    NMs{j} = ['S',num2str(i),' ',DSPLNM_strct.BUF,'']; COLs{j} = colour_BUF; LSs{j} = LS_tests;
    j=orEBu*nbstats+i; % save error between the previous and LARGE
    Zamp(j,:) =  errorFactor* abs(Zamp(iLa,:) - Zamp(iBUF,:));
    NMs{j} = ['S',num2str(i),' $',num2str(errorFactor),'\times|$',DSPLNM_strct.LAR,'$-$',DSPLNM_strct.BUF,'$|$']; COLs{j} = colour_BUF; LSs{j} = LS_errors;
  end
end

