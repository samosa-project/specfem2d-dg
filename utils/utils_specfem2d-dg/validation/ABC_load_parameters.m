function [INFO_all] = ABC_load_parameters(largeRunOF, farFieldRunOF, bufferRunOF, curCase, BUF_params)
  % safeguard
  if(numel(bufferRunOF)>0); if(not(strcmp(bufferRunOF(end),filesep))); bufferRunOF=[bufferRunOF,filesep]; end; end;
  if(numel(largeRunOF)>0); if(not(strcmp(largeRunOF(end),filesep))); largeRunOF=[largeRunOF,filesep]; end; end;
  if(numel(farFieldRunOF)>0); if(not(strcmp(farFieldRunOF(end),filesep))); farFieldRunOF=[farFieldRunOF,filesep]; end; end;
  % retrieve parameters
  OF=largeRunOF; [xmM_La, zmM_La, ~] = readExampleFiles([OF,'input_parfile'], [OF,'SOURCE'], [OF,'input_interfaces']);
  OF=bufferRunOF; [xmM_Bu, zmM_Bu, ~] = readExampleFiles([OF,'input_parfile'], [OF,'SOURCE'], [OF,'input_interfaces']); LBUF=readExampleFiles_extractParam([OF,'input_parfile'], 'ABC_STRETCH_TOP_LBUF', 'float');
  OF=farFieldRunOF; [xmM_Fa, zmM_Fa, ~] = readExampleFiles([OF,'input_parfile'], [OF,'SOURCE'], [OF,'input_interfaces']);
  % prepare information strings
  domain_LAR_tex = ['$[',num2str(min(xmM_La)),', ',num2str(max(xmM_La)),']\times[',num2str(min(zmM_La)),', ',num2str(max(zmM_La)),']$'];
  domain_BUF_tex = ['$[',num2str(min(xmM_Bu)),', ',num2str(max(xmM_Bu)),']\times[',num2str(min(zmM_Bu)),', ',num2str(max(zmM_Bu)),']$'];
  domain_FAF_tex = ['$[',num2str(min(xmM_Fa)),', ',num2str(max(xmM_Fa)),']\times[',num2str(min(zmM_Fa)),', ',num2str(max(zmM_Fa)),']$'];
  domain_LAR = [xmM_La, zmM_La];
  domain_BUF = [xmM_Bu, zmM_Bu];
  domain_FAF = [xmM_Fa, zmM_Fa];
  % tit_La = ['$z\in[',num2str(min(zmM_La)),', ',num2str(max(zmM_La)),']$'];
  % tit_Bu = ['$z\in[',num2str(min(zmM_Bu)),', ',num2str(max(zmM_Bu)),']$, Buffer=',num2str(LBUF),', ',suffixBu];
  % tit_Fa = ['$z\in[',num2str(min(zmM_Fa)),', ',num2str(max(zmM_Fa)),']$'];
  % get energyboxes
  ENERGYBOX_La = [readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_XMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_XMAX', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMAX', 'float')];
  ENERGYBOX_Bu = [readExampleFiles_extractParam([bufferRunOF,'input_parfile'], 'ENERGYBOX_XMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_XMAX', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMAX', 'float')];
  ENERGYBOX_Fa = [readExampleFiles_extractParam([farFieldRunOF,'input_parfile'], 'ENERGYBOX_XMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_XMAX', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMIN', 'float'),readExampleFiles_extractParam([largeRunOF,'input_parfile'], 'ENERGYBOX_ZMAX', 'float')];
  if(not(all(ENERGYBOX_La==ENERGYBOX_Bu & ENERGYBOX_Bu==ENERGYBOX_Fa)))
    error(['ENERGY BOXES WERE NOT THE SAME, ABORTING']);
  else
    ENERGYBOX = ENERGYBOX_La;
    globxmin=min([min(xmM_La),min(xmM_Bu),min(xmM_Fa)]);
    globxmax=max([max(xmM_La),max(xmM_Bu),max(xmM_Fa)]);
    globzmin=min([min(zmM_La),min(zmM_Bu),min(zmM_Fa)]);
    globzmax=max([max(zmM_La),max(zmM_Bu),max(zmM_Fa)]);
    if(ENERGYBOX(1)<globxmin)
      ENERGYBOX(1)=globxmin;
    end
    if(ENERGYBOX(2)>globxmax)
      ENERGYBOX(2)=globxmax;
    end
    if(ENERGYBOX(3)<globzmin)
      ENERGYBOX(3)=globzmin;
    end
    if(ENERGYBOX(4)>globzmax)
      ENERGYBOX(4)=globzmax;
    end
  end
  
%   INFO_all = {[DSPLNM_strct.LAR,'=[',domain_LAR,'], ',DSPLNM_strct.FAF,'=[',domain_FAF,']'],[DSPLNM_strct.BUF,'=[',domain_BUF,'], '],['energy box: ',['$[',num2str(min(ENERGYBOX(1:2))),', ',num2str(max(ENERGYBOX(1:2))),']\times[',num2str(min(ENERGYBOX(3:4))),', ',num2str(max(ENERGYBOX(3:4))),']$']]};
  
  INFO_all.curCase = curCase;
  INFO_all.LAR.domain = domain_LAR;
  INFO_all.FAF.domain = domain_FAF;
  INFO_all.BUF.domain = domain_BUF;
  INFO_all.BUF.length = LBUF;
  INFO_all.BUF.epsilon = BUF_params(1);
  INFO_all.BUF.p = BUF_params(2);
  INFO_all.BUF.q = BUF_params(3);
  INFO_all.EBox = ENERGYBOX_La;
end

