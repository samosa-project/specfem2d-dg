function [outputArg1,outputArg2] = ABC_print_params(curCase_INFOS, outputDir)
  fp = fopen([outputDir,curCase_INFOS.curCase,'_parameters.txt'], 'w');
  
  fprintf(fp, ['domain of interest : ', sprintf('%7.2f ', curCase_INFOS.EBox),' (format is [xmin xmax zmin zmax])\n']);
  fprintf(fp, ['LAR domain:          ', sprintf('%7.2f ', curCase_INFOS.LAR.domain),'\n']);
  fprintf(fp, ['FAF domain:          ', sprintf('%7.2f ', curCase_INFOS.FAF.domain),'\n']);
  fprintf(fp, ['BUF domain:          ', sprintf('%7.2f ', curCase_INFOS.BUF.domain),'\n']);
  fprintf(fp, ['BUF length:          ', sprintf('%7.2f ', curCase_INFOS.BUF.length),'\n']);
  fprintf(fp, ['BUF epsilon:            ', sprintf('%.6e ', curCase_INFOS.BUF.epsilon),'\n']);
  fprintf(fp, ['BUF p:                  ', sprintf('%.6e ', curCase_INFOS.BUF.p),'\n']);
  fprintf(fp, ['BUF q:                  ', sprintf('%.6e ', curCase_INFOS.BUF.q),'\n']);
  
  fclose('all');
end

