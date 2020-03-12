function [] = write_bg_model_from_dumps(OFD, IT, outputFolder)
  % Eventually re-interpolate on an uniform grid. Not advised.
  uniform.do = 0;
  uniform.nx = 4;
  uniform.nz = uniform.nx;
  
  % Load dumps.
  [ROWS] = dumps_to_bgmodel(OFD, IT, uniform);
  pause;
  
  header = {};
  header{1} = 'custom/default LNS generalised background model\n';
  header{2} = ['using dumps from IT=',num2str(IT),' in ',OFD,'\n'];
  
  % Actually write.
  write_bg_model(ROWS, 'fileType', 'bin', 'outputFolder', outputFolder, 'header', header);
end

