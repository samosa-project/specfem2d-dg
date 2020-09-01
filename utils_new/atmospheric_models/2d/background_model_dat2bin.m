% This script is aimed at converting a background model (supposedly under
% the right format) from ASCII to binary. Prefer using the other scripts,
% this method seems quite hefty.

%function background_model_dat2bin()
  input = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_lns_custom_wavefield/background_model.dat';
  output = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_lns_custom_wavefield/background_model.bin';
  
  prec = 'real*4';
  
  fin = fopen(input, 'r');
  fout = fopen(output, 'w');
  
  background_model = textscan(fin, '%f %f %f %f %f %f %f %f %f %f', 'CollectOutput', 1, 'headerlines', 3);
  background_model = background_model{1};
  fclose(fin);
  
  % TODO: replicate header
  
  % Print numbers to file.
  
  fwrite(fout, background_model, prec);
  
  fclose(fout);
  
  % test output
  fbin = fopen(output, 'r');
  A = fread(fbin, size(background_model), prec);
  reldiff = max(max(abs(A-background_model)./abs(background_model)));
  disp(['[',mfilename,'] Relative difference between input and output = ',sprintf('%.6f',reldiff*100),' %.']);
  
  
  fclose('all');
%end