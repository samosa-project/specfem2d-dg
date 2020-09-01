prec = 'real*8';

output = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_lns_custom_wavefield_4elements_test_binary/background_model.bin';
output_head = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_lns_custom_wavefield_4elements_test_binary/background_model_bin_header.dat';


% values = [[1,2,3,4,5,6,7,8,9,10]]
values = [[1,2,3,4,5,6,7,8,9,10];...
          (10+[1,2,3,4,5,6,7,8,9,10])]
numlines = size(values, 1);

% first number should be an integer containing the number of "lines"
fout_hed = fopen(output_head, 'w');
fprintf(fout_hed, '%d', numlines);
fclose(fout_hed);

% then, write all values
fout = fopen(output, 'w');
fwrite(fout, values', prec);
fclose(fout);

% close all
fclose('all');