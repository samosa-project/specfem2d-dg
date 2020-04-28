DATAFILE = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB/background_model.bin';
HEADER = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB/background_model_header.dat';

% Extract number of lines.
fhead = fopen(HEADER,'r');
NLINES = textscan(fhead, '%f');
NLINES = NLINES{1};
fclose(fhead);

% Extract data.
fbody = fopen(DATAFILE, 'r');
[asc__nb_significantdigits, bin__precision] = bg_model_parameters();
data = fread(fbody, bin__precision);
data = reshape(data, [numel(data)/NLINES, NLINES])';
fclose(fbody);

% Put in nice format.
[order, ~] = order_bg_model();
nb_qty = size(order, 1);

BGMODEL = struct();
for iqty = 1:nb_qty
  qty = order(iqty, :);
  BGMODEL.(qty) = data(:, iqty);
end