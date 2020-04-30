% BGMODEL = load_bg_model(DATAFILE, HEADER)

function BGMODEL = load_bg_model(DATAFILE, HEADER)
%   DATAFILE = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB/background_model.bin';
%   HEADER = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB/background_model_header.dat';
  
  if(not(exist('HEADER', 'var')))
    % try and guess header file.
    spl = split(DATAFILE, filesep);
    spl{end} = regexprep(spl{end}, '.bin', '_header.dat');
    HEADER = join(spl, filesep);
    HEADER = HEADER{1};
    if(not(exist(HEADER, 'file')))
      error('could not guess header file, or did not find it');
    end
  end
  
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
end