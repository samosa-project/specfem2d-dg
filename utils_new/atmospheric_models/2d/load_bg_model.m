% BGMODEL = load_bg_model(DATAFILE, header_or_nlines)

function BGMODEL = load_bg_model(DATAFILE, header_or_nlines)
%   DATAFILE = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB/background_model.bin';
%   header_or_nlines = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNS_GWB/background_model_header.dat';
  
  if(not(exist('header_or_nlines', 'var')))
    % try and guess header file.
    spl = split(DATAFILE, filesep);
    spl{end} = regexprep(spl{end}, '.bin', '_header.dat');
    header_or_nlines = join(spl, filesep);
    header_or_nlines = header_or_nlines{1};
    if(not(exist(header_or_nlines, 'file')))
      error('could not guess header file, or did not find it');
    end
  end
  if(isfloat(header_or_nlines) & numel(header_or_nlines)==1)
    NLINES = header_or_nlines;
  elseif(ischar(header_or_nlines))
    % Extract number of lines.
    fhead = fopen(header_or_nlines,'r');
    NLINES = textscan(fhead, '%f');
    NLINES = NLINES{1};
    fclose(fhead);
  else
    error('kek');
  end

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