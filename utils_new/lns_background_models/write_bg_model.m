% Author:        Léo Martire.
% Description:   Writes a LNS generalised background model array (which
%                shape is dictated by the script 'order_bg_model.m') to a
%                file for reading by SPECFEM-LNS.
% Notes:         The array should be prepared according to the order
%                dictacted by the script 'order_bg_model.m', preferably by
%                using it, or coding somthing similar to what is
%                implemented in 'dumps_to_bgmodel.m'.
%                In the case of binary format, a second file used for
%                header is produced. It is necessary for SPECFEM-LNS to
%                know how many points are exported.
%
% Usage:
%   write_bg_model(ROWS, 'fileType', 'bin', 'outputFolder', outputFolder, 'header', header);
% with:
%   ROWS                 a LNS generalised background model array (which
%                        shape is dictated by the script
%                        'order_bg_model.m')
%   [optional arguments],
%      fileType          the type of file (either 'binary' or 'ascii',
%                        defaults to 'binary'),
%      outputFolder      the folder in which to produce the files
%                        (defaults to ./),
%      header            some header information under cell format (should
%                        not contain more than 2 cells for ascii format,
%                        is unlimited for binary format, defaults to {}).

function [outpath] = write_bg_model(varargin)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PARAMETER READING AND SETTING.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nargs_needed = 1;
  if(nargin<nargs_needed)
    help(mfilename);
    error('');
  else
    narg_option_provided = nargin - nargs_needed;
    if(mod(narg_option_provided, 2))
      narg_option_provided
      error('must provide optional arguments as couples');
    else
      narg_option_couples_provided = narg_option_provided/2;
    end
  end
  % Read main parameters.
  i = 1;
  ROWS = varargin{i}; i=i+1;
  % Read optional parameters.
  fileType_fulfilled = 0;
  outputFolder_fulfilled = 0;
  header_fulfilled = 0;
  dummy2_fulfilled = 0;
  for c = 1:narg_option_couples_provided
    cur_name = varargin{nargs_needed + 2*(c-1) + 1};
    cur_val = varargin{nargs_needed + 2*c};
    switch(lower(cur_name))
      case 'filetype'
        fileType = cur_val;
        fileType_fulfilled = 1;
      case 'outputfolder'
        outputFolder = cur_val;
        outputFolder_fulfilled = 1;
      case 'header'
        header = cur_val;
        header_fulfilled = 1;
      case 'dummy2'
        dummy2 = cur_val;
        dummy2_fulfilled = 1;
      otherwise
        error('parameter unknown');
    end
  end
  % Default values.
  if(not(fileType_fulfilled))
    fileType = 'bin';
  end
  if(not(outputFolder_fulfilled))
    outputFolder = './';
  end
  if(not(header_fulfilled))
    header = {};
  end
  if(not(dummy2_fulfilled))
    dummy2 = 'dumdum';
  end
  % Check.
  switch(fileType)
    case {'asc', 'ascii', 'bin', 'binary'}
      % Ok.
    otherwise
      error(['[',mfilename,', ERROR] fileType unknown.']);
  end
  if(not(outputFolder(end) == filesep))
    outputFolder = [outputFolder, filesep];
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % START ACTUAL FUNCTION.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Precision.
  [asc__nb_significantdigits, bin__precision] = bg_model_parameters();
  
  % File names. See 'specfem2D_par_lns.f90'.
  asc__output_filename = 'background_model.dat';
  bin__output_filename = 'background_model.bin';
  bin__output_header_filename = 'background_model_header.dat';
  asc__output_path = [outputFolder, asc__output_filename];
  bin__output_path = [outputFolder, bin__output_filename];
  bin__output_header_path = [outputFolder, bin__output_header_filename];
  
  switch(fileType)
    case {'asc', 'ascii'}
      nb_qty = size(ROWS, 2);
      
      % Prepare ASCII format.
      asc__spacing = asc__nb_significantdigits+7;
      asc__format_number = ['%',num2str(asc__spacing),'.',num2str(asc__nb_significantdigits),'e '];
      asc__format_line = repmat(asc__format_number, [1,nb_qty]);
      asc__format_line = [asc__format_line(1:end-1), '\n'];

      % Open file and write in it.
      foutput = fopen(asc__output_path, 'w');

      % header
      if(numel(header)>2)
        error(['[',mfilename,', ERROR] With ASCII files, you must provide exactly 2 lines of header.']);
      end
      fprintf(foutput, [header{1}, '\n']);
      fprintf(foutput, [header{2}, '\n']);
      [~, lab] = order_bg_model();
      for q = 1:nb_qty
        fprintf(foutput, pad(lab{q}, asc__spacing+1));
      end
      fprintf(foutput, '\n');
      fprintf(foutput, asc__format_line, ROWS');
      fclose(foutput);
      disp(['[',mfilename,'] Finished writing background model to ASCII file (''',asc__output_path,''').']);
      
    case {'bin', 'binary'}
      numlines = size(ROWS, 1);

      % First line should be an integer containing the number of "lines".
      fout_hed = fopen(bin__output_header_path, 'w');
      fprintf(fout_hed, '%d\n', numlines);
      for i = 1:numel(header)
        fprintf(fout_hed, [header{i}, '\n']);
      end
      fclose(fout_hed);

      % Then, write all values sequentially.
      outpath = bin__output_path;
      fout = fopen(bin__output_path, 'w');
      fwrite(fout, ROWS', bin__precision);
      fclose(fout);
      disp(['[',mfilename,'] Finished writing background model to BIN file (''',bin__output_path,'''), with separate header (''',bin__output_header_path,''').']);
      
    otherwise
      error(['[',mfilename,', ERROR] File type unknown.']);
  end
  
  fclose('all'); % safety
  
end

