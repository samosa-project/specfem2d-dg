% Author:        Léo Martire.
% Description:   Wrapper for loading dumps from SPECFEM simulations and
%                converting those to LNS generalised background model
%                files, in binary format.
% Notes:         N. A.
%
% Usage:
%   write_bg_model_from_dumps(OFD, IT, outputFolder)
% with:
%   TODO.

function [] = write_bg_model_from_dumps(OFD, IT, outputFolder)
%   fileType = 'asc';
  fileType = 'bin';
  
  % Eventually re-interpolate on an uniform grid. Not advised.
  uniform.do = 0;
  uniform.nx = 4;
  uniform.nz = uniform.nx;
  
  % Load dumps.
  [ROWS, header] = dumps_to_bgmodel(OFD, IT, uniform);
  
  switch(fileType)
    case 'asc'
      header{1} = [header{1}, ' This file was generated using the Matlab script ''',mfilename('fullpath'),'.m''.'];
    case 'bin'
      header{numel(header) + 1} = ['This file was generated using the Matlab script ''',mfilename('fullpath'),'.m''.'];
    otherwise
  end
  
  % Actually write.
  write_bg_model(ROWS, 'fileType', fileType, 'outputFolder', outputFolder, 'header', header);
end

