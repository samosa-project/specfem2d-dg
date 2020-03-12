% Author:        Léo Martire.
% Description:   Wrapper for loading dumps from SPECFEM simulations and
%                converting those to LNS generalised background model
%                files, in binary format.
% Notes:         Needs scripts:
%                  utils_new/lns_background_models/dumps_to_bgmodel.m
%                  utils_new/lns_background_models/write_bg_model.m
%
% Usage:
%   write_bg_model_from_dumps(OFD, IT, outputFolder)
% with:
%   TODO.

function [] = write_bg_model_from_dumps(OFD, IT, outputFolder)
  % Eventually re-interpolate on an uniform grid. Not advised.
  uniform.do = 0;
  uniform.nx = 4;
  uniform.nz = uniform.nx;
  
  % Load dumps.
  [ROWS] = dumps_to_bgmodel(OFD, IT, uniform);
  
  header = {};
  header{1} = ['This file was generated using the Matlab script ''',mfilename('fullpath'),'.m''.\n'];
  header{2} = ['This background model is generated from SPECFEM dumps located in ',OFD,' at IT=',num2str(IT),'.\n'];
  
  % Actually write.
  write_bg_model(ROWS, 'fileType', 'bin', 'outputFolder', outputFolder, 'header', header);
end

