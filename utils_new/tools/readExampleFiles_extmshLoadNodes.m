% Author:        Léo Martire.
% Description:   Loads nodes from a 'Nodes_extMesh' file.
% Notes:         The 'Nodes_extMesh' file only exists after having ran the
%                GMSH to SPECFEM mesh converter.
%
% Usage:
%   [P] = extMesh_loadNodes(Nodes_extMesh)
% with:
%   TODO.
% yields:
%   TODO.

function [P] = readExampleFiles_extmshLoadNodes(Nodes_extMesh)
  if(not(exist('Nodes_extMesh')))
    Nodes_extMesh = input(['[',mfilename,'] Path to Nodes_extMesh? > '],'s');
  end
  Nodes_extMesh_f = fopen(Nodes_extMesh,'r');
  lines = textscan(Nodes_extMesh_f, '%f %f');
%   X = lines{1}(2:end); Z = lines{2}(2:end); P = [X,Z];
  P = [lines{1}(2:end), lines{2}(2:end)];
  if( not(lines{1}(1)==size(P, 1)))
    error(['[',mfilename,', ERROR] Issue while reading Nodes_extMesh: number of points read (number of lines) does not correspond to number of points specified (first line).']);
  end
end