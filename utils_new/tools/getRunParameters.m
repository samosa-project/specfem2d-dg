% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   [xminmax, zminmax, Xsource] = getRunParameters(parfile, sourcefile, interfacesfile)
% with:
%   TODO.
% yields:
%   TODO.

function [xminmax, zminmax, Xsource] = getRunParameters(parfile, sourcefile, interfacesfile)
  externalMesh = extractParamFromInputFile(parfile, 'read_external_mesh', 'bool');
  if(externalMesh)
    disp(['[',mfilename,'] Need to read an external mesh.']);
    P=extMesh_loadNodes();
    xminmax = [min(P(:,1)), max(P(:,1))];
    zminmax = [min(P(:,2)), max(P(:,2))];
  else
    xminmax = [extractParamFromInputFile(parfile, 'xmin', 'float'), extractParamFromInputFile(parfile, 'xmax', 'float')];
    [zmin, zmax] = extractZminZmaxFromInterfacesFile(interfacesfile);
    zminmax = [zmin,zmax];
  end
  Xsource = [extractParamFromInputFile(sourcefile, 'xs', 'float'), extractParamFromInputFile(sourcefile, 'zs', 'float')];
  disp(['[',mfilename,'] [xminmax zminmax] found to be [',num2str([xminmax, zminmax]),'].']);
end