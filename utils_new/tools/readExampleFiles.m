% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   [xminmax, zminmax, Xsource] = readExampleFiles(parfile, sourcefile, interfacesfile)
% with:
%   TODO.
% yields:
%   TODO.

function [xminmax, zminmax, Xsource] = readExampleFiles(parfile, sourcefile, interfacesfile)
  externalMesh = readExampleFiles_extractParam(parfile, 'read_external_mesh', 'bool');
  
  if(externalMesh)
    disp(['[',mfilename,'] Need to read an external mesh.']);
    P=extMesh_loadNodes();
    xminmax = [min(P(:,1)), max(P(:,1))];
    zminmax = [min(P(:,2)), max(P(:,2))];
  else
    xminmax = [readExampleFiles_extractParam(parfile, 'xmin', 'float'), readExampleFiles_extractParam(parfile, 'xmax', 'float')];
    [zmin, zmax] = readExampleFiles_zMinMaxInterfacesFile(interfacesfile);
    zminmax = [zmin, zmax];
  end
  
  Xsource = [readExampleFiles_extractParam(sourcefile, 'xs', 'float'), readExampleFiles_extractParam(sourcefile, 'zs', 'float')];
  disp(['[',mfilename,'] [xminmax zminmax] found to be [',num2str([xminmax, zminmax]),'].']);
end