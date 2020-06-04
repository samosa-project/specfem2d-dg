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

function [xminmax, zminmax, Xsource] = readExampleFiles(parfile, sourcefile, interfacesfile, Nodes_extMesh, verbose)
  if(not(exist('Nodes_extMesh','var')))
    Nodes_extMesh=[];
  end
  if(not(exist('verbose','var')))
    verbose=0;
  end
  
  if(not(isempty(parfile)) && not(exist(parfile,'file')))
    error(['[',mfilename,', ERROR] Parfile ''',parfile,''' does not exist.']);
  end
  if(not(isempty(sourcefile)) && not(exist(sourcefile,'file')))
    error(['[',mfilename,', ERROR] Source file ''',sourcefile,''' does not exist.']);
  end
  if(not(isempty(interfacesfile)) && not(exist(interfacesfile,'file')))
    error(['[',mfilename,', ERROR] Interfaces file ''',interfacesfile,''' does not exist.']);
  end
  
  % default values for cases in which the user only asks one of the things
  xminmax=nan;
  zminmax=nan;
  Xsource=nan;
  
  if(not(isempty(parfile)))
    externalMesh = readExampleFiles_extractParam(parfile, 'read_external_mesh', 'bool');

    if(externalMesh)
      if(verbose)
        disp(['[',mfilename,'] Need to read an external mesh.']);
      end
      P=readExampleFiles_extmshLoadNodes(Nodes_extMesh);
      xminmax = [min(P(:,1)), max(P(:,1))];
      zminmax = [min(P(:,2)), max(P(:,2))];
    else
      xminmax = [readExampleFiles_extractParam(parfile, 'xmin', 'float'), readExampleFiles_extractParam(parfile, 'xmax', 'float')];
  %     [zmin, zmax] = readExampleFiles_zMinMaxInterfacesFile(interfacesfile);
  %     zminmax = [zmin, zmax];
%       [zminmax(1), zminmax(2)] = readExampleFiles_zMinMaxInterfacesFile(interfacesfile);
      [~, ~, zminmax(1), zminmax(2)] = readExampleFiles_meshfem_mesh(interfacesfile);
    end
  end
  
  Xsource = [readExampleFiles_extractParam(sourcefile, 'xs', 'float'), readExampleFiles_extractParam(sourcefile, 'zs', 'float')];
  
  if(verbose)
    disp(['[',mfilename,'] [xminmax zminmax] found to be [',num2str([xminmax, zminmax]),'].']);
  end
end