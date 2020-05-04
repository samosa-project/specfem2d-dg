function [layers, interfaces, zmin, zmax] = readExampleFiles_meshfem_mesh(intfile)

  % intfile = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test__LNS_generalised__using_custom_fields/interfaces_input';

  ff = fopen(intfile);

  line = {};
  i = 1;
  finished = 0;
  while(not(finished))
    line{i} = fgetl(ff);
    if(numel(line{i})==0)
      % emtpy line, skip
      continue;
    end
    if(line{i}(1)=='#')
      % commented line, ignore
      continue;
    end
    if(isfloat(line{i}) & numel(line{i})==1)
      % EOF
      finished = 1;
    end
    i = i+1;
  end

  fclose(ff);

  line(end)=[]; % remove EOF
  
  for ll=1:numel(line)
    line{ll} = rmComEndl(line{ll});
  end
  
  nint = str2double(line{1});
  line(1)=[]; % remove

  interfaces = {};
  for i=1:nint
    interfaces{i}.npts = str2double(line{1}); line(1)=[];
    for pt = 1:interfaces{i}.npts
      interfaces{i}.xz(pt,1:2) = str2num(line{1}); line(1)=[];
    end
  end

  nlay = nint-1;

  layers = {};
  for l=1:nlay
    layers{l}.nz = str2double(line{1}); line(1)=[];
  end
  
  % find zmin zmax
  zmin=Inf;
  zmax=-Inf;
  for i=1:nint
    curz = interfaces{i}.xz(:,2);
    curmin = min(curz);
    curmax = max(curz);
    if(curmin<zmin)
      zmin = curmin;
    end
    if(curmax>zmax)
      zmax = curmax;
    end
  end
end

function out = rmComEndl(str)
  out = regexprep(str, '#.*', '');
end