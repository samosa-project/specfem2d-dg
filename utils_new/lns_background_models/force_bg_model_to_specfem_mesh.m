function [newmodel, oldmodel] = force_bg_model_to_specfem_mesh(EXAMPLE_DIR)

  % EXAMPLE_DIR = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test__LNS_generalised__using_custom_fields';

  if(not(EXAMPLE_DIR(end)==filesep))
    EXAMPLE_DIR = [EXAMPLE_DIR, filesep];
  end
  oldmodel = [EXAMPLE_DIR,filesep,'background_model.bin'];

  % Load it.
  bgm = load_bg_model(oldmodel);
  [bgm, isMatrix] = try_make_bg_model_matrix(bgm);

  % Grab SPECFEM2D mesh.
  parfile = [EXAMPLE_DIR, 'parfile_input'];
  intfile = [EXAMPLE_DIR, 'interfaces_input'];
  read_external_mesh = readExampleFiles_extractParam(parfile, 'read_external_mesh', 'bool');
  if(read_external_mesh)
    error('not implemented');
  else
    xmin = readExampleFiles_extractParam(parfile, 'xmin', 'float');
    xmax = readExampleFiles_extractParam(parfile, 'xmax', 'float');
    NX = readExampleFiles_extractParam(parfile, 'nx', 'int') + 1;
    [layers, ~, zmin, zmax] = readExampleFiles_meshfem_mesh(intfile);
    if(numel(layers)>1)
      error('not implemented');
    end
    NZ = layers{1}.nz + 1;
    x = gllify(linspace(xmin, xmax, NX)); z = gllify(linspace(zmin, zmax, NZ));
    z = z - min(z); % model should start at z=0
    [X, Z] = meshgrid(x, z);

    if(isMatrix)
      [order, ~, ~, ~] = order_bg_model();
      nb_qty = size(order, 1);
%       [min(X(:)), max(X(:)), min(Z(:)), max(Z(:))]
      for iqty = 3:nb_qty
        bgm.(order(iqty, :)) = interp2(bgm.(order(1,:)), bgm.(order(2,:)), bgm.(order(iqty,:)), X, Z);
      end
      bgm.(order(1,:)) = X;
      bgm.(order(2,:)) = Z;

      % Format under large unit table.
      npts = numel(bgm.(order(1,:))(:));
      ROWS = zeros(npts, nb_qty);
      for iqty = 1:nb_qty
        curQty = order(iqty,:);
        id = find(all((order==curQty)'));
        ROWS(:, id) = reshape(bgm.(curQty), npts, 1);
      end
      
    else
      error('not implemented');
    end
  end
  
  % Ask for confirmation before writing.
  estimate_size = numel(ROWS)*8/1048576;
  if(estimate_size>50)
    userIsSure = -1;
    while(not(numel(userIsSure)==1 && ismember(userIsSure, [0,1])))
      userIsSure = input(['[',mfilename,'] You are about to write ',num2str(estimate_size),' Mb (',num2str(numel(ROWS)),' floating point numbers) to disk. Are you sure? > ']);
    end
  else
    userIsSure = 1;
  end
  if(userIsSure)
    copyfile(oldmodel, [regexprep(oldmodel, '\.bin', '_ORIGINAL.bin')]); % save original
    newmodel = write_bg_model(ROWS, 'outputFolder', EXAMPLE_DIR);
  else
    error('aborted');
  end
end