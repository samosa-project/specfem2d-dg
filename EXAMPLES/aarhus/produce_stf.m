% see also /home/l.martire/Documents/work/mars/attenuation_aarhus

function produce_stf(f0)
  do_compare = 0;
  do_plot = 1;
  
  addpath(genpath('/home/l.martire/Documents/SPECFEM/specfem-dg-master/utils_new/tools'));
  thisFolder = [regexprep(mfilename('fullpath'),mfilename,'')];
  parfile = [thisFolder,filesep,'parfile_input'];
  sourcefile = [thisFolder,filesep,'source_input'];
  
  outfile = [thisFolder,filesep,'source.stf'];
  
%   f0 = 2090;
  % shift = 9625e-8;
%   shift = 0;
%   shift = 1/f0;
%   shift = 0.5/f0;
  shift = 0.25/f0;
  % shift = 0.2/f0;

  NSTEP = readExampleFiles_extractParam(parfile, 'NSTEP', 'int');
  dt = readExampleFiles_extractParam(parfile, 'DT', 'float');
  f0source = readExampleFiles_extractParam(sourcefile, 'f0', 'float');

  it = 0:NSTEP;
  t = it*dt;

  % produce sinus
  sinus = sin((t-shift)*2*pi*f0);

  % tukey
  shift_i = floor((shift)/dt);
  L = ceil((3/f0)/dt);
  tuk = 0*t;
  tuk(shift_i+[0:L-1]) = tukeywin(L, 0.6);

  % combine
  v = sinus.*tuk;

  if(do_compare)
    % compare
    load('source.stf_old.mat'); % for comparison
    figure();
    plot(t, v); hold on;
    plot(t_old, v_old, ':');
    xlim([0.4,2.4]*1e-3);
  end

  if(do_plot)
    % plot
    figure();
    plot(t, sinus); hold on;
    plot(t, tuk); hold on;
    plot(t, v);
    xlim([0, 1/f0 + shift+3/f0 + 1/f0]);
  end
  
  t0_specfem = -0.000574162679426;
  t0_specfem = -1.2/f0source;
  
  fid = fopen(outfile, 'w');
  prec = 15;
  fprintf(fid, ['%.',num2str(prec),'g %.',num2str(prec),'g\n'], [t+t0_specfem; v]);
  fclose(fid);
end
