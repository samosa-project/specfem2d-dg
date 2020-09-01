function pos_sources = loadSources(OFdir)
  pos_sources = [inf, inf]; % Allocate a row for the first source's position.
  % fid = fopen([rootd, '/DATA/SOURCE']);
  fid = fopen([OFdir, 'SOURCE']);
  if (fid == - 1)
    fid = fopen([OFdir, 'input_source']);
    if (fid == - 1)
      error(['[',mfilename,', ERROR] Cannot open SOURCE file (', OFdir, 'SOURCE', ').']);
    end
  end
  line = 0; xfound = 0; zfound = 0; stop = 0;
  while (stop == 0)
    % TODO: Loop on source number.
    line = fgetl(fid);
    if length(line) > 0
      if (line == - 1)
        stop = 1;
      end
      line = regexprep(regexprep(line, ' +', ' '), '^ ', ''); % Remove multiple spaces, and then eventually remove space if it there is one as first character.
      if strcmp(line(1:2), 'xs')
        xfound = 1; pos_sources(1, 1) = str2num(regexprep(regexprep(line(3:end), ' *#.*', ''), ' * = * *', '')); % Remove comments (everything after a '#'), remove the equals sign and spaces around it, and cast it as source position.
      end
      if strcmp(line(1:2), 'zs')
        zfound = 1; pos_sources(1, 2) = str2num(regexprep(regexprep(line(3:end), ' *#.*', ''), ' * = * *', ''));
      end
    end
    if (xfound && zfound)
      break
    end
  end
  fclose('all');
end