% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   [var] = readExampleFiles_extractParam(path2file, varName, varType)
% with:
%   TODO.
% yields:
%   TODO.

function [var] = readExampleFiles_extractParam(path2file, varName, varType, verbose)
  if(not(exist('verbose')))
    verbose=0;
  end

%   path2parfile='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/validation_lns_gravito/OUTPUT_FILES/input_parfile';
%   varName='NPROC';
%   varType='int';
  
  if(not(exist(path2file,'file')))
%     if(verbose)
      %error(['[',mfilename,', ERROR] File ''',path2file,''' does not exist. Run ''help ',mfilename,'''.']);
      disp(['[',mfilename,', ERROR] ''',path2file,''' is not a path to a file. Setting outputs to NaNs.']);
%     end
    var = nan;
    
  else
    command=['grep -r "^ *',varName,' *=" ',path2file];
    [~,res] = system(command);
    if(isempty(res))
      disp(['[',mfilename,', ERROR] Found nothing.']);
      switch(varType)
        case{'int', 'float'}
          disp([blanks(length(mfilename)+9),' Setting requested int/float to NaN.']); var = nan;
        case{'bool', 'boolean'}
          disp([blanks(length(mfilename)+9),' Setting requested boolean to 0.']); var = 0;
        case{'string'}
          disp([blanks(length(mfilename)+9),' Setting requested string to ''''.']); var = '';
        otherwise
          error('kek');
      end
    end
    res = regexprep(res,'#[^\n]*',''); % remove comments
  %   res = regexprep(res,'[\n\r]+',' ') % remove endline
    res = splitlines(res); % split occurences
  %   grp_remove_comments=['grep -v "#.*"'];

    for l=1:numel(res)
      resl=res{l};
      if(not(isempty(resl)))
        command=['echo "',resl,'" | grep -oP "=.*"'];
        if(ismac)
          command = regexprep(command,'-oP','-oe');
          command = regexprep(command,'\+','\\+');
        end
        [~, resl] = system(command);
        command=['echo "',resl,'" | grep -oP "[^=]"'];
        if(ismac)
          command = regexprep(command,'-oP','-oe');
          command = regexprep(command,'\+','\\+');
        end
        [~, resl] = system(command);

        switch(varType)
          case{'int', 'float'}
            var(l) = str2num(regexprep(resl,'[\n\r]+','')); % remove endline

          case{'bool', 'boolean'}
            resl=regexprep(resl,'[\n\r]+',''); % remove endlines
            resl=regexprep(resl,' +',''); % remove spaces
            resl=regexprep(resl,'\.+',''); % remove dots
            if(strcmp(resl,'true'))
              var=1;
            else
              var=0;
            end

          case{'string'}
            resl=regexprep(resl,'[\n\r]+',''); % remove endlines
            resl=regexprep(resl,' +',''); % remove spaces
            var = resl;

          otherwise
            error('kek');
        end

      end
    end
  end
end

