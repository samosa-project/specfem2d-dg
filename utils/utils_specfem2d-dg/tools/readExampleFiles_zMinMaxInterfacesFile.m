% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   [ymin, ymax] = readExampleFiles_zMinMaxInterfacesFile(interfaces_file)
% with:
%   TODO.
% yields:
%   TODO.

function [ymin, ymax] = readExampleFiles_zMinMaxInterfacesFile(interfaces_file, verbose)
  error('now, prefer the use of readExampleFiles_meshfem_mesh');
  if(not(exist('verbose')))
    verbose=0;
  end
  
  if(not(exist(interfaces_file)==2))
    if(verbose)
      disp(['[',mfilename,', ERROR] ''',interfaces_file,''' is not a path to a file. Setting outputs to NaNs.']);
    end
    ymin = nan;
    ymax = nan;
  else
  %   interfaces_file=[OFd,'input_interfaces'];
    grep_remove_comments=['grep -v "^#.*"'];
    numberregexp='\-?[0-9]+\.?[0-9]*(e|d)?\+?[0-9]*';
    grep_find_interface_lines=['grep -P "',[numberregexp, ' ',numberregexp],'"'];

    [~,y]=system(['cat ',interfaces_file,' | ',grep_remove_comments,' | ',grep_find_interface_lines,' | head -1']);
    y=str2num(y);
    ymin=y(end);

    [~,y]=system(['cat ',interfaces_file,' | ',grep_remove_comments,' | ',grep_find_interface_lines,' | tail -1']);
    y=str2num(y);
    ymax=y(end);
  end
end

