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

function [ymin, ymax] = readExampleFiles_zMinMaxInterfacesFile(interfaces_file)
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

