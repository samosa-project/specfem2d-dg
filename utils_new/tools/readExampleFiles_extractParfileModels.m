% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [models] = readExampleFiles_extractParfileModels(parfile)
  tag_before='^nbmodels';
  tag_after='^TOMOGRAPHY';
  hugenumberoflines=1e5;
  grp_find_after=['grep "',tag_before,'" "',parfile,'" -A',sprintf('%d',hugenumberoflines)];
  grp_find_before=['grep "',tag_after,'" -B',sprintf('%d',hugenumberoflines)];
  grp_remove_commentlines=['grep -v "^#.*"'];
  grp_remove_tags=['grep -v "',tag_before,'" | grep -v "',tag_after,'"'];
  [stat,res]=system([grp_find_after,' | ',grp_find_before,' | ',grp_remove_commentlines,' | ',grp_remove_tags]);
  if(isempty(res))
    error(['[',mfilename,', ERROR] grep returned nothing. Make sure you gave a path to a valid parfile.']);
  end
  returnsfoundIDs=regexp(res,'\n');
  models=[];
  strt=1;
  j=1;
  format shortg
  for i=1:numel(returnsfoundIDs)
    % Parse line.
    txt=res(strt:returnsfoundIDs(i)-1);
    % Save line if not empty.
    if(not(isempty(txt)))
      % Find comments.
      comments=regexp(txt,'(#.*)');
      % Remove comments.
      if(not(isempty(comments)))
        txt=txt(1:comments-1);
      end
      % Change all "d" to "e", and multiple spaces by single spaces.
      txt=regexprep(strrep(txt,'d','e'),' +',' ');
      % If last character is space, remove it (can happen if there were comments).
      if(txt(end)==' ')
        txt=txt(1:end-1);
      end
      % Convert to array.
      models(j,:)=str2double(split(txt,' '));
      j=j+1;
    end
    strt=returnsfoundIDs(i)+1;
  end
end

