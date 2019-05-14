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

function niceFormatForYTickLabels(fighandle, maxprecision)
  childs=fighandle.Children;
  
  if(not(exist('maxprecision')))
    maxprecision=1;
  end
  
  for i=1:numel(childs)
    child = childs(i);
    if(strcmp(child.Type,'axes'))
      yticks = get(child,'ytick');
      
%       txt = sprintf(['%.',num2str(maxprecision),'e '],[0.2e13,1.123456789e13,0.4e13,-1e12,0,-3e-15]); % put into nice parseable format
      txt = sprintf(['%.',num2str(maxprecision),'e '],yticks); % put into nice parseable format
      txt = regexprep(txt,'e+','\{\\times\}10^{'); % build latex exponent
      txt = regexprep(txt,'+',''); % remove + (which should only appear at exponent)
      txt = regexprep(txt,' ','}$ $'); % add closing }$ and opening $
      txt = ['$',txt]; % add opening latex
      
      while(not(isempty(regexp(txt,'0\{\', 'once'))))
        txt = regexprep(txt,'0\{\','{\'); % remove as many unnecessary zeros as needed
      end
      txt = regexprep(txt,'\.\{\','{\'); % remove unnecessary dots
      
      while(not(isempty(regexp(txt,'\^\{0', 'once'))))
        txt = regexprep(txt,'\^\{0','^{'); % remove as many unnecessary '^{0' as needed
      end
      while(not(isempty(regexp(txt,'\^\{-0', 'once'))))
        txt = regexprep(txt,'\^\{-0','^{-'); % idem with '^{-0'
      end
      
      txt = regexprep(txt,'\{\\times\}10\^\{\}',''); % replace empty 10^{} by nothing
      
      txt = txt(1:end-2); % remove ending $ and ending space
      spl = split(txt);
%       spl = split(string(txt));
      
      set(child, 'yticklabels', spl);
    end
  end
end