% Author:        Léo Martire.
% Description:   Draws nice brackets onton some plot.
% Notes:         To remove the brackets:
%                  1) Run:
%                     kk=findall(gca,'type','line')
%                  2) Find the IDs of the bracket lines, let's call them IDzzz.
%                  3) Run:
%                     delete(kk(IDzzz));
%                To change the color of the brackets, run:
%                  set(linehandles, 'color', 'k')
%                or anything else.
%
% Usage:
%   linehandles = drawBrackets(axxx, X)
%   linehandles = drawBrackets(axxx, X, topbottom)
%   linehandles = drawBrackets(axxx, X, topbottom, opening0closing1)
%   linehandles = drawBrackets(axxx, X, topbottom, opening0closing1, liplength)
% with:
%   axxx             some already existing plot axis,
%   X                the X position(s) of the brackets (either a length 1 or 2 vector)
%   topbottom        (optionnal) a vector of Y positions of vertical span
%                    of brackets (or 'auto'),
%   opening0closing1 (optionnal, used only if numel(X)==1) code for
%                    opening, closing, or simple line (0 for an opening
%                    bracket, 1 for a closing bracket, anything else for a
%                    simple line),
%   liplength        (optionnal, only used if opening0closing1 is neither
%                    0 nor 1) length of the lips of the brackets,
% yields:
%   linehandles      an array of handles to the created brackets.

function linehandles = drawBrackets(axxx, X, topbottom, opening0closing1, liplength)
  if(not(exist('topbottom')) || strcmp(topbottom, 'auto'))
    topbottom = get(axxx, 'ylim') * 0.5;
  end
  if(not(exist('liplength')))
    liplength_given = 0;
  else
    liplength_given = 1;
  end
  
  % parameters
  oneBracketPercentageForLip = 5; % lip in the 1-bracket case will be #% of xlims
  twoBracketPercentageForLip = 5; % lip in the 2-bracket case will be #% of bracket span
  
%   topbottom
  
  linehandles=[];
  switch(numel(X))
    case(1)
      if(not(exist('opening0closing1')) || not(ismember(opening0closing1, [0, 1])))
        % draws a simple line
        linehandles = [linehandles, line(X*[1, 1], topbottom, 'HandleVisibility', 'off')];
      else
        % draw a closing (opening0closing1=0) or opening bracket (opening0closing1=1)
        if(liplength_given)
          lip = liplength;
        else
          lip = diff(get(axxx,'xlim')) * oneBracketPercentageForLip;
        end
        switch(opening0closing1)
          case(0)
            linehandles=[linehandles, drawOneBracket(X, topbottom, lip,  1)]; % lips towards right (opening)
          case(1)
            linehandles=[linehandles, drawOneBracket(X, topbottom, lip, -1)]; % lips towards left (closing)
          otherwise
            error('kek');
        end
      end
    case(2)
      if(liplength_given)
        disp(['[',mfilename,'] Lip length given in a 2-bracket case, overriding default automatic value.']);
        lip = liplength;
      else
        lip = diff(X) * twoBracketPercentageForLip/100;
      end
      linehandles=[linehandles, drawOneBracket(X(1), topbottom, lip,  1)]; % lips towards right (opening)
      linehandles=[linehandles, drawOneBracket(X(2), topbottom, lip, -1)]; % lips towards left (closing)
    otherwise
      error('kek');
  end
end

function hh = drawOneBracket(x, Y, lip, openingClosingSign)
  hh = [];
  hh = [hh, line(x*[1, 1], Y, 'HandleVisibility', 'off')];
  hh = [hh, line(x+[0, lip]*openingClosingSign, Y(2)*[1, 1], 'HandleVisibility', 'off')];
  hh = [hh, line(x+[0, lip]*openingClosingSign, Y(1)*[1, 1], 'HandleVisibility', 'off')];
end