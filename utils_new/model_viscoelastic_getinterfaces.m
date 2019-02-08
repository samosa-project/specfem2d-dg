% Author:        Léo Martire.
% Description:   Parses a parfile, deduce optimal DZ from user-provided
%                thicknesses, a main frequency, and a number of points per
%                wavelength.
% Notes:         Needs to be on Unix to access GNU grep (version 3.1 works
%                fine) and parse the parfile.
%
% Usage:
%   [interf_s, nelts_s, IDparf] = model_viscoelastic_getinterfaces(f0, np)
% with:
%   f_0 a main frequency,
%   np a number of points per wavelength,
% yields:
%   interf_s the depths of the interfaces (top one being z=0),
%   nelts_s the number of elements per layer,
%   IDparf the corresponding model IDs for each layer in the parfile.

function [interf_s, nelts_s, IDparf] = model_viscoelastic_getinterfaces(f0, np)
  if(nargin~=2)
    error(['[',mfilename,', ERROR] Not enough input arguments. Needs ''f0, np''.']);
  end

  parfile=input(['[',mfilename,'] path to parfile > '],'s');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Parsing parfile to find     %
  % models.                     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tag_before='^nbmodels';
  tag_after='^TOMOGRAPHY';
  hugenumberoflines=1e5;
  grp_find_after=['grep "',tag_before,'" "',parfile,'" -A',sprintf('%d',hugenumberoflines)];
  grp_find_before=['grep "',tag_after,'" -B',sprintf('%d',hugenumberoflines)];
  grp_remove_commentlines=['grep -v "^#.*"'];
  grp_remove_tags=['grep -v "',tag_before,'" | grep -v "',tag_after,'"'];
  [stat,res]=system([grp_find_after,' | ',grp_find_before,' | ',grp_remove_commentlines,' | ',grp_remove_tags]);
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

  % Crop unnecessary quantities.
  interestingquantitiesIDs=[1,3,4,5,8,9];
  %disp(['[',mfilename,'] Models as loaded from parfile (number, 1, rho, vp, vs, 0, 0, QKappa, Qmu, 0, 0, 0, 0, 0, 0): ']);
  disp(['[',mfilename,'] Models as loaded from parfile (number, rho, vp, vs, QKappa, Qmu): ']);
  disp(models(:,interestingquantitiesIDs));

  % Ask user for which models to use.
  IDparf=[];
  while (not(length(IDparf)>=1 && min(IDparf)>1 && max(IDparf)<=size(models,1)))
    IDparf=input(['[',mfilename,'] Model numbers to consider as viscoelastic (Matlab format, between 1 and <= ',num2str(size(models,1)),')? > ']);
  end

  % Crop unnecessary quantities.
  models=models(IDparf,interestingquantitiesIDs);

  Nlayerz_solid=size(models,1);

  % Ask user for thicknesses.
  thicks=[];
  for i=1:Nlayerz_solid
    % Ask thicknesses.
    thick=[];
    while (not(length(thick)==1 && thick>0))
      thick=input(['[',mfilename,'] Thickness of layer with model ',num2str(i),' (',num2str(models(i,:)),') [m]? > ']);
    end
    thicks(i)=thick;
  end
  thicks=thicks';
  disp(['[',mfilename,'] Chosen thicknesses: ']);
  disp(thicks);

  % Generate interfaces.
  interf_s=[0;0-cumsum(thicks)];

  % Generate optimal number of elements with formula f0*np*thickness/vp.
  recompute=1;
  while(recompute)
    nelts_s = ceil(recompute * f0 * np * thicks ./ models(:,3));
    disp(['[',mfilename,'] Computed number of elements: ']);
    disp(nelts_s);
    recompute=input(['[',mfilename,'] Recompute (0 for no, mutiplying factor for yes)? > ']);
  end
  
  % Make bottom-most one first.
  interf_s=flip(interf_s);
  nelts_s=flip(nelts_s);
  IDparf=flip(IDparf);
  % Safeguard.
  if(not(all(sort(interf_s)==interf_s)))
    error(['[',mfilename,', ERROR] Something is wrong with interfaces storing. Check ''',mfilename,'''.']);
  end
end