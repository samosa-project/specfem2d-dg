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
%   f_0      a main frequency,
%   np       a number of points per wavelength,
% yields:
%   interf_s the depths of the interfaces (top one being z=0),
%   nelts_s  the number of elements per layer,
%   IDparf   the corresponding model IDs for each layer in the parfile.

function [interf_s, nelts_s, IDparf] = model_viscoelastic_getinterfaces(f0, np)
  if(nargin~=2)
    help(mfilename)
    error(['[',mfilename,', ERROR] Not enough input arguments.']);
  end

  parfile=input(['[',mfilename,'] Path to parfile > '],'s');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Parsing parfile to find     %
  % models.                     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  addpath([regexprep(mfilename('fullpath'),mfilename,''),'tools']);
  models = readExampleFiles_extractParfileModels(parfile);

  % Crop unnecessary quantities.
  interestingquantitiesIDs=[1,3,4,5,8,9];
  %disp(['[',mfilename,'] Models as loaded from parfile (number, 1, rho, vp, vs, 0, 0, QKappa, Qmu, 0, 0, 0, 0, 0, 0): ']);
  disp(['[',mfilename,'] Models as loaded from parfile (number, rho, vp, vs, QKappa, Qmu): ']);
  disp(models(:,interestingquantitiesIDs));

  % Ask user for which models to use.
  IDparf=[];
  while (not(length(IDparf)>=1 && min(IDparf)>=1 && max(IDparf)<=size(models,1)))
    IDparf=input(['[',mfilename,'] Model numbers to consider as viscoelastic (Matlab format, >=1 and <= ',num2str(size(models,1)),')? > ']);
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
%   disp(['[',mfilename,'] Chosen thicknesses: ']);
%   disp(thicks);

  % Generate interfaces.
  interf_s=[0;0-cumsum(thicks)];

  % Generate optimal number of elements with formula f0*np*thickness/vp.
  recompute_in=1;
  recompute=1;
  while(recompute)
    if(numel(recompute_in)==1 && recompute_in==0)
%       disp(['[',mfilename,']   > ==0: ']);
      recompute=recompute_in;
      override=0;
      break;
    elseif(numel(recompute_in)==1 && recompute_in~=0)
%       disp(['[',mfilename,']   > mult. factor: ']);
      recompute=recompute_in;
      override=0;
    elseif(numel(recompute_in)==numel(nelts_s))
%       disp(['[',mfilename,']   > override: ']);
      recompute=1;
      override=1;
      nelts_s=recompute_in;
    else
      disp(['[',mfilename,'] > Wrong input, going back to default.']);
      recompute=1;
      override=0;
    end
    if(not(override))
      nelts_s = ceil(recompute * f0 * np * thicks ./ models(:,3));
    end
    nelts_s=ceil(reshape(nelts_s,[numel(nelts_s),1]));
    disp(['[',mfilename,'] Computed [thicks, nelts_s, dz]: ']);
    disp([thicks, nelts_s, thicks./nelts_s]);
    recompute_in=input(['[',mfilename,'] Recompute (0 for no, mutiplying factor for yes, or overriding new array for yes)? > ']);
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