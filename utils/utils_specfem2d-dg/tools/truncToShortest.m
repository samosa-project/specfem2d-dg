% Author:        Léo Martire.
% Description:   TODO.
% Notes:         TODO.
%
% Usage:
%   [Ztime, Zamp] = truncToShortest(raw_t, raw_s)
% with:
%   TODO.
% yields:
%   TODO.

function [Ztime, Zamp] = truncToShortest(raw_t, raw_s, verbose)
  if(not(all(size(raw_t)==size(raw_s))))
    error('must have same sizes');
  end
  if(not(exist('verbose', 'var')))
    verbose = 1;
  end
  Nsamples = size(raw_t, 2);
  difftimes = raw_t(:,2:end)-raw_t(:,1:end-1);
  relevantdifftimes = (difftimes>0);
  
%   minimum_last_relevant = +Inf;
%   for i = 1:size(raw_t, 1)
%     locminim = find(relevantdifftimes(i,:)>0, 1, 'last'); % Find last index having difftime > 0.
%     if(locminim<minimum_last_relevant)
%       minimum_last_relevant = locminim; % Save smallest index having difftime > 0 across all data.
%     end
%   end
% %   minimum_last_relevant
%   if(minimum_last_relevant<Nsamples)
%     Ztime = raw_t(:,1:minimum_last_relevant);
%     Zamp  = raw_s(:,1:minimum_last_relevant);
%     disp(['[',mfilename,'] Data in data pack was not of same length. Truncated to shortest, t is now [',num2str(min(min(Ztime))),', ',num2str(max(max(Ztime))),'] (was [',num2str(min(min(raw_t))),', ',num2str(max(max(raw_t))),']).']);
%   else
%     Ztime = raw_t;
%     Zamp  = raw_s;
%   end
  
  % build selection by consecutive AND vector operations.
  select = relevantdifftimes(1,:);
  for i = 2:size(relevantdifftimes,1)
    select = and(select, relevantdifftimes(i,:));
  end
  if(not(any(select)))
    error(['[',mfilename,', ERROR] Nothing in selection, stopping script.']);
  end
  % Here, select is true only where difftimes is >0 for ALL lines.
  Ztime = raw_t(:, select);
  Zamp  = raw_s(:, select);
  if(verbose)
    disp(['[',mfilename,'] Truncated to shortest reliable time vector. t is now [',num2str(min(min(Ztime))),', ',num2str(max(max(Ztime))),'] (was [',num2str(min(min(raw_t))),', ',num2str(max(max(raw_t))),']).']);
  end
end

