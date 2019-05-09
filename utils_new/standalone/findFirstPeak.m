% Author:        Léo Martire.
% Description:   Tries to find the first peak in some data by a very rough
%                dichotomy algorithm.
% Notes:         TODO.
%
% Usage:
%   [found_x] = findFirstPeak(X, Y, start_rate, maxiter, verbose, processPlot)
% with:
%   TODO
%   start_rate a float in [0, 1] representing the minimum positive relative
%              amplitude the peak must have to be found,
%   TODO
% yields:
%   found_x    the abscissa (from X) where the peak was found.

% X=linspace(0,4*pi,1000); Y=0.5*abs(0.1*treat_t+sin(treat_t)); % test data

function [found_x] = findFirstPeak(X, Y, start_rate, maxiter, verbose, processPlot)
  if(not(exist('start_rate')))
    start_rate = 0.5;
  end
  if(not(exist('maxiter')))
    maxiter = 10;
  end
  if(not(exist('verbose')))
    verbose = 0;
  end
  if(not(exist('processPlot')))
    processPlot = 0;
  end
  
  if(processPlot)
    figure();
    plot(X, Y); hold on;
%     xlim([53,58]);
  end

  mask = ones(size(X)); % mask will do nothing
  rate = start_rate; % starting point
  rate_upp = 1; % upper bound
  rate_low = 0.1; % lower bound
  iter = 0; 
  
  while(iter < maxiter)
    locsel = (Y > rate*max(Y));
    locsel = (locsel & mask);
    
    if(verbose)
      disp(['Iteration ', num2str(iter)]);
    end
    
    if(not(any(locsel)))
      if(verbose)
        disp(['  No peaks found with amplitude over ', num2str(rate), ' *amplitude']);
      end
      
      rate_old = rate;
      rate = (rate_low+rate)*0.5; % update rate upwards
      rate_upp = rate_old;  % updat lower bound
      
      if(verbose)
        disp(['    Updating rate downwards: ', num2str(rate), ' (lower bound ',num2str(rate_low),')']);
      end
    else
      if(numel(find(locsel==1))==1)
        if(verbose)
          disp(['    One unique point found: t=', num2str(found_x), '. Stopping algorithm.']);
        end
        
        found_x = X(locsel);
        break;
      else
        Npeaksfound=numel(find(diff(locsel)<0)); % count the drops in diff to count the number of windows
        
        if(verbose)
          disp(['  Found ',num2str(Npeaksfound), ' peak(s) with amplitude over ',num2str(rate), ' *amplitude']);
          disp(['    Selecting first one of those, updating mask.']);
        end
        
        bounds_ids=find(abs(diff(locsel))>0,2,'first');
        mask = 0*mask; mask(bounds_ids(1):bounds_ids(2))=1;
        
        found_x = mean(X(mask==1));
        
        if(verbose)
          disp(['    Midpoint of interval previously found: t=',num2str(found_x),' (interval is [',num2str(min(X(mask==1))),', ',num2str(min(X(mask==1))),'])']);
        end

        rate_old=rate;
        rate = (rate + rate_upp)*0.5; % update rate upwards
        rate_low = rate_old; % updat lower bound
      end
    end
    
    if(processPlot)
      plot(X, locsel);
    end
    
    iter = iter+1;
  end
  legend();
end

