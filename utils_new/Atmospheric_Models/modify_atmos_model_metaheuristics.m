% Author:        Léo Martire.
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         Modules used by 'modify_atmos_model.m'.

% clear all;
% close all;
% clc;
set(0, 'DefaultLineLineWidth', 1.5); % Default at 0.5.
set(0, 'DefaultLineMarkerSize', 6); % Default at 6.
%set(0, 'defaultTextFontSize', 20);
set(0, 'defaultAxesFontSize', 12); % Default at 10.
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

format longG

global D Z G % Not very good pratice, to remove.
rng default % For reproducibility, remove before flight.

if(not(exist('callFromOutside')))
  % If this script is not called from outside script, choose an algorithm.
  % If it is, the variables here should have been initialised.
  % algo = 'sa';
  % algo = 'ga';
  algo = 'ps';
  ftol = 1e-12;
  maxit = Inf;
end

useTestData=0;
if(not(exist('Z') & exist('P') & exist('RHO')) | useTestData)
  % If data does not exist, or if the user want test data, use test data.

  % Test data.
  Z=linspace(0,60000,111)';
  D = differentiation_matrix(Z, 0);
  div=0.1;
  P=101325*exp(-Z/10000);
  G=0*P+9.81;
  RHO=-D*P./G;
  P=P.*(1+div*sin(20*Z/Z(end)));
  RHO=RHO.*(1+div*sin(90*Z/Z(end)));
  % HR0=hydrostatic_ratio(D, P, RHO, G); plot(HR0(1:end), Z(1:end), 'DisplayName', 'original');
else
  % Else, we anyway need D.
  D = differentiation_matrix(Z, 0);
end

n=length(P);

% Weights for objective function.
% w = logspace(0,-5,n)';
w = logspace(2,-2,n)'/100; %100% to 0.1%
% w = P/P(1);
% w = RHO/RHO(1);

% Objective function only depending on X.
fun = @(X) (objective_function(X, w));

% Initial state (only used by SA).
x0 = [P;RHO];

% Bounds.
% lb = [0*ones(n,1);0*ones(n,1);0*ones(n,1)]+1e-5;
% ub = 1.5*[max(P)*ones(n,1);max(RHO)*ones(n,1);max(G)*ones(n,1)];
dpc=10/100;
lb = (1-dpc)*x0;
ub = (1+dpc)*x0;
% ub = lb+max((1+dpc)*x0-lb,0.0021); % Only like this to make Matlab's SA work (because it approximates derivatives and complains if ub-lb is too small.
% ub = lb+max((1+dpc)*x0-lb,0.0002); % Only like this to make Matlab's SA work (because it approximates derivatives and complains if ub-lb is too small.

% In case we deal with negative values, flip upper/lower bounds.
oub=ub; olb=lb; nub=lb(oub<olb); nlb=ub(oub<olb); ub(oub<olb)=nub; lb(oub<olb)=nlb;

tic;
if(strcmp(algo,'sa'))
  disp(["Launching optimisation by simulated annealing."]);
  % options = saoptimset('PlotFcns',{@saplotbestx,@saplotbestf,@saplotx,@saplotf}, 'MaxIter', 1000);%, 'AnnealingFcn', @test_evo_newpt);
  options = saoptimset('PlotFcns',{@saplotbestf} ...
                       , 'MaxIter', maxit ...
                       , 'InitialTemperature', 100 ...
                       , 'TolFun', ftol ...
                       , 'PlotInterval', 1000 ...
                       );
                       %, 'AnnealingFcn', @test_evo_newpt
                       %, 'MaxIter', 10000
                       %, 'AnnealingFcn', @test_evo_newpt
  [x, fval, exitflag, output] = simulannealbnd(fun,x0,lb,ub,options);
elseif(strcmp(algo,'ga'))
  disp(["Launching optimisation by genetic algorithm."]);
  options = gaoptimset('PlotFcns', {@gaplotbestf} ...
                       , 'PopulationSize', 200 ...
                       , 'Generations', maxit ...
                       , 'TolFun', ftol ...
                       , 'PlotInterval', 100 ...
                       );
                       %, 'UseParallel', true
  [x,fval,exitflag] = ga(fun,numel(x0),[],[],[],[],lb,ub,[], options);
elseif(strcmp(algo,'ps'))
  disp(strcat(['Launching particle swarm optimisation (ftol=',num2str(ftol),', maxit=',num2str(maxit),').']));
  options = optimoptions('particleswarm' ...
                         , 'PlotFcns', {@pswplotbestf} ...
                         , 'SwarmSize', 200 ...
                         , 'InitialSwarmSpan', x0 ...
                         , 'MaxIterations', maxit ... %200*numel(x0) ...
                         , 'TolFun', ftol ...
                         );
  [x,fval,exitflag] = particleswarm(fun,numel(x0),lb,ub,options);
end
exectime=toc;
h=floor(exectime/3600); m=floor((exectime-3600*h)/60); s=floor(exectime-3600*h-60*m);
disp(['Execution time: ',num2str(h), 'h ',num2str(m),'m ',num2str(s),'s.']);

if(exitflag~=-1)
  disp(['Exit flag: ', num2str(exitflag), '.']);
  x=reshape(x,length(x),1);
  hydrostatic_ratio = @(D, P, RHO, G) (D * P) ./ (-RHO .* G);

  nP=x(1:n);
  nRHO=x(n+1:2*n);
  % nG=x(2*n+1:3*n);
  nG=G;

  % TRICK
%   FINALX=x;
%   originalBETTER=(abs(hydrostatic_ratio(D, P, RHO, G)-1)<abs(hydrostatic_ratio(D, nP, nRHO, nG)-1));
%   originalBETTER=[originalBETTER;originalBETTER];
%   FINALX(originalBETTER)=x0(originalBETTER);
%   FnP=FINALX(1:n);
%   FnRHO=FINALX(n+1:2*n);
%   FnG=G;

  figure();
  HR0=hydrostatic_ratio(D, P, RHO, G);
  HR=hydrostatic_ratio(D, nP, nRHO, nG);
  plot(HR0(1:end), Z(1:end), 'DisplayName', 'original'); hold on;
  plot(HR(1:end), Z(1:end), 'DisplayName', 'after optimisation');
  plot(ones(size(Z)), Z, 'k:', 'DisplayName', 'one');
  title('hydrostatic ratio');
  legend('Location', 'best');
%   FHR=hydrostatic_ratio(D, FnP, FnRHO, FnG);
%   plot(FHR(1:end), Z(1:end), 'DisplayName', 'after trick');

  figure();
  HG0=abs(D*P+RHO.*G)./(RHO.*G);
  HG=abs(D*nP+nRHO.*nG)./(nRHO.*nG);
  % FHG=abs(D*FnP+FnRHO.*FnG)./(FnRHO.*FnG);
  semilogx(HG0(1:end), Z(1:end), 'DisplayName', 'original'); hold on;
  semilogx(HG(1:end), Z(1:end), 'DisplayName', 'after optimisation');
  % plot(FHG(1:end), Z(1:end), 'DisplayName', 'after trick');
  title('relative gap');
  legend('Location', 'best');
  
  chngP=100*abs(P-nP)./P;
  chngRHO=100*abs(RHO-nRHO)./RHO;
  disp(['changes (%):   P: min: ', sprintf('%1.1e', min(chngP)), ', avg: ', sprintf('%1.1e', mean(chngP)), ', max: ', sprintf('%1.1e', max(chngP))]);
  disp(['             RHO: min: ', sprintf('%1.1e', min(chngRHO)), ', avg: ', sprintf('%1.1e', mean(chngRHO)), ', max: ', sprintf('%1.1e', max(chngRHO))]);
  figure();
  plot(chngP, Z, chngRHO, Z);
  title('Percent change from model');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Objective function.
function y = objective_function(PRHO, w)
  global D Z G
%   tn=length(PRHO);
%   PRHO=reshape(PRHO,tn,1);
%   P = PRHO(1:tn/2);
%   RHO = PRHO(tn/2+1:end);
%   epsilon=w.*metric(PRHO);
%   y=sqrt(trapz(Z, epsilon.^2)/(Z(end)-Z(1)));
%   y=trapz(Z, abs(epsilon))/(Z(end)-Z(1));
  y=trapz(Z, abs(w.*metric(PRHO)));
end
% Metric.
function epsilon = metric(PRHO)
  global D G
  tn=length(PRHO);
  PRHO=reshape(PRHO,tn,1);
%   P = PRHO(1:tn/2);
%   RHO = PRHO(tn/2+1:end);
%   epsilon=(D * P) ./ (-RHO .* G)-1;
  epsilon=(D * PRHO(1:tn/2)) ./ (-PRHO(tn/2+1:end) .* G)-1;
end

% Objective function.
function y = objective_function_wind(NW, w)
  global D Z G
%   tn=length(PRHO);
%   PRHO=reshape(PRHO,tn,1);
%   P = PRHO(1:tn/2);
%   RHO = PRHO(tn/2+1:end);
%   epsilon=w.*metric(PRHO);
%   y=sqrt(trapz(Z, epsilon.^2)/(Z(end)-Z(1)));
%   y=trapz(Z, abs(epsilon))/(Z(end)-Z(1));
  y=trapz(Z, abs(metric_wind(NW)));
end
% Metric.
function epsilon = metric_wind(NW)
  global D G
  tn=length(NW);
  NW=reshape(NW,tn,1);
%   P = PRHO(1:tn/2);
%   RHO = PRHO(tn/2+1:end);
%   epsilon=(D * P) ./ (-RHO .* G)-1;
  epsilon= ((NW(1:tn/2) ./ (D * NW(tn/2+1:end))) .^ 2.0) <1;
end