clear all;
% close all;
clc;

usegll = 1;
Nexact = 2000+1;
% N2test = [1:1:9,logspace(1,2,3)];
N2test = [1, 100];

testCases = {'inviscid', 'kappa', 'mu'};

x_exact = linspace(0,1,Nexact);
y_exact = linspace(0,1,Nexact);
[Xe, Ye] = meshgrid(x_exact, y_exact);

for tci = 1:numel(testCases)
  testCase = testCases{tci};
  Ze = func2test(testCase, Xe, Ye);
  % figure(); pcolor(Xe,Ye,Ze); shading flat; title(testCase); pause;
  for ni = 1:numel(N2test)
    n = N2test(ni);
    x = linspace(0,1,n+1);
    y = linspace(0,1,n+1);

    if(usegll)
      xgll = sort(lglnodes(5)+1)*0.5;
      dx = mean(diff(x)); x = x(1:end-1); x = x+xgll * dx; x= unique(x(:));
      dy = mean(diff(y)); y = y(1:end-1); y = y+xgll * dy; y= unique(y(:));
    end

    [X, Y] = meshgrid(x, y);
    Z.rho = [];
    Z.vel = [];
    Z.pre = func2test(testCase, X, Y); Z.pre = Z.pre(:);
    [Xi, Yi, Vi] = interpDumps(X(:), Y(:), Z, Nexact, Nexact);
  %   figure(); surf(Xi,Yi,Vi.pre); shading flat;
    zerr = Ze-Vi.pre;
  %   figure(); surf(Xi,Yi,zerr); shading flat;
    L2NormOfError{tci}(ni) = (trapz(x_exact,trapz(y_exact,zerr.^2)))^0.5;
  end
end

figure();
for tci = 1:numel(testCases)
  [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, colourPlot, ~, ~] = MMS_constants(testCases{tci}); % Get plot parameters as function of test case.
  switch(testCases{tci})
    case 'inviscid'
      LW = 5;
    otherwise
      LW = 3;
  end
%   LW, testCases{tci}
  loglog(N2test, L2NormOfError{tci}, 'color', colourPlot, 'linewidth', LW, 'displayname', ['', testCases{tci}, ', theoretical error']); hold on;
end
legend()

function Z = func2test(testCase, X, Y)
%   Z = -sin(2*pi*X/1);
  RHO0 = 1.4;
  V0_x = 10;
  C = 340;
  GAM = 1.401;
  
  P0 = C^2*RHO0/GAM;
  E0 = P0/(GAM-1) + 0.5*RHO0*V0_x^2;
  [drho_th, dvx_th, dvz_th, dE_th, dp_th] = MMS_analytic(testCase, X, Y, RHO0, V0_x, E0, P0, GAM);
  Z = dp_th;
end