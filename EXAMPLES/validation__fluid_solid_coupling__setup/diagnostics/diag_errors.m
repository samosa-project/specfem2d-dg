function [fig] = diag_errors(sf_ortho_err, sf_slant_err, fs_ortho_err, fs_slant_err)
  setupLocal;

  fig = figure('units','pixels','position', [0,0,1300,700]);
  ax = tight_subplot(1,1,0,[.14,.04],[.1,.01]);
  l0 = 370/100; % smallest wavelength
  x = l0./(200./nvals);
  col = get(0,'defaultAxesColorOrder');
  ls = ':';
  stric = [sprintf('%.0f',ic_deg),'^\circ'];
  tsfpp = '$T^\mathrm{SF}_\mathrm{PP}$';
  rsfpp = '$R^\mathrm{S}_\mathrm{PP}$';
  rsfps = '$R^\mathrm{S}_\mathrm{PS}$';
  tfspp = '$T^\mathrm{FS}_\mathrm{PP}$';
  tfsps = '$T^\mathrm{FS}_\mathrm{PS}$';
  rfspp = '$R^\mathrm{F}_\mathrm{P}$';
  % dnam = {['S2F, ',tsfpp,', $\theta_\mathrm{i}=0$'],         ['S2F, ',rsfpp,', $\theta_\mathrm{i}=0$'],         ['S2F, ',rsfps,', $\theta_\mathrm{i}=0$'], ...
  %         ['S2F, ',tsfpp,', $\theta_\mathrm{i}=',stric,'$'], ['S2F, ',rsfpp,', $\theta_\mathrm{i}=',stric,'$'], ['S2F, ',rsfps,', $\theta_\mathrm{i}=',stric,'$'];
  %         ['F2S, ',tfspp,', $\theta_\mathrm{i}=0$'],         ['F2S, ',tfsps,', $\theta_\mathrm{i}=0$'],         ['F2S, ',rfspp,', $\theta_\mathrm{i}=0$'], ...
  %         ['F2S, ',tfspp,', $\theta_\mathrm{i}=',stric,'$'], ['F2S, ',tfsps,', $\theta_\mathrm{i}=',stric,'$'], ['F2S, ',rfspp,', $\theta_\mathrm{i}=',stric,'$']};
  dnam = {['',tsfpp,', $\theta_\mathrm{i}=0$'],         ['',rsfpp,', $\theta_\mathrm{i}=0$'],         ['',rsfps,', $\theta_\mathrm{i}=0$'], ...
          ['',tsfpp,', $\theta_\mathrm{i}=',stric,'$'], ['',rsfpp,', $\theta_\mathrm{i}=',stric,'$'], ['',rsfps,', $\theta_\mathrm{i}=',stric,'$'];
          ['',tfspp,', $\theta_\mathrm{i}=0$'],         ['',tfsps,', $\theta_\mathrm{i}=0$'],         ['',rfspp,', $\theta_\mathrm{i}=0$'], ...
          ['',tfspp,', $\theta_\mathrm{i}=',stric,'$'], ['',tfsps,', $\theta_\mathrm{i}=',stric,'$'], ['',rfspp,', $\theta_\mathrm{i}=',stric,'$']};
  mark = {'o', '^', 'v', 'o', '^', 'v'; 'o', 's', '^', 'o', 's', '^'};
  f = {1,1,1,0,0,0;1,1,1,0,0,0};

  i=1;
  for j=[1:2]; cplot(x, sf_ortho_err(j, :), mark{i,0+j}, col(5,:), dnam{i,0+j}, f{i,0+j}); hold on; end
  for j=[1:3]; cplot(x, sf_slant_err(j, :), mark{i,3+j}, col(5,:), dnam{i,3+j}, f{i,3+j}); hold on; end
  i=2;
  for j=[1,3]; cplot(x, fs_ortho_err(j, :), mark{i,0+j}, col(4,:), dnam{i,0+j}, f{i,0+j}); hold on; end
  for j=[1:3]; cplot(x, fs_slant_err(j, :), mark{i,3+j}, col(4,:), dnam{i,3+j}, f{i,3+j}); hold on; end

  ll = get(gca,'children');for i=1:numel(ll); set(ll(i), 'linestyle', ls, 'linewidth', 2); end
  plot([min(x),max(x)], [min(x),max(x)].^(-5)*0.03, '-', 'displayname', ['$\propto N^{-5}$'], 'color', col(1,:),'linewidth',4);
  plot([min(x),max(x)], [min(x),max(x)].^(-3)*30, '--', 'displayname', ['$\propto N^{-3}$'], 'color', col(1,:),'linewidth',4);
  set(gca,'xscale','log','yscale','log');
  xlabel('elements per smallest P wavelength ($\lambda_0/\Delta x$)');
  ylabel('relative error [\%]');
  xlim([.9,10]); xticks([.9,1:10]);
  ylim([1e-7,1e2]); yticks(logspace(-7,2,10));
  legend('NumColumns',1,'location','eastoutside');
end

function cplot(x,y,m,c,n,f)
  h = plot(x, y, 'marker', m, 'color', c, 'displayname', n);
  if(f)
    set(h, 'markerfacecolor', 'none', 'markeredgecolor', h.Color);
  else
    set(h, 'markerfacecolor', h.Color, 'markeredgecolor', 'none');
  end
end