function plot_hydrostat_unbalance(Z, RHO, TEMP, P, G, KAPPA, VISCMU, W, tit_plus, save_plots)
  figure();
  [term1, term2, term3]=hydrostat_unbalance(Z, RHO, TEMP, P, G, KAPPA, VISCMU, W);
  semilogx(term1, Z, 'r'); hold on;
  semilogx(term2, Z, 'g');
  semilogx(term3, Z, 'b');
  tmp_valmat=cell2mat(get(get(gca, 'children'), 'XData')); tmp_plot_maxval=max(max(tmp_valmat)); tmp_valmat(tmp_valmat==0)=Inf; tmp_plot_minval=min(min(tmp_valmat));
  xlim([0.5*tmp_plot_minval, 2*tmp_plot_maxval]);
  xlabel({'amplitude of hydrostatic unbalance terms', '(projected wind)'}); ylabel('altitude (m)');
  legend('$\left|\partial_z\left(\mu\partial_zw\right)\right|$', '$\left|\partial_z\left(\kappa\partial_zT+w\mu\partial_zw\right)\right|$', '$\left|\partial_zp + \rho g_z\right|$', 'Location', 'best');
  title(tit_plus);
  if save_plots == 1
    saveas(gcf, strcat(DATAFILE,'__unbalance_terms.png'));
  end
end