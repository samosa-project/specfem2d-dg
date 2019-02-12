% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         It might not be a good idea to run this script alone. At the very least, anareg_atmospheric_model.m should have been ran before.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bruteforce Density.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: bruteforce density such that hydrostatic gap is zero.
if(strcmp(method, 'bruteforce_rho'))
  D=differentiation_matrix(Z, 0);
  bruteforced_RHO = - (D * P) ./ G;
end
if(strcmp(method, 'bruteforce_rho_log'))
  D=differentiation_matrix(Z, 0);
  bruteforced_RHO_log = - ((D * (log(P))) .* P ) ./ G;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate Pressure.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Integrate the model's RHS ($-\rho g$) along z using an
% interative scheme in order to retrieve p. Do this smartly, see 
% reference below.
% Ref.: http://acoustique.ec-lyon.fr/publi/haniquecockenpot_thesis.pdf.
if(strcmp(method, 'integrate'))
  D=differentiation_matrix(Z, 1);
  nz=numel(Z);
  itermax = 1e5;
  err = ones(nz, 1);
  residual = zeros(nz, 1);
  model = log(P);
  iter = 0;
  % res_evol = [];
  err_evol = [];
  hydrostatic_RHS_log = - RHO .* G ./ P;
  hydrostatic_RHS_log(1) = log(P(1)); % Has to be in accordance with conditionning in the differentiation matrix.
  % hydrostatic_RHS_log(1) = log(P(end)); % Has to be in accordance with conditionning in the differentiation matrix.
  tic
  Dml = inv(D);
  while norm(err) / norm(model) > 10 ^ (- 15) && iter <= itermax
    residual = hydrostatic_RHS_log - (D * model);
    err = Dml * residual;
    model = model + err;
    iter = iter + 1;
  %   res_evol = [res_evol; norm(residual)];
    err_evol = [err_evol; norm(err)];
  end
  toc

  figure();
  loglog(1:iter, err_evol);
  limsy = get(gca, 'YLim'); ylim([limsy(1), 2 * max(err_evol)]); xlim([1, iter]);
  xlabel('iterations'); ylabel('optimisation process error');
  title(tit_plus);

  optimised_P = exp(model);

  % figure();
  % ylabel('altitude (m)');
  % semilogx(hydrostatic_ratio(D, P, RHO, G), Z, hydrostatic_ratio(D, optimised_P, RHO, G), Z, ones(size(Z)), Z, 'k:');
  % xlabel('ratio'); ylabel('altitude (m)');
  % legend('before optimisationon $p$', 'after optimisation on $p$', 'Location', 'best');
  % title('Hydrostatic Ratio $-\partial_zp/\rho g_z$');

  figure();
  semilogx(abs(D * P + RHO .* G) ./ (RHO .* G), Z, abs(D * optimised_P + RHO .* G) ./ (RHO .* G), Z);
  xlabel('relative gap $|\partial_z p + \rho g|$/$\rho g$'); ylabel('altitude (m)');
  legend('model', 'after integration', 'Location', 'best');
  title(tit_plus);
  if save_plots == 1
    saveas(gcf, strcat(DATAFILE,'__regularisation_by_integration.png'));
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use a metaheuristic.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: bruteforce density such that hydrostatic gap is zero.
if(strcmp(method, 'metaheuristic'))
  disp(['[',mfilename,', WARNING] Metaheuristic methods can be very long.']);
  callFromOutside=1;
  algo = 'ps';
  ftol = 1e-16;
  maxit = 10000;
  modify_atmos_model_metaheuristics;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%