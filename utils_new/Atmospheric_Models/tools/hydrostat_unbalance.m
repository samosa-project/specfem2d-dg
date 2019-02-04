function [term1, term2, term3]=hydrostat_unbalance(Z, RHO, TEMP, P, G, KAPPA, VISCMU, W)
  D=differentiation_matrix(Z, 0);
  term1 = abs(D * (VISCMU .* (D * W)));
  term2 = abs(D * (KAPPA .* (D * TEMP) + W .* VISCMU .* (D * W)));
  term3 = abs(D * P + RHO .* G);
end