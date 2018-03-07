function [] = rewrite_model(NEWDATAFILE, OLDDATAFILE, Z, RHO, TEMP, SOUNDSPEED, P, LOCALPRESSURESCALE, G, NBVSQ, KAPPA, VISCMU, MUVOL, WNORTH, WEAST, W, CP, CV, GAMMA)
  dimensions=[size(Z);size(RHO);size(TEMP);size(SOUNDSPEED);size(P);size(LOCALPRESSURESCALE);size(G);size(NBVSQ);size(KAPPA);size(VISCMU);size(MUVOL);size(WNORTH);size(WEAST);size(W);size(CP);size(CV);size(GAMMA)];
  if(any(dimensions(:,1)/dimensions(1,1)~=1) || any(dimensions(:,2)/dimensions(1,2)~=1))
    error("Dimensions of all datasets must agree.");
  end
  
  decision=-1;
  if(exist(NEWDATAFILE)~=0)
    while(decision~=0 && decision~=1)
      decision=input([char(strcat("File (", NEWDATAFILE,') exists.\nOverwrite? (0 for no, 1 for yes) >')),' ']);
    end
    if(decision==0)
      disp("Overwrite cancelled, stopping script."); return;
    end
  end
  
  f_old = fopen(OLDDATAFILE, 'r');
  if(f_old==-1)
    error(strcat("Cannot open file ", OLDDATAFILE,').'))
  end
  f_new = fopen(NEWDATAFILE, 'w');
  if(f_new==-1)
    error(strcat("Cannot open file ", NEWDATAFILE,').'))
  end
  
  count=1;
  while(count<=3) % Rewrite header.
    count=count+1;
    fprintf(f_new, fgetl(f_old));
    fprintf(f_new, "\n");
  end
  
%z[m] rho[kg/(m^3)] T[K] c[m/s] p[Pa] H[m] g[m/(s^2)] N^2[1/s] kappa[J/(s.m.K)] mu[kg(s.m)] mu_vol[kg/(s.m)] w_M[m/s] w_Z[m/s] w_P[m/s] c_p[J/(mol.K)] c_v[J/(mol.K)] gamma
  for(i=1:numel(Z))
    line_numbers=[Z(i), RHO(i), TEMP(i), SOUNDSPEED(i), P(i), LOCALPRESSURESCALE(i), G(i), NBVSQ(i), KAPPA(i), VISCMU(i), MUVOL(i), WNORTH(i), WEAST(i), W(i), CP(i), CV(i), GAMMA(i)];
    fprintf(f_new,'%15.8e ',line_numbers);
    fprintf(f_new, "\n");
  end
  
  fclose('all');
  disp(strcat("Model output to file: '", NEWDATAFILE, "'."));
end