figcheckhandle=figure();
% colours=jet(nstat);
% colours=prism(nstat);
colours=winter(nstat);
ratio_anal = peak2peak(synf(end,:))/peak2peak(synf(1,:));
ratio_exp = peak2peak(Zamp(end,:))/peak2peak(Zamp(1,:));
if(externalDGAtmosModel)
  error('no th value for external atmos model');
else
  if(USE_ISOTHERMAL_MODEL)
    th_vz_v0 = exp( range(zstattab)/(2*H) );
    
    % Look at agreement with exponential decay.
    factor = 1; facttxt = '1';
%     factor = exp(2*zstattab/(SOUNDSPEED^2)); facttxt = 'exp(2*z/c$^2$)';
%     factor = exp(1.446e-5 * zstattab); facttxt = 'exp(1.446e-5 * z)';
    isothDecay = [zstattab, factor.* peak2peak(synf(:,:),2)./peak2peak(synf(1,:)), exp( zstattab/(2*H) ), peak2peak(Zamp(:,:),2)/peak2peak(Zamp(1,:))];
    x=zstattab; y=isothDecay(:,4)./isothDecay(:,3); lel = fit(x,y,'exp1');
%     if(log10(abs(lel.b))>=-6)
    if(abs(ratio_anal-ratio_exp)/ratio_anal>0.01)
      % If a non-negligible multiplying factor is found, plot and display.
%       lel
%       sprintf('%.16f',lel.b)
      figure();
      plot(zstattab/1e3, isothDecay(:,3), '-o', 'markersize', 10, 'displayname', ['exponential decay $\exp(z/(2H))$ $H=',num2str(H/1e3),'$~km']); hold on;
      kek = fit(x,isothDecay(:,4),'exp1');
      plot(zstattab/1e3, isothDecay(:,4), '-o', 'markersize', 10,  'displayname', ['synthetics, fit $\exp(z/(2H))$ yield $H=',num2str((0.5/kek.b)/1e3),'$~km']); hold on;
      kuk = fit(x,isothDecay(:,2),'exp1');
      plot(zstattab/1e3, isothDecay(:,2), '-o', 'markersize', 10,  'displayname', ['analytic, fit $\exp(z/(2H))$ yield $H=',num2str((0.5/kuk.b)/1e3),'$~km']); hold on;
      legend();
      ylabel('$v_z(z)/v_z(z=0)$ ratio [1]');
      xlabel('$z$ [km]');
      legend('location', 'northoutside');
    end
  else
    th_vz_v0 = 1;
  end
end
figure(figcheckhandle);
for i = 1:nstat;
  plot(t,synf(i,:),'displayname',num2str(zstattab(istattab(i))),'color',colours(i,:)); hold on;
end
title({['This analytical solution: $v(z=z_\mathrm{max})/v(z=0)$=',num2str(ratio_anal),'.'],...
       ['Theoretical value: =$e^{z/(2H)}$=',num2str(th_vz_v0),'.'],...
       ['Current simulation: =',num2str(ratio_exp),'.']},'fontsize',14);
legend();
prettyAxes(figcheckhandle);
isThisOk=-1;
while(not(ismember(isThisOk,[0,1])))
  isThisOk=input(['[',mfilename,'] Does this analytical solution look ok (0 for no, 1 for yes)? > ']);
end
if(not(isThisOk))
  error('stop kek');
end