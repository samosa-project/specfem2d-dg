figcheckhandle=figure();
% colours=jet(nstat);
% colours=prism(nstat);
colours=winter(nstat);
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
    x=zstattab; y=isothDecay(:,3)./isothDecay(:,2); lel = fit(x,y,'exp1');
    if(log10(abs(lel.b))>=-6)
      % If a non-negligible multiplying factor fit is found, plot and display.
      lel
      sprintf('%.16f',lel.b)
      figure();
      plot(zstattab, isothDecay(:,3), 'displayname', ['exponential decay']); hold on;
      plot(zstattab, isothDecay(:,4), 'displayname', ['synthetics']); hold on;
      plot(zstattab, isothDecay(:,2), 'displayname', ['analytic $\times$',facttxt]); hold on;
      legend();
    end
  else
    th_vz_v0 = 1;
  end
end
figure(figcheckhandle);
for i = 1:nstat;
  plot(t,synf(i,:),'displayname',num2str(zstattab(istattab(i))),'color',colours(i,:)); hold on;
end
title({['This analytical solution: $v(z=z_{max})/v(z=0)$=',num2str(peak2peak(synf(end,:))/peak2peak(synf(1,:))),'.'],...
       ['Theoretical value: =$e^{z/(2H)}$=',num2str(th_vz_v0),'.'],...
       ['Current simulation: =',num2str(peak2peak(Zamp(end,:))/peak2peak(Zamp(1,:))),'.']},'fontsize',14);
legend();
prettyAxes(figcheckhandle);
isThisOk=-1;
while(not(ismember(isThisOk,[0,1])))
  isThisOk=input(['[',mfilename,'] Does this analytical solution look ok (0 for no, 1 for yes)? > ']);
end
if(not(isThisOk))
  error('stop kek');
end