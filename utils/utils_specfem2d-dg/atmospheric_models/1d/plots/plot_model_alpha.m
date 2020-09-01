% Author:        Léo Martire.
% Description:   Plots an atmospheric model absorption coefficients.
% Notes:         Needs:
%                a) .m scripts and functions (if not alongside this
%                   script, recover via Léo):
%                  1) frsvib2tausigtaueps.m
%                  2) relax_alphavib.m
%                  3) relax_alphavol.m
%
% Usage:
%   TODO.
% with:
%   TODO.
% yields:
%   TODO.

function [fh] = plot_model_alpha(freq, Z, RHO, C, MUvol, FR, SVIB, zminmax, Nsel)
  [SPCFMEXloc] = setup_overall();
  
  zminmax = sort(zminmax);
  [TAU_SIG, TAU_EPS] = frsvib2tausigtaueps(FR, SVIB);
  ALPHAVIB = relax_alphavib(freq, C, TAU_EPS, TAU_SIG);
  ALPHAVOL = relax_alphavol(freq, RHO, C, MUvol);
  selalt = ceil(linspace(find(Z >= zminmax(1), 1, 'first'),find(Z <= zminmax(end), 1, 'last'), Nsel));
  colours = winter(numel(selalt));
  fh=figure('units','normalized','outerposition',[0 0 0.5 1]);
  for i=1:Nsel
    ialt = selalt(i);
    loglog(freq,ALPHAVIB(ialt,:),'color',colours(i,:),'displayname',['$\alpha_{vib}\left(z=',num2str(Z(ialt)),'\right)$']); hold on;
  end
  for i=1:Nsel
    ialt = selalt(i);
    loglog(freq,ALPHAVOL(ialt,:),':','color',colours(i,:),'displayname',['$\alpha_{vol}\left(z=',num2str(Z(ialt)),'\right)$']); hold on;
  end
%   legend('location', 'best');
  legend('location', 'southeast');
  miny=min(min(min(ALPHAVIB)),min(min(ALPHAVOL)));
  miny=10^floor(log10(miny));
  xlim([min(freq),max(freq)]); ylim([miny,1]);
  xlabel('$f$ [Hz]'); ylabel('$\alpha$ [Np/m]');
  set(gca,'ticklabelinterpreter','latex');
  grid on;
  box on;
  title(['Absorption Coefficient']);
end

