% Prepare colorbar.
colorOK = [1,1,1]; % colour for |value|<thresh
colorZero = [255,228,181]/255; % moccasin
% CMAP = 'jet'; % default colorbar
CMAP = [0,0,0.25;
        0,0,0.5;
        0,0,1;
        colorOK; colorZero; colorOK;
        1,0,0;
        0.5,0,0;
        0.25,0,0];
xcm=[-1, -0.85, -0.15, plotFields_CMAP_thresh]; xcm=[xcm, 0, -fliplr(xcm)];
N=1000; % number of points for color values
CMAP = interp1(xcm, CMAP, linspace(min(xcm),max(xcm),N)');
colormap(CMAP);

% Colorbar scale.
all_caxx = [caxxxxiiss_th;caxxxxiiss_exp;caxxxxiiss_err];
glob_caxx = caxxxxiiss_th;
if(numel(unique(glob_caxx))==1)
  % If zth = cste, choose glob_caxx based on zexp
  glob_caxx = max(max(abs(zexp)))*[-1,1];
  % If again zexp = cste, choose glob_caxx arbitrarily
  if(numel(unique(glob_caxx))==1)
    if(glob_caxx(1)==0)
      % If cste=0, choose arbitrarily a maximum and minimum
      glob_caxx = [-1,1]*0.1;
    else
      glob_caxx = glob_caxx(1)*[0.9, 1.1];
    end
  end
end

% Figure.
if(plotFields_addErrorPanel)
  nbSubplots = 3;
else
  nbSubplots = 2;
end
feg = figure('units','normalized','outerposition',[0 0.4 1 0.6]);
axxx = tight_subplot(1, nbSubplots, 0.012, [0.14,0.08], [0.05, 0.12]); % not mandatory, but prettier

% Analytic. %%%%
axes(axxx(1));
surf(Xe,Ye,zth);
shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx); %colorbar;
xlabel(plotFields_xlab); ylabel(plotFields_ylab);
title('Analytic');

% Synthetic. %%%
axes(axxx(2));
surf(Xe,Ye,zexp);
shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx); %colorbar;
yticklabels([]);
xlabel(plotFields_xlab);
%       title([synthname,' ($n=',num2str(IT),'$)']);
title([synthname,' ($t=',scientific_latex_notation(DT*IT,2),'$~s)']);

if(plotFields_addErrorPanel)
% Error. %%%%%%%
  axes(axxx(3));
  surf(Xe,Ye,zerr);
  shading interp; view([0,0,1]); colormap(CMAP); caxis(glob_caxx);
  xlabel(plotFields_xlab);
  yticklabels([]);
  title({[errName,'=',scientific_latex_notation(L2NormOfError,2),'$']});
end

% Color bar. %%%
posAxxx3 = get(axxx(nbSubplots),'Position');
h_cb = colorbar('Position', [posAxxx3(1)+posAxxx3(3)+0.03  posAxxx3(2)  0.03  posAxxx3(4)]);
colorbarticks=linspace(-1,1,plotFields_NColorbarTicks); %colorbarticks=sort([colorbarticks, [-1,1]*max(thresh)]);
if(max(max(abs(zth)))>0)
  colorbarticks = colorbarticks*max(max(abs(zth)));
else
  colorbarticks = colorbarticks*max(max(abs(zexp)));
end
colorbarticks = colorbarticks;
set(h_cb, 'ytick', colorbarticks);
set(h_cb,'tickdir','both');

% Link = linkprop(axxx,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
Link = linkprop(axxx,{'CameraPosition', 'XLim', 'YLim'});
setappdata(gcf, 'StoreTheLink', Link);
prettyAxes(feg);
customSaveFig(feg, savefigfullpath,plotFields_extsToSave);