function feg = MMS_plotFields(err_l2, L2NormOfError, Xe, Ye, zth, zexp, zerr, IT, DT, NX, synthname, YLAB_CB)
  % Parameters
%   TIT_SYNTH = [synthname,' ($t=',scientific_latex_notation(DT*IT,2),'$~s)'];
  TIT_SYNTH = [synthname,' ($t=',sprintf('%.2f',DT*IT),'$~s)'];
  XLAB = '$x$ [m]';
  YLAB = '$z$ [m]';
  NColorbarTicks = 9;
  addErrorPanel = 1; % Set to 1 for debug, set to 0 for printing nice plots
  CMAP_thresh = 0.01; % relative to max
  
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
  xcm=[-1, -0.85, -0.15, CMAP_thresh]; xcm=[xcm, 0, -fliplr(xcm)];
  N=1000; % number of points for color values
  CMAP = interp1(xcm, CMAP, linspace(min(xcm),max(xcm),N)');
%   colormap(CMAP);

  % Colorbar scale.
  caxis_exp = [min(min(zexp)), max(max(zexp))];
  caxis_the = max(max(abs(zth)))*[-1,1];
  caxis_err = [min(min(zerr)), max(max(zerr))];
  all_caxx = [caxis_the; caxis_exp; caxis_err];
  glob_caxx = caxis_the;
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
  if(addErrorPanel)
    nbSubplots = 3;
  else
    nbSubplots = 2;
  end
  
  feg = figure('units','normalized','outerposition',[0 0 1 0.72]);
  axxx = tight_subplot(1, nbSubplots, 0.012, [0.12,0.06], [0.06, 0.12]); % not mandatory, but prettier

  % Analytic. %%%%
  axes(axxx(1));
%   surf(Xe,Ye,zth); view([0,0,1]);
  pcolor(Xe,Ye,zth);
  shading interp; colormap(CMAP); caxis(glob_caxx); %colorbar;
  xlabel(XLAB); ylabel(YLAB);
  title('Analytic');
  daspect([1,1,1]);

  % Synthetic. %%%
  axes(axxx(2));
%   surf(Xe,Ye,zexp); view([0,0,1]);
  pcolor(Xe,Ye,zexp);
  shading interp; colormap(CMAP); caxis(glob_caxx); %colorbar;
  yticklabels([]);
  xlabel(XLAB);
  %       title([synthname,' ($n=',num2str(IT),'$)']);
  title(TIT_SYNTH);
  daspect([1,1,1]);
  
  if(addErrorPanel)
  % Error. %%%%%%%
    axes(axxx(3));
%     surf(Xe,Ye,zerr); view([0,0,1]);
    pcolor(Xe,Ye,zerr);
    shading interp; colormap(CMAP); caxis(glob_caxx);
    xlabel(XLAB);
    yticklabels([]);
    title({[err_l2.prefix,'(',num2str(NX),')=',scientific_latex_notation(L2NormOfError,2),'$~',err_l2.unit]});
    daspect([1,1,1]);
  end
  
  pause(0.1);
  
  % Color bar. %%%
  posAxxx3 = get(axxx(nbSubplots),'Position');
  h_cb = colorbar('Position', [posAxxx3(1)+posAxxx3(3)+0.01  posAxxx3(2)  0.03  posAxxx3(4)]);
  colorbarticks=linspace(-1,1,NColorbarTicks); %colorbarticks=sort([colorbarticks, [-1,1]*max(thresh)]);
  if(max(max(abs(zth)))>0)
    colorbarticks = colorbarticks*max(max(abs(zth)));
  else
    colorbarticks = colorbarticks*max(max(abs(zexp)));
  end
%   colorbarticks = colorbarticks;
  set(h_cb, 'ticks', colorbarticks);
  
  YTL = split(sprintf('%.3f|',colorbarticks),'|'); YTL(end)=[];
  set(h_cb, 'ticklabels',YTL);
  
  ylabel(h_cb,YLAB_CB,'interpreter','latex', 'fontsize',26);
  
  set(h_cb,'tickdir','both');
  
  ll = add_labels_subplots(feg, 0.9, 0, [0, -0.0375]);

  % Link = linkprop(axxx,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
%   Link = linkprop(axxx,{'CameraPosition', 'XLim', 'YLim'});
%   setappdata(gcf, 'StoreTheLink', Link);
%   prettyAxes(feg);
end