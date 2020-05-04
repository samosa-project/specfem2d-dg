function [bgm, isMatrix] = try_make_bg_model_matrix(bgm)
  [order, ~, ~, ~] = order_bg_model();
  nb_qty = size(order, 1);
  dbgmx = diff(bgm.xx);

  % Other method, but too lax.
%   NZ = numel(unique(round(bgm.zz,12)));
%   NX = numel(unique(round(bgm.xx,12)));
%   if(numel(bgm.xx)==NX*NZ)
%     for iqty = 1:nb_qty
%       bgm.(order(iqty,:)) = reshape(bgm.(order(iqty,:)), NZ, NX);
%     end
%     isMatrix = 1;
%   else
%     isMatrix = 0;
%   end
  
  if(not(any(diff(diff(find(dbgmx>0))))))
    % bgm.xx is a series of repeating xx values
    NZ = find(dbgmx>0, 1, 'first');
    NX = numel(bgm.xx)/NZ;
    for iqty = 1:nb_qty
      bgm.(order(iqty,:)) = reshape(bgm.(order(iqty,:)), NZ, NX);
    end
    isMatrix = 1;
  else
    isMatrix = 0;
  end
end

