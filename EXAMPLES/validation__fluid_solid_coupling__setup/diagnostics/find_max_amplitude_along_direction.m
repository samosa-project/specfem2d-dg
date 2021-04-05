function [vals] = find_max_amplitude_along_direction(vx, vz, angs)
  do_debug_fig = 1;
  
  vals = [];
  [th, r] = cart2pol(vx, vz);

  if(do_debug_fig)
    figure();
    subplot(211);
    plot(th); hold on;
    subplot(212);
    plot(r); hold on;
  end
  
  for i=1:numel(angs)
    ang = angs(i);
    cr = r;

    ang1 = wrapToPi(ang);
    ang2 = wrapToPi(ang+pi);
    zero_outside_angle_of_interest = (log(min([abs(th-ang1), abs(th-(ang2))]')) >= -2);
    cr(zero_outside_angle_of_interest) = 0;
    
    if(do_debug_fig)
      subplot(211);
      plot(zero_outside_angle_of_interest);
      subplot(212);
      plot(cr);
    end

    [~, vals(i)] = findFirstPeak(1:numel(cr), cr);
  end
end