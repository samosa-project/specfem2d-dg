function [hh] = draw_ampl_circle(ampl, ang_centre, factor)
  ang_width = factor * 2/ampl;
  
  ang = linspace(ang_centre-0.5*ang_width, ang_centre+0.5*ang_width, 10);
%   ang = [ang, ang+pi]; % both sides
  
  [x, z] = pol2cart(ang, ang*0 + ampl);
  hold on;
  hh = plot(x, z);
end

