function [i1_i2_j2t_j2r] = get_predicted_angles_deg(ic, c, vp, vs)
%   i1_i2_j2t_j2r = asin([c*sin(ic)/vp, vp*sin(ic)/c, vs*sin(ic)/c, vs*sin(ic)/vp])* 180/pi;
  
  i1_i2_j2t_j2r = [snells(vp, c, ic), snells(c, vp, ic), snells(c, vs, ic), snells(vp, vs, ic)]* 180/pi;
  
end