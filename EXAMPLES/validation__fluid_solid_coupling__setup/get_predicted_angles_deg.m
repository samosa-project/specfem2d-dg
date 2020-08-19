function [i1_i2_j2t_j2r] = get_predicted_angles_deg(ic_rad, alpha__1, alpha__2, beta__2)
%   i1_i2_j2t_j2r = asin([c*sin(ic)/vp, vp*sin(ic)/c, vs*sin(ic)/c, vs*sin(ic)/vp])* 180/pi;
  
  i1_i2_j2t_j2r = [snells(alpha__2, alpha__1, ic_rad), snells(alpha__1, alpha__2, ic_rad), snells(alpha__1, beta__2, ic_rad), snells(alpha__2, beta__2, ic_rad)]* 180/pi;
  
end