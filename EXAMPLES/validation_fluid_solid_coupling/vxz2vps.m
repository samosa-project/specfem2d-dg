function [vp, vs] = vxz2vps(vx, vz, ig, i2, j2)
  tn2xz = [cos(ig), -sin(ig);
           sin(ig), cos(ig)];
  xz2tn = inv(tn2xz);
  
  vt = xz2tn(1,1)*vx + xz2tn(1,2)*vz; % motion tangent to interface (towards positive x)
  vn = xz2tn(2,1)*vx + xz2tn(2,2)*vz; % motion normal to interface (towards fluid)
  
  [argv, modv] = cart2pol(vt, vn);
  
  vp = modv .* cos(argv-0.5*pi+i2); % motion along i2 positive towards fluid
  vp = -vp; % motion along i2 positive towards solid
  
%   vs = modv .* cos(argv-0.5*pi+j2); % motion along j2 positive towards fluid
%   vs = -vs; % motion along j2 positive towards solid
  vs = modv .* cos(argv-0.5*pi+j2+0.5*pi); % motion transverse to j2 (j2+90) positive towards fluid
  vs = -vs; % motion along j2 positive towards solid
  
%   vp = sin(i2)*vx + cos(i2)*vz; % ground velocity along P-wave direction
%   vs = cos(j2)*vx + sin(j2)*vz; % ground velocity 90° from S-wave direction % 90°
%   vs = sin(j2)*vx + cos(j2)*vz; % ground velocity 90° from S-wave direction % 0°
end

