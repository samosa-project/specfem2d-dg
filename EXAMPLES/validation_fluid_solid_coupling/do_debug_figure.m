figure('units','normalized','outerposition',[0,0,1,1]);
subplot(131);
plot(incoming_reflected_v*1e-6,'displayname', '$p''$ downscaled *1e-6'); hold on;
plot(transmitted_vx, 'displayname', '$v_x$'); plot(transmitted_vz, 'displayname', '$v_z$');
plot(transmitted_vp, ':', 'displayname', 'P-wave'); plot(transmitted_vs,':', 'displayname', 'S-wave');
legend('location', 'northwest');
title(['i=',num2str(incident_angle)]);
subplot(132);
scatter(transmitted_vx, transmitted_vz, 50, 1:numel(transmitted_vx), 'filled'); hold on;
colormaps_fromPython('hsv', 1); colorbar;
plot([-1,1]*1e-6*cos(incident_angle), [-1,1]*1e-6*sin(incident_angle), 'displayname', 'interface angle');
plot([-1,1]*1e-6*cos(incident_angle-90*pi/180-i2), [-1,1]*1e-6*sin(incident_angle-90*pi/180-i2), ':', 'displayname', 'P angle');
plot([-1,1]*1e-6*cos(incident_angle-90*pi/180-j2), [-1,1]*1e-6*sin(incident_angle-90*pi/180-j2), ':', 'displayname', 'S angle');
daspect([1,1,1]);
legend('location', 'northoutside')
xlabel('vx'); ylabel('vz');
subplot(133);
[th, r] = cart2pol(transmitted_vx, transmitted_vz);
scatter(th*180/pi, r, 50, 1:numel(transmitted_vx), 'filled'); hold on;
colormaps_fromPython('hsv', 1);
plot((incident_angle*180/pi-90)*[1,1],ylim, 'displayname', 'solid normal');
plot((incident_angle*180/pi-90-i2*180/pi)*[1,1],ylim, ':', 'displayname', 'P angle');
plot((incident_angle*180/pi-90-j2*180/pi+90)*[1,1],ylim, ':', 'displayname', '$90^\circ$ anticlockwise from S angle');
legend('location', 'northeast')