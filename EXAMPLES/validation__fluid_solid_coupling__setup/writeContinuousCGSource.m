
function [N] = writeContinuousCGSource(source_out, xmin, xmax, zmin, zmax, f0, endpoints_x, endpoints_z)
  N = ceil(max(xmax-xmin, zmax-zmin)/0.5);
%   disp(['[',mfilename,'] ',num2str(N),' sources generated. Set parfile accordingly.']);
  
  xsarr = linspace(endpoints_x(1), endpoints_x(end), N);
  zsarr = linspace(endpoints_z(1), endpoints_z(end), N);
  
%   amp = 1;
%   amp = exp(-(xsarr/((xmax-xmin)/6)).^2);
  l=50; amp = 0.5*(cos(2*pi*(xsarr-xmin+l)/(2*l))+1).*(xsarr<=xmin+l) + (xsarr>=xmin+l).*(xsarr<=xmax-l) + 0.5*(cos(2*pi*(xsarr-xmax+l)/(2*l))+1).*(xsarr>=xmax-l);
  
  amp = 1e3*amp;
  
%   source_out = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_transmission/source_input_solid';
  
  fid = fopen(source_out, 'w');
  
  for i=1:N
    fprintf(fid,'source_surf                     = .false.    # Source at the surface? If set to false, source is inside the medium.\n');
    fprintf(fid,'xs                              = %7.2fd0  # Source x location (m).\n', xsarr(i));
    fprintf(fid,'zs                              = %7.2fd0  # Source y location (m).\n', zsarr(i));
    fprintf(fid,'source_type                     = 1          # Elastic force or acoustic pressure = 1, moment tensor = 2. For the moment, only 1 is supported by DG.\n');
    fprintf(fid,'time_function_type              = 3          # Source time function type: Ricker = 1, Gaussian first derivative = 2, Gaussian = 3, Dirac = 4, Heaviside = 5, 6 = ??, 7 = ??, 8 = read from file, 9 = burst, 10 = Gaussian primitive.\n');
    fprintf(fid,'name_of_source_file             = "/path/to/file"\n');
    fprintf(fid,'burst_band_width                = 200.d0     # Only for time_function_type 9. Band width of the burst.\n');
    fprintf(fid,'f0                              = %7.2fd0\n', f0);
    fprintf(fid,'tshift                          = 0.d0\n');
    fprintf(fid,'anglesource                     = 0.d0       # Only for a force. Angle of the source (degrees). 0 is upwards, 180 is downwards.\n');
    fprintf(fid,'Mxx                             = 1.d0       # Only for moment tensor source. Mxx component.\n');
    fprintf(fid,'Mzz                             = 1.d0       # Only for moment tensor source. Mzz component.\n');
    fprintf(fid,'Mxz                             = 0.d0       # Only for moment tensor source. Mxz component.\n');
    fprintf(fid,'factor                          = %7.2fd0     # Amplification factor.\n', amp(i));
    fprintf(fid,'###############################\n');
  end
  
  fclose(fid);
end