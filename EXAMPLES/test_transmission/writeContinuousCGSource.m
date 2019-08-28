
function [] = writeContinuousCGSource(xmin, xmax, zmin, zmax, f0)
  N = ceil(max(xmax-xmin,zmax-zmin)/0.75);
  xarr = linspace(xmin, xmax, N);
  zarr = linspace(zmin, zmax, N);
  
  source_out = '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_transmission/source_input_solid_test';
  
  fid = fopen(source_out, 'w');
  
  for i=1:N
    fprintf(fid,'source_surf                     = .false.    # Source at the surface? If set to false, source is inside the medium.\n');
    fprintf(fid,'xs                              = %7.2fd0  # Source x location (m).\n', xarr(i));
    fprintf(fid,'zs                              = %7.2fd0  # Source y location (m).\n', zarr(i));
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
    fprintf(fid,'factor                          = -1.d3      # Amplification factor.\n');
    fprintf(fid,'###############################\n');
  end
  
  fclose(fid);
end