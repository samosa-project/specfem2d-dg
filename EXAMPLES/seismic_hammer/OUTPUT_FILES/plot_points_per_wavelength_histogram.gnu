 set term wxt
 #set term gif
 #set output "points_per_wavelength_histogram_S_in_solid.gif"

 set boxwidth    1.55795789    
 set xlabel "Range of min number of points per S wavelength in solid"
 set ylabel "Percentage of elements (%)"
 set loadpath "./OUTPUT_FILES"
 plot "points_per_wavelength_histogram_S_in_solid.txt" with boxes
 pause -1 "hit any key..."
 #set term gif
 #set output "points_per_wavelength_histogram_P_in_fluid.gif"

 set boxwidth   0.173799202    
 set xlabel "Range of min number of points per P wavelength in fluid"
 set ylabel "Percentage of elements (%)"
 set loadpath "./OUTPUT_FILES"
 plot "points_per_wavelength_histogram_P_in_fluid.txt" with boxes
 pause -1 "hit any key..."
