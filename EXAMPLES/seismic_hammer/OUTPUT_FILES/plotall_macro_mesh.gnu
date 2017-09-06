 set term wxt
 # set term postscript landscape color solid "Helvetica" 22
 # set output "macro_mesh.ps"
 set xlabel "X"
 set ylabel "Y"
 set title "Spectral Element (Macrobloc) Mesh"
 set size ratio -1
 set loadpath "./OUTPUT_FILES"
 plot "macros2.gnu" title '' w l lc 2, "macros1.gnu" title '' w l lc 3
 pause -1 "Hit any key to exit..."
