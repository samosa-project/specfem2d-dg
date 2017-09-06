 set term wxt
 # set term postscript landscape color solid "Helvetica" 22
 # set output "gll_mesh.ps"
 set xlabel "X"
 set ylabel "Y"
 set title "Gauss-Lobatto-Legendre Mesh"
 set size ratio -1
 set loadpath "./OUTPUT_FILES"
 plot "gllmesh1.gnu" title '' w l lc 2, "gllmesh2.gnu" title '' w l lc 3
 pause -1 "Hit any key to exit..."
