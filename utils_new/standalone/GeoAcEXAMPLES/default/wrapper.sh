#!/bin/bash
# Author:        Léo Martire.
# Description:   TODO.
# Last modified: TODO.
# Usage:         N/A.
# Notes:         N/A.

path2GeoAc='/home/l.martire/TMP/GeoAc'
version='3D'
#echo $version

if [[ $# -lt 5 ]]; then
  echo
  echo "> 5 arguments needed at least:"
  echo "  - wind projection angle (clockwise from North) [°],"
  echo "  - altitude of source [km],"
  echo "  - minimum launch angle (from ground) [°],"
  echo "  - maximum launch angle (from ground) [°],"
  echo "  - launch angle step [°]."
  echo "> Optionnal arguments:"
  echo "  - minimum number of bounces."
  echo "> Script will now exit."
  echo
  exit -1
fi

# bring plotting tools
cp ../plot_utils/* .

# parse
proj=$(($1 % 360))
zs=$2
tmi=$3
tma=$4
ts=$5

b=$6

echo " "
echo "  wind projection angle (clockwise from North) $proj [°],"
echo "  altitude of source                 $zs [km],"
echo "  minimum launch angle (from ground) $tmi [°],"
echo "  maximum launch angle (from ground) $tma [°],"
echo "  launch angle step                  $ts [°]."
echo "  minimum number of bounces          $b."
echo " "

$path2GeoAc/GeoAc$version -prop AMGeoAc.met z_src=$zs theta_step=$ts theta_min=$tmi theta_max=$tma bounces=$b phi_min=$proj phi_max=$proj phi_step=180
gnuplot -p "plotRays$version"

#echo $((($proj+180) % 360))

newproj=$((($proj+180) % 360))
$path2GeoAc/GeoAc$version -prop AMGeoAc.met z_src=$zs theta_step=$ts theta_min=$tmi theta_max=$tma bounces=$b phi_min=$newproj phi_max=$newproj phi_step=180
gnuplot -p "plotRays$version"
