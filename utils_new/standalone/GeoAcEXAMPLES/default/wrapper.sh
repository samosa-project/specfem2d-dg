#!/bin/bash
# Author:        Léo Martire.
# Description:   TODO.
# Last modified: TODO.
# Usage:         1) Download GeoAc from https://github.com/LANL-Seismoacoustics/GeoAc.
#                2) 'make all' GeoAc.
#                3) Set path2GeoAc below to your installation path.
#                4) Set version as needed.
#                5) Call using parameters described below.
# Notes:         N/A.

path2GeoAc='/home/l.martire/TMP/GeoAc'
version='3D'
#echo $version

atmFileName='AMGeoAc'

neededArgs=6
if [[ $# -lt $neededArgs ]]; then
  echo
  echo "> $neededArgs arguments needed at least:"
  echo "  - wind projection angle (clockwise from North) [°],"
  echo "  - altitude of source [km],"
  echo "  - minimum launch angle (from ground) [°],"
  echo "  - maximum launch angle (from ground) [°],"
  echo "  - launch angle step [°],"
  echo "  - frequency [Hz]."
  echo "> Optionnal arguments:"
  echo "  - minimum number of bounces."
  echo "> Script will now exit."
  echo
  exit -1
fi

# Bring plotting tools.
cp ../plot_utils/* .

# Parse inputs.
proj=$(($1 % 360))
zs=$2
tmi=$3
tma=$4
ts=$5
b=$6

# Echo out parsed input for user checking.
echo " "
echo "  wind projection angle (clockwise from North) $proj [°],"
echo "  altitude of source                 $zs [km],"
echo "  minimum launch angle (from ground) $tmi [°],"
echo "  maximum launch angle (from ground) $tma [°],"
echo "  launch angle step                  $ts [°]."
echo "  minimum number of bounces          $b."
echo " "

# Run in the direction asked.
#$path2GeoAc/GeoAc$version -prop $atmFileName.met z_src=$zs theta_step=$ts theta_min=$tmi theta_max=$tma bounces=$b phi_min=$proj phi_max=$proj phi_step=180
$path2GeoAc/GeoAc$version -prop $atmFileName.met z_src=$zs theta_step=$ts theta_min=$tmi theta_max=$tma bounces=$b azimuth=$proj
# Copy out output.
cp ${atmFileName}_raypaths.dat ${atmFileName}_raypaths_$proj.dat
# Plot.
gnuplot -p "plotRays$version"

# Run in the opposite direction of the direction asked.
#echo $((($proj+180) % 360))
newproj=$((($proj+180) % 360))
#$path2GeoAc/GeoAc$version -prop $atmFileName.met z_src=$zs theta_step=$ts theta_min=$tmi theta_max=$tma bounces=$b phi_min=$newproj phi_max=$newproj phi_step=180
echo "$path2GeoAc/GeoAc$version -prop $atmFileName.met z_src=$zs theta_step=$ts theta_min=$tmi theta_max=$tma bounces=$b azimuth=$newproj"
$path2GeoAc/GeoAc$version -prop $atmFileName.met z_src=$zs theta_step=$ts theta_min=$tmi theta_max=$tma bounces=$b azimuth=$newproj
# Copy out output.
cp ${atmFileName}_raypaths.dat ${atmFileName}_raypaths_$newproj.dat
# Plot.
gnuplot -p "plotRays$version"
