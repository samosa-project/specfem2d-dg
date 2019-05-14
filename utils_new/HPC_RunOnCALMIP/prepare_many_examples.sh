#!/bin/bash
# Author:        Léo Martire.
# Description:   Prepares many EXAMPLES at once.
# Last modified: See file metadata.
# Usage:         N/A.
# Notes:         N/A.


nargneed=1
if [[ $# -lt $nargneed ]]; then
  echo
  echo "> This script needs at least $nargneed EXAMPLE folder, and may take more."
  echo "> Script will now exit."
  echo
  exit -1
fi

echo " "
echo " "
echo " "

echo "This will prepare the following EXAMPLES:"
for var in "$@"
do
  var=$(basename $var)
  echo "  $var"
done

echo " "

read -n 1 -s -r -p "Press any key to continue, or <CTRL+C> to stop. > "
echo ""

TEMPDIR=$(mktemp -d -t) # prepare a temporary directory
#echo $TEMPDIR

for var in "$@"
do
  var=$(basename $var)
  echo " "
  echo "################################"
  echo "# Treating $var"
  echo "################################"
  tmpname=$TEMPDIR"/"$var
  outfile=$tmpname"_out"
  errfile=$tmpname"_err"
  ./prepare_example_new.sh $var > $outfile 2>$errfile
  #echo "Printed output in '"$tmpname"_{out,err}'."
  cat $errfile # show eventual errors
done

echo " "
echo " "
echo " "
