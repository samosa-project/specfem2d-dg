#!/bin/bash
# Author:        Léo Martire.
# Description:   Replaces a value in all databases files in a folder.
# Last modified: See file metadata.
# Usage:         N/A.
# Notes:         N/A.

if [[ $# -ne 2 ]]; then
  echo
  echo "> 2 arguments needed:"
  echo "  - path to an OUTPUT_FILES folder already containing 'Database*' files."
  echo "  - the new value."
  echo "> Script will now exit."
  echo
  exit -1
fi

header="NSTEP DT"
format=" *[0-9]+ +"

OUTPUT_FILES=$1
new_value=$2

curentfolder=$(pwd)

cd $OUTPUT_FILES


for file in Database*
do
  perl -0777 -i.ORIGINAL -pe "s/$header\n$format/$header\n $new_value /igs" $file
done

cd $curentfolder
