#!/bin/bash
# Author:        Léo Martire.
# Description:   Replaces NSTEP in all databases files in a folder.
# Last modified: See file metadata.
# Usage:         N/A.
# Notes:         N/A.

if [[ $# -ne 2 ]]; then
  echo
  echo "> 2 arguments needed:"
  echo "  - path to an OUTPUT_FILES folder already containing 'Database*' files."
  echo "  - the new NSTEP value."
  echo "> Script will now exit."
  echo
  exit -1
fi

OUTPUT_FILES=$1
newNSTEP=$2

curentfolder=$(pwd)

cd $OUTPUT_FILES


for file in Database*
do
  perl -0777 -i.ORIGINAL -pe 's/NSTEP DT\n *[0-9]+ +/NSTEP DT\n '$newNSTEP' /igs' $file
done

cd $curentfolder
