#!/bin/bash

curDir=$(pwd)

cd $1

tar -cvzf OUTPUT_FILES_X.tar.gz input_* image* AA.*$2.* STATIONS PBSJOB*

cd $curDir

echo ""
echo "rsync -avP -e 'ssh' lmartire@gattaca:$(pwd)/$1/OUTPUT_FILES_X.tar.gz $1/"
echo ""
