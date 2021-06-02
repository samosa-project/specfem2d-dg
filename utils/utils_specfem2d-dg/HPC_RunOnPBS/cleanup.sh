#!/bin/bash

JOBTAG="SPECFEM_JOB"
slurmFile="batchJobScript.sh"

curDir=$(pwd)
cd $1

if [[ -d "OUTPUT_FILES" ]]; then
id=$(ls ${JOBTAG}* | grep -oe "[0-9]\+" | head -n 1)
mv $slurmFile OUTPUT_FILES
mv ${JOBTAG}* OUTPUT_FILES
mv OUTPUT_FILES OUTPUT_FILES_$id
fi

cd $curDir
