#!/bin/bash
# Author:        Léo Martire.
# Description:   Reads an output folder and extracts all available information.
# Last modified: See file metadata.
# Usage:         N/A.
# Notes:         N/A.

if [[ $# -ne 1 ]]; then
  echo
  echo "> 1 argument needed:"
  echo "  - path to an OUTPUT_FILES folder containing a 'slurm' (CALMIP) or '.output' (CEA) file."
  echo "> Script will now exit."
  echo
  exit -1
fi

folder=$1

echo "  Analysing folder '$folder'."

slurm=$(ls $folder/slurm*)
if [ $? -ne 0 ]; then
  slurm=$(ls $folder/*.output)
  if [  $? -ne 0 ]; then
    exit
  fi
fi

parfile=$(ls $folder/* | grep -P --regexp=".*[Pp]ar[_]?file") # Attemps to locate a parfile.

echo "  SLURM file: '$slurm'."
echo "  Parfile:    '$parfile'."

nelemtot=$(grep "Total number of elements:" $slurm | grep -o "[0-9]*")
nelemacous=$(grep "Total number of acoustic elements:" $slurm | grep -o "[0-9]*")

echo "  Total elements:    $nelemtot."
echo "  Acoustic elements: $nelemacous (included above)."

echo "> trying to find CFL <"
cfl=$(grep "Max CFL stability condition" $slurm -A2 | grep -oP "[0-9]\.[0-9]*[E|e|D|d]?\-?\+?[0-9]*" | head -2 | tail -1)
echo "  CFL:               $cfl."

echo "> trying to find last time step <"
lasttimestep=$(grep -e "Timestep" $slurm | tail -1 | grep -o "[0-9]* (" | grep -o "[0-9]*")
echo "  Last time step:    $lasttimestep."

echo "> trying to find NPROC <"
NPROC=$(grep -e "NPROC *= *[0-9]*" $parfile | grep -oe "[0-9]*")
echo "  CPUs:              $NPROC."

echo "> trying to find NSTEP_BETWEEN_OUTPUT_IMAGES <"
snapfreq=$(grep -e "NSTEP_BETWEEN_OUTPUT_IMAGES *= *[0-9]*" $parfile | grep -oe "[0-9]*" | head -1)
echo "  Snapshots saving:  $snapfreq timesteps."


echo "> trying to find number of stations <"
nbstations=$(grep "Station" $slurm | tail -1 | grep -P -o "[0-9]?[0-9]?[0-9]" | head -1)
echo "  Receivers:         $nbstations."

echo "> trying to find NSTEP_BETWEEN_OUTPUT_SEISMOS <"
synthfreq=$(grep -e "NSTEP_BETWEEN_OUTPUT_SEISMOS *= *[0-9]*" $parfile | grep -oe "[0-9]*")
echo "  Synthetics saving: $synthfreq timesteps."

dateregexep="D a t e : [0-9][0-9] - [0-9][0-9] - [0-9][0-9][0-9][0-9]"
timeregexp="T i m e : [0-9][0-9]:[0-9][0-9]:[0-9][0-9]"
slurmregexp="[0-9]?[0-9]?[0-9]?[0-9]?[0-9][0-9][0-9][0-9][0-9]"

DATELINE1=$(grep -e "$dateregexep" $slurm | head -1)
DATELINE2=$(grep -e "$dateregexep" $slurm | tail -1) # If job did not terminate completely, this is the same as DATELINE1.
#echo $DATELINE1
#echo $DATELINE2

DATESTART=$(echo $DATELINE1 | grep -oe "$dateregexep")
#echo $DATESTART
YEARSTART=$(echo $DATESTART | grep -oe "[0-9][0-9][0-9][0-9]")
#echo $YEARSTART
MONTHSTART=$(echo $DATESTART | grep -o "[0-9][0-9] -" | tail -1 | grep -o "[0-9][0-9]")
#echo $MONTHSTART
DAYSTART=$(echo $DATESTART | grep -o "[0-9][0-9]" | head -1)
#echo $DAYSTART
TIMESTART=$(echo $DATELINE1 | grep -oe "$timeregexp" | grep -o "[0-9][0-9]:[0-9][0-9]:[0-9][0-9]")
#echo $TIMESTART
DATESTART="$YEARSTART-$MONTHSTART-$DAYSTART $TIMESTART"
STAMPSTART=$(date --utc --date "$DATESTART" +%s)
#echo "  Start date:  $DATESTART (timestamp $STAMPSTART)."

if grep -q "CANCELLED AT" $slurm; # This tries to find the cancelling message in a CALMIP slurm file.
then
  DATEEND=$(grep "CANCELLED AT" $slurm | grep -oe "[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]" | head -1)
  #STAMPEND=$(date --utc --date "$DATEEND" +%s)
  #echo "    End date:    $DATECANCELLED (timestamp $STAMPCANCEL, job cancelled)."
  #elapsed=$(($STAMPCANCEL-$STAMPSTART))
  #echo "  Elapsed:           $elapsed seconds."
elif grep -q -E "CANCELLED|FAILED" $slurm; # This tries to find the cancelling/failing message in a CEA .output file.
then
  LINEWITHDATES=$(grep -e "end =" $slurm)
  TIMEEND=$(echo $LINEWITHDATES | grep -o "[0-9][0-9]:[0-9][0-9]:[0-9][0-9]" | tail -1)
  DATEEND=$(echo $LINEWITHDATES | grep -o "[0-9][0-9]/[0-9][0-9]/[0-9][0-9][0-9][0-9]" | tail -1)
  YEAREND=$(echo $DATEEND | grep -oe "[0-9][0-9][0-9][0-9]")
  MONTHEND=$(echo $DATEEND | grep -o "[0-9][0-9]/" | tail -1 | grep -o "[0-9][0-9]")
  DAYEND=$(echo $DATEEND | grep -o "[0-9][0-9]" | head -1)
  DATEEND="$YEAREND-$MONTHEND-$DAYEND $TIMEEND"
  #STAMPEND=$(date --utc --date "$DATEEND" +%s)
  #echo "    End date:    $DATEEND."
else # No cancel message found proceed to seek SPECFEM final message.
  DATEEND=$(echo $DATELINE2 | grep -oe "$dateregexep")
  #echo $DATEEND
  YEAREND=$(echo $DATEEND | grep -oe "[0-9][0-9][0-9][0-9]")
  #echo $YEAREND
  MONTHEND=$(echo $DATEEND | grep -o "[0-9][0-9] -" | tail -1 | grep -o "[0-9][0-9]")
  #echo $MONTHEND
  DAYEND=$(echo $DATEEND | grep -o "[0-9][0-9]" | head -1)
  #echo $DAYEND
  TIMEEND=$(echo $DATELINE2 | grep -oe "$timeregexp" | grep -o "[0-9][0-9]:[0-9][0-9]:[0-9][0-9]")
  #echo $TIMEEND
  DATEEND="$YEAREND-$MONTHEND-$DAYEND $TIMEEND"
  #STAMPEND=$(date --utc --date "$DATEEND" +%s)
  #elapsed=$(($STAMPEND-$STAMPSTART))
  #echo "    End date:    $DATEEND (timestamp $STAMPEND)."
  #echo "  Run duration:      $elapsed seconds."
fi
STAMPEND=$(date --utc --date "$DATEEND" +%s)
elapsed=$(($STAMPEND-$STAMPSTART))
echo "  Run duration:      $elapsed seconds."

cflrounded=$(printf %.3f $cfl)
slurmid=$(echo $slurm | grep -Po $slurmregexp | tail -1)
foldername=$(echo $folder | grep -o "/[^/]*$" | grep -o "[^/]*")
echo "  One-liner to copy-paste in Matlab script './utils_new/estimate_run_time.m':"
echo "  RUN_RAWDATA(i,:)=[$nelemtot $nelemacous $cflrounded $lasttimestep $NPROC $snapfreq $nbstations $synthfreq $elapsed]; RUNINFO{i}={$slurmid,'$foldername'}; i=i+1;"
