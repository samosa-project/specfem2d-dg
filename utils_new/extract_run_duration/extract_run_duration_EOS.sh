#!/bin/bash
# Reads a slurm file and extracts run duration.

slurm='/home/l.martire/Documents/SPECFEM/Ongoing_Work/2018_GRL/SH/OUTPUT_FILES_668888_stopped_12kit/slurm-668888.out'
#slurm='/home/l.martire/Documents/SPECFEM/Ongoing_Work/2018_GRL/SH/SH_soft_final/OUTPUT_FILES_593959/slurm-593959.out'

echo "Analysing SLURM file '$slurm'."

dateregexep="D a t e : [0-9][0-9] - [0-9][0-9] - [0-9][0-9][0-9][0-9]"
timeregexp="T i m e : [0-9][0-9]:[0-9][0-9]:[0-9][0-9]"

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
echo "  Start date:  $DATESTART (timestamp $STAMPSTART)."

if grep -q "CANCELLED AT" $slurm;
then
  DATECANCELLED=$(grep "CANCELLED AT" $slurm | grep -oe "[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]" | head -1)
  STAMPCANCEL=$(date --utc --date "$DATECANCELLED" +%s)
  echo "  Job cancelled."
  echo "  Cancel date: $DATECANCELLED (timestamp $STAMPCANCEL)."
  echo "  Elapsed seconds: $(($STAMPCANCEL-$STAMPSTART))."
else
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
  STAMPEND=$(date --utc --date "$DATEEND" +%s)
  echo "  End date:    $DATEEND (timestamp $STAMPEND)."
  echo "  Elapsed seconds: $(($STAMPEND-$STAMPSTART))."
fi
