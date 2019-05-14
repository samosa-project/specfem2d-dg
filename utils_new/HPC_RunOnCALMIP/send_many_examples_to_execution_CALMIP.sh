#!/bin/bash
# Author:        Léo Martire.
# Description:   Sends many EXAMPLES to execution at once.
# Last modified: See file metadata.
# Usage:         N/A.
# Notes:         N/A.

nargneed=4
if [[ $# -lt $nargneed ]]; then
  echo
  echo "> This script needs at least $nargneed arguments:"
  echo "  - the number of nodes,"
  echo "  - the number of cores,"
  echo "  - the time limit (format hours:minutes:seconds), and"
  echo "  - at least one EXAMPLE folder."
  echo "> Script will now exit."
  echo
  exit -1
fi

nodes=$1
cores=$2
timelimit=$3

echo " "
echo " "
echo " "

echo "This will send to execution the following EXAMPLES:"
for var in "${@:4}"
do
  var=$(basename $var)
  echo "  $var"
done
echo "Each EXAMPLE will be ran on $nodes nodes, and each node will have $cores cores, for a total of "$(( $1*$2 ))" CPUs."
echo "Each EXAMPLE will have a time limit set as $timelimit."

echo " "

read -n 1 -s -r -p "Press any key to continue, or <CTRL+C> to stop. > "
echo ""

for var in "${@:4}"
do
  var=$(basename $var)
  echo " "
  echo "################################"
  echo "# Treating $var"
  echo "################################"
  ./send_example_to_execution_CALMIP.sh $nodes $cores $var $timelimit
done

# Show queue.
  echo " "
echo "################################"
echo "# Showing queue for user $USER (you)"
echo "################################"
squeue -lu $USER

echo " "
echo " "
echo " "
