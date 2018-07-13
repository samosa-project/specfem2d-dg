#!/bin/bash
# Backups simulation configuration files (except external_bottom_forcing files, since they are huge).

SOURCE="."
TARGET="/users/p1404/martire/save/specfem-dg-master/EXAMPLES"

echo ""
echo "# Quotas before: ###############"
quota
echo "################################"
echo ""

echo "Saving."
for file in "Par_file" "SOURCE" "atmospheric_model.dat" "parfile_input" "interfaces_input" "source_input" "slurm"
do
  find $SOURCE -name $file -exec cp -bru --parents \{\} $TARGET \;
done

echo ""
echo "# Quotas after: ################"
quota
echo "################################"
echo ""
