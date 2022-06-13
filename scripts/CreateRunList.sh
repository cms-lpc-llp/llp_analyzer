#!/bin/bash

inputDir=$1
echo ${inputDir}


#Make List of Runs First
echo "Making List of Runs"
rm -v tmpRuns.txt
for i in ${inputDir}/*.*; do
   echo "Processing file ${i}"
   cat ${i} | awk '{print $1}' | sort -n | uniq >> tmpRuns.txt
done
cat tmpRuns.txt | sort -n | uniq > tmp; mv tmp tmpRuns.txt
echo "Completed Making List of Runs"
echo ""
