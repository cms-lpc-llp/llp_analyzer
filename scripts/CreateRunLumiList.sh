#!/bin/bash

inputDir=$1
outputDir=$2

echo ${inputDir}
echo ${outputDir}

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


#For each run, make list of lumi
echo "Making list of run:lumi"
rm -v tmpLumis.txt
for r in `cat tmpRuns.txt`; do
   echo "Processing Run "${r}
   for i in ${inputDir}/*.*; do
       cat ${i} | grep "^${r} " | awk '{print $1 " " $2}' | sort -n | uniq >> tmpLumis.txt
   done
done
cat tmpLumis.txt | sort -n | uniq > tmp; mv tmp tmpLumis.txt



