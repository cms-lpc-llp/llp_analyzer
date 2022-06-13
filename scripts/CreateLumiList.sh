#!/bin/bash

inputDir=$1
inputRunList=$2

echo ${inputDir}
echo ${inputRunList}

#For each run, make list of lumi
echo "Making list of run:lumi"
rm -v tmpLumis.txt
for r in `cat ${inputRunList}`; do
   echo "Processing Run "${r}
   for i in ${inputDir}/*.*; do
       cat ${i} | grep "^${r} " | awk '{print $1 " " $2}' | sort -n | uniq >> tmpLumis.txt
   done
done
cat tmpLumis.txt | sort -n | uniq > tmp; mv tmp tmpLumis.txt



