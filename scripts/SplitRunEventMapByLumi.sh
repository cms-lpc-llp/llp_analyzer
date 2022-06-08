#!/bin/bash

inputLumisFile=$1
inputDir=$2
outputDir=$3

echo ${inputDir}
echo ${outputDir}

IFS=' '
while read -r line;
do
   read -a strarr <<< "${line}"
   echo "Processing (run,lumi) (${strarr[0]} , ${strarr[1]})"

   for i in ${inputDir}/RunDataRunEventIndexing_*.*; do
       #echo "${i}"
       cat ${i} | grep "${line}" >> ${outputDir}/RunEventMap_${strarr[0]}_${strarr[1]}.txt
   done

done < $inputLumisFile
