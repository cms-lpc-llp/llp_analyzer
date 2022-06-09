#!/bin/bash

ntupleMapDir=$1
muonHitsMapByLumiDir=$2
outputDir=$3

echo ${ntupleMapDir}
echo ${muonHitsMapByLumiDir}

for i in ${ntupleMapDir}/*.*; do
    echo $i
    index=`echo "$i"| awk -F"RunDataRunEventIndexing_option0_" '{print $2}' | sed "s/.root//"`
    rm -f tmpFileList.txt
    for runlumi in `cat $i | awk '{print $1","$2}' | sort -n | uniq `; do
	#echo ${runlumi}
	run=`echo ${runlumi}|awk -F"," '{print $1}'`
	lumi=`echo ${runlumi}|awk -F"," '{print $2}'`
	#echo "${run} ${lumi}"
	if test -f "${muonHitsMapByLumiDir}/RunEventMap_${run}_${lumi}.txt"; then
	    cat ${muonHitsMapByLumiDir}/RunEventMap_${run}_${lumi}.txt | awk '{print $5}' | sort -n | uniq 
	    cat ${muonHitsMapByLumiDir}/RunEventMap_${run}_${lumi}.txt | awk '{print $5}' | sort -n | uniq >> tmpFileList.txt
	fi	

	if test -f tmpFileList.txt; then
	    cat tmpFileList.txt | sort -n | uniq > ${outputDir}/NtupleToMuonHitsFileMap_${index}.txt
	fi

    done
    
done
