#!/bin/bash
root_dir=/mnt/hadoop/store/group/phys_exotica/delayedjets/${version}
list_dir=$CMSSW_BASE/src/llp_analyzer/lists/${version}
echo $list_dir
mkdir -p $list_dir

for year in \
Data2016 \
Data2017 \
Data2016
do
	version=displacedJetMuonNtuple/V1p7/Data2016/v12/sixie/
	for sample in \
	SingleElectron \
	EGamma \
	SingleMuon
	do
	        echo "${list_dir}${sample}.txt"
	        rm -f ${list_dir}${sample}.txt
	        find ${root_dir}${sample} -name "*.root" >> ${list_dir}${sample}.txt
	        sed -i '/failed/d' ${list_dir}${sample}.txt
	        echo "input list created for $sample"
	
	done
done

