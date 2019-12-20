#!/bin/bash
years='Data2018'
for year in ${years}
do

	version=displacedJetMuonNtuple/V1p12/${year}/
	root_dir=/mnt/hadoop/store/group/phys_exotica/delayedjets/${version}/v1/sixie/ # for 2018ABD
	root_dir=/mnt/hadoop/store/group/phys_exotica/delayedjets/${version}/v2/sixie/MET/ # for 2018C
	list_dir=$CMSSW_BASE/src/llp_analyzer/lists/${version}
	echo $list_dir
	mkdir -p $list_dir
	#Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_MET_Run2018A-HighMET-17Sep2018-v1_v1 \
	#Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_MET_Run2018B-HighMET-17Sep2018-v1_v1 \
	#Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_MET_Run2018D-HighMET-PromptReco-v2_v1
	for sample in \
	Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018-v1_v2_v3
	do
	        echo "${list_dir}${sample}.txt"
	        rm -f ${list_dir}${sample}.txt
	        find ${root_dir}${sample} -name "*.root" >> ${list_dir}${sample}.txt
	        sed -i '/failed/d' ${list_dir}${sample}.txt
	        echo "input list created for $sample"
	done
done






















