#!/bin/bash
years='Data2018'
for year in ${years}
do
	version=displacedJetMuonNtuple/V1p12/${year}/v2/sixie/
	root_dir=/mnt/hadoop/store/group/phys_exotica/delayedjets/${version}
	list_dir=$CMSSW_BASE/src/llp_analyzer/lists/${version}
	echo $list_dir
	mkdir -p $list_dir
	for sample in \
	Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_Run2018C-ZElectron-17Sep2018-v1 \
	Run2_displacedJetMuonNtupler_V1p12_Data2018_EGamma_17Sept2018_Run2018D-ZElectron-PromptReco-v2 \
	Run2_displacedJetMuonNtupler_V1p12_Data2018_SingleMuon_17Sept2018_Run2018D-ZMu-PromptReco-v2
	do
	        echo "${list_dir}${sample}.txt"
	        rm -f ${list_dir}${sample}.txt
	        find ${root_dir}${sample} -name "*.root" >> ${list_dir}${sample}.txt
	        sed -i '/failed/d' ${list_dir}${sample}.txt
	        echo "input list created for $sample"
	done
done

