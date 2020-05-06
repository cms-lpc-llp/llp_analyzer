#!/bin/sh

for sample in \
Run2016B-07Aug17 \
Run2016C-07Aug17 \
Run2016D-07Aug17 \
Run2016E-07Aug17 \
Run2016F-07Aug17 \
Run2016G-07Aug17 \
Run2016H-07Aug17 \

#QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \

do
	#echo "Sample" ${sample}
	inputDir=/mnt/hadoop/store/group/phys_exotica/jmao/susy_llp/llp_analyzer/${sample}/
	#echo "InputDir" ${inputDir}
	echo "hadd -f ${sample}.root ${inputDir}/*.root"
	hadd -f ${sample}.root ${inputDir}/*.root

done
