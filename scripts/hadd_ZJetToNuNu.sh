#!/bin/sh 
outputRoot=ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_1pb_weighted.root
dir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p16/v1/MC_Summer16/normalized/
eval `scram runtime -sh`


hadd $outputRoot /mnt/hadoop/$dir/ZJetsToNuNu_HT-100To200_13TeV-madgraph_1pb_weighted.root /mnt/hadoop/$dir/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph_1pb_weighted.root /mnt/hadoop/$dir/ZJetsToNuNu_HT-200To400_13TeV-madgraph_1pb_weighted.root /mnt/hadoop/$dir/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph_1pb_weighted.root /mnt/hadoop/$dir/ZJetsToNuNu_HT-400To600_13TeV-madgraph_1pb_weighted.root /mnt/hadoop/$dir/ZJetsToNuNu_HT-600To800_13TeV-madgraph_1pb_weighted.root /mnt/hadoop/$dir/ZJetsToNuNu_HT-800To1200_13TeV-madgraph_1pb_weighted.root
if [ -f $outputRoot ]
then
	echo "hadd done"
fi

eval `scram unsetenv -sh`
LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

gfal-copy -f $outputRoot gsiftp://transfer.ultralight.org/${dir}/$outputRoot

if [ -f /mnt/hadoop/$dir/$outputRoot ]
then
	echo "copy succeed"
	rm $outputRoot
fi
