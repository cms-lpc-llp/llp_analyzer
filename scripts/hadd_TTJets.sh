#!/bin/sh

outputRoot=TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root
dir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p16/v1/MC_Summer16/normalized/
eval `scram runtime -sh`
#hadd $outputRoot /mnt/hadoop/$dir/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root /mnt/hadoop/$dir/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root /mnt/hadoop/$dir/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root
hadd $outputRoot /mnt/hadoop/$dir/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root /mnt/hadoop/$dir/TTJets_SingleLeptFromTbar_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root

if [ -f $outputRoot ]
then
        echo "hadd done"
fi

eval `scram unsetenv -sh`
LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

gfal-copy --checksum-mode=both $outputRoot gsiftp://transfer.ultralight.org/${dir}/$outputRoot

if [ -f /mnt/hadoop/$dir/$outputRoot ]
then
	echo "copy succeed"
	rm $outputRoot
fi
