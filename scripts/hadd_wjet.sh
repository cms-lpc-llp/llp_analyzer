outputRoot=WJetsToLNu_HT-70ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root
dir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p16/v1/MC_Summer16/normalized/
eval `scram runtime -sh`

hadd $outputRoot /mnt/hadoop/$dir/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root /mnt/hadoop/$dir/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root /mnt/hadoop/$dir/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root /mnt/hadoop/$dir/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root /mnt/hadoop/$dir/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root /mnt/hadoop/$dir/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root /mnt/hadoop/$dir/WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root /mnt/hadoop/$dir/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root

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
