ver=V1p15
outputRoot=Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018-HighMET-17Sep2018.root
dir=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/driftTube/${ver}/Data2018/v4/v4/normalized/
eval `scram runtime -sh`


hadd $outputRoot /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018A-HighMET-17Sep2018.root /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018.root /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018.root /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018D-HighMET-PromptReco.root



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
