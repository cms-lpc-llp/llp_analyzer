ver=V1p15
outputRoot=Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016-HighMET-07Aug17.root
dir=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/driftTube/${ver}/Data2016/v4/v4/normalized/
eval `scram runtime -sh`
hadd $outputRoot /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016B-HighMET-07Aug17_ver1.root /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016B-HighMET-07Aug17_ver2.root /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016C-HighMET-07Aug17.root /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016D-HighMET-07Aug17.root /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016E-HighMET-07Aug17.root /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016F-HighMET-07Aug17.root /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016G-HighMET-07Aug17.root /mnt/hadoop/$dir/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016H-HighMET-07Aug17.root



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
