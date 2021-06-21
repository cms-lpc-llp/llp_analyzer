ver=V1p17
anaVer=/v5/v83/
outputRoot=(
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016-ZMu-07Aug17.root
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017-ZMu-17Nov2017.root
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018-ZMu-17Sep2018.root
)
dir=(
/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/Data2016/${anaVer}/normalized/
/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/Data2017/${anaVer}/normalized/
/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/Data2018/${anaVer}/normalized/
)
eval `scram runtime -sh`


###############
# hadd by year
##############

hadd ${outputRoot[0]} /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-ZMu-07Aug17_ver1.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-ZMu-07Aug17_ver2.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016C-ZMu-07Aug17.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016D-ZMu-07Aug17.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016E-ZMu-07Aug17.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016F-ZMu-07Aug17.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016G-ZMu-07Aug17.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016H-ZMu-07Aug17.root

hadd ${outputRoot[1]} /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017B-ZMu-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017C-ZMu-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017D-ZMu-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017E-ZMu-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017F-ZMu-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017G-ZMu-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017H-ZMu-17Nov2017.root

hadd ${outputRoot[2]} /mnt/hadoop/${dir[2]}/Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_Run2018D-ZMu-PromptReco.root /mnt/hadoop/${dir[2]}/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018B-ZMu-17Sep2018.root /mnt/hadoop/${dir[2]}/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018A-ZMu-17Sep2018.root /mnt/hadoop/${dir[2]}/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018C-ZMu-17Sep2018.root


eval `scram unsetenv -sh`
LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
for i in ${!outputRoot[@]}
do
	echo $i
	if [ -f ${outputRoot[$i]} ]
	then
		echo " ${outputRoot[$i]} hadd done"
	fi
	
	gfal-copy --checksum-mode=both  ${outputRoot[$i]} gsiftp://transfer.ultralight.org/${dir[$i]}/${outputRoot[$i]}
	
	if [ -f /mnt/hadoop/${dir[$i]}/${outputRoot[$i]} ]
	then
		echo "copy succeed"
		rm  ${outputRoot[$i]}
	fi
done
