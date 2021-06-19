ver=V1p17
anaVer=/v5/v63/
outputRoot=(
Run2_displacedJetMuonNtupler_${ver}_Data2016_AOD_Run2016-07Aug17.root
Run2_displacedJetMuonNtupler_${ver}_Data2017_AOD_Run2017-17Nov2017.root
Run2_displacedJetMuonNtupler_${ver}_Data2018_AOD_17Sept2018_Run2018-17Sep2018.root
)
dir=(
/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/Data2016_AOD/${anaVer}/normalized/
/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/Data2017_AOD/${anaVer}/normalized/
/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/Data2018ABC_AOD/${anaVer}/normalized/
/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/Data2018D_AOD/${anaVer}/normalized/
)
eval `scram runtime -sh`


###############
# hadd by year
##############

hadd ${outputRoot[0]} /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver1.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver2.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016C-07Aug17.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016D-07Aug17.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016E-07Aug17.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016F-07Aug17.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016G-07Aug17.root /mnt/hadoop/${dir[0]}/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17.root

hadd ${outputRoot[1]} /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017C-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017D-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017F-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017G-17Nov2017.root /mnt/hadoop/${dir[1]}/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017H-17Nov2017.root

hadd ${outputRoot[2]} /mnt/hadoop/${dir[2]}/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018A-17Sep2018.root /mnt/hadoop/${dir[2]}/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018.root /mnt/hadoop/${dir[2]}/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018.root /mnt/hadoop/${dir[3]}/Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_AOD_Run2018D-PromptReco.root

eval `scram unsetenv -sh`
LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
for i in ${!outputRoot[@]}
do
	echo $i
	if [ -f ${outputRoot[$i]} ]
	then
		echo " ${outputRoot[$i]} hadd done"
	fi
	
	gfal-copy -f --checksum-mode=both  ${outputRoot[$i]} gsiftp://transfer.ultralight.org/${dir[$i]}/${outputRoot[$i]}
	
	if [ -f /mnt/hadoop/${dir[$i]}/${outputRoot[$i]} ]
	then
		echo "copy succeed"
		rm  ${outputRoot[$i]}
	fi
done
