ver=V1p17
anaVer=/v5/v179/
outputRoot=(
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016-HighMET-07Aug17.root
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017-HighMET-17Nov2017.root
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018-HighMET-17Sep2018.root
)
dir=(
/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/Data2016/${anaVer}/normalized/
/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/Data2017/${anaVer}/normalized/
/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/Data2018/${anaVer}/normalized/
)
dir=(
/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/V1p17/Data2016/${anaVer}/normalized/
/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/V1p17/Data2017/${anaVer}/normalized/
/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/V1p17/Data2018/${anaVer}/normalized/
)
eval `scram runtime -sh`


###############
# hadd by year
##############
hadd ${outputRoot[0]} ${dir[0]}/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016B-HighMET-07Aug17_ver1.root ${dir[0]}/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016B-HighMET-07Aug17_ver2.root ${dir[0]}/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016C-HighMET-07Aug17.root ${dir[0]}/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016D-HighMET-07Aug17.root ${dir[0]}/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016E-HighMET-07Aug17.root ${dir[0]}/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016F-HighMET-07Aug17.root ${dir[0]}/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016G-HighMET-07Aug17.root ${dir[0]}/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016H-HighMET-07Aug17.root
hadd ${outputRoot[1]} ${dir[1]}/Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017B-HighMET-17Nov2017.root ${dir[1]}/Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017C-HighMET-17Nov2017.root ${dir[1]}/Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017D-HighMET-17Nov2017.root ${dir[1]}/Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017E-HighMET-17Nov2017.root ${dir[1]}/Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017F-HighMET-17Nov2017.root


hadd ${outputRoot[2]} ${dir[2]}/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018A-HighMET-17Sep2018.root ${dir[2]}/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018.root ${dir[2]}/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018.root ${dir[2]}/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018Dv1-HighMET-PromptReco.root ${dir[2]}/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018Dv2-HighMET-PromptReco.root ${dir[2]}/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018E-HighMET-PromptReco.root 

eval `scram unsetenv -sh`
LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
for i in ${!outputRoot[@]}
do
	echo $i
	if [ -f ${outputRoot[$i]} ]
	then
		echo " ${outputRoot[$i]} hadd done"
	fi
	
	#gfal-copy -f --checksum-mode=both  ${outputRoot[$i]} gsiftp://transfer-lb.ultralight.org/${dir[$i]}/${outputRoot[$i]}
	cp ${outputRoot[$i]} ${dir[$i]}/${outputRoot[$i]}	
	if [ -f ${dir[$i]}/${outputRoot[$i]} ]
	then
		echo "copy succeed"
		rm  ${outputRoot[$i]}
	fi
done
