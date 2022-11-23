ver=V1p17
ver2=/v5/v179/
outputRoot=Run2_displacedJetMuonNtupler_${ver}_Data2016_Data2017_Data2018-HighMET_goodLumi.root
input_path=/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/
input_path=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/
dir1=${input_path}/${ver}/Data2016/${ver2}/normalized/
dir2=${input_path}/${ver}/Data2017/${ver2}/normalized/
dir3=${input_path}/${ver}/Data2018/${ver2}/normalized/

rm -f ${outputRoot}
eval `scram runtime -sh`

hadd $outputRoot ${dir1}/Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016-HighMET-07Aug17_goodLumi.root ${dir2}/Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017-HighMET-17Nov2017_goodLumi.root ${dir3}/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018-HighMET-17Sep2018_goodLumi.root 


#/mnt/hadoop/${dir3}/Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018AB-HighMET-17Sep2018_goodLumi.root

if [ -f $outputRoot ]
then
        echo "hadd done"
fi

eval `scram unsetenv -sh`
LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

#gfal-copy -f --checksum-mode=both $outputRoot gsiftp://transfer-lb.ultralight.org/${dir3}/$outputRoot
cp ${outputRoot} ${dir3}/$outputRoot
if [ -f $dir/$outputRoot ]
then
	echo "copy succeed"
	rm ${outputRoot}
fi
