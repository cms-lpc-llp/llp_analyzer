ver=V1p17
ver2=/v1/v4
dir1=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/MC_*/${ver2}/normalized/
outDir=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/MC_Fall18/${ver2}/normalized/
outputRoot=ZH_HToSSTobbbb_ms55_1pb_weighted.root

rm -f $outputRoot
eval `scram runtime -sh`
#hadd $outputRoot /mnt/hadoop/${dir1}/Z*root /mnt/hadoop/${dir1}/W*root

hadd $outputRoot /mnt/hadoop/${dir1}/Z*root

if [ -f $outputRoot ]
then
        echo "hadd done"
fi

eval `scram unsetenv -sh`
LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

gfal-copy --checksum-mode=both $outputRoot gsiftp://transfer.ultralight.org/${outDir}/$outputRoot

if [ -f /mnt/hadoop/${outDir}/$outputRoot ]
then
	echo "copy succeed"
	rm $outputRoot
fi
