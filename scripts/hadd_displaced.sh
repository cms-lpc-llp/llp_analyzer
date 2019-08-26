dir=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/V1p7/MC_Summer16/v3/signals/wH/normalized/
eval `scram runtime -sh`

for sample in \
HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \

do
eval `scram runtime -sh`
outputRoot=WH_${sample}_1pb_weighted.root

hadd $outputRoot /mnt/hadoop/$dir/WplusH_${sample}_1pb_weighted.root /mnt/hadoop/$dir/WminusH_${sample}_1pb_weighted.root

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
	rm ${outputRoot}
	gfal-rm gsiftp://transfer.ultralight.org//$dir/WplusH_${sample}_1pb_weighted.root
	gfal-rm gsiftp://transfer.ultralight.org//$dir/WminusH_${sample}_1pb_weighted.root
fi
done

