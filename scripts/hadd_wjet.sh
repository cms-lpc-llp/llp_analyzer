dir=/store/group/phys_exotica/delayedjets/llp_analyzer/V1p6/MC_Summer16/v3/bkg/wH/normalized/
#outputRoot=WJetsToLNu_Pt-100ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root
outputRoot=WJetsToLNu_Pt-100ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root

eval `scram runtime -sh`

hadd $outputRoot /mnt/hadoop/$dir/WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root /mnt/hadoop/$dir/WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root /mnt/hadoop/$dir/WJetsToLNu_Pt-400To600_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root /mnt/hadoop/$dir/WJetsToLNu_Pt-600ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root

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
