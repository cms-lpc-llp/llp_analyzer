#!/bin/sh 
dir=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/V1p17/MC_Fall17/v1/v83/normalized/
#outputRoot=ZToMuMu_NNPDF31_13TeV-powheg_M_50_Inf_1pb_weighted.root
eval `scram runtime -sh`


#hadd $outputRoot /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200_1pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300_1pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400_1pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500_1pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500_1pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000_1pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120_1pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_6000_Inf_1pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800_1pb_weighted.root

#outputRoot=ZToMuMu_NNPDF31_13TeV-powheg_M_50_800_1pb_weighted.root
#hadd $outputRoot /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200_1pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400_1pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120_1pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800_1pb_weighted.root

outputRoot=ZToMuMu_NNPDF31_13TeV-powheg_M_120_800_137000pb_weighted.root
hadd $outputRoot /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200_137000pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400_137000pb_weighted.root /mnt/hadoop/$dir/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800_137000pb_weighted.root


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
