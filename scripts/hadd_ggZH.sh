#!/bin/sh 
dir=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/V1p17/MC_Fall18/v2/v147/normalized/
#outputRoot=ZToMuMu_NNPDF31_13TeV-powheg_M_50_Inf_1pb_weighted.root
eval `scram runtime -sh`



outputRoot=ggZHToSS_SToBB_ms40_pl1000_137000pb_weighted.root
hadd ${outputRoot} $dir/ggZHToSS_SToBB_ZToLL_ms40_pl1000_137000pb_weighted.root $dir/ggZHToSS_SToBB_ZToNuNu_ms40_pl1000_137000pb_weighted.root $dir/ggZHToSS_SToBB_ZToQQ_ms40_pl1000_137000pb_weighted.root


if [ -f $outputRoot ]
then
	echo "hadd done"
fi

cp ${outputRoot} ${dir}/${outputRoot}

if [ -f /$dir/$outputRoot ]
then
	echo "copy succeed"
	rm $outputRoot
fi
