ver=V1p17
ver2=/v1/v157
#sample=ggH_HToSSTodddd_MH-125
#VBFH_HToSSTo4b_MH-125
for sample in \
ggH_HToSSTodddd_MH-125
do
	for year in \
	MC_Fall17
	do
		outDir=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/${year}/${ver2}/normalized/
		outputRoot=${sample}.root

		rm -f $outputRoot
		eval `scram runtime -sh`

		hadd $outputRoot ${outDir}/${sample}*ctau-*root
		if [ -f $outputRoot ]
		then
		        echo "hadd done"
		fi


		cp $outputRoot ${outDir}/$outputRoot
		if [ -f ${outDir}/$outputRoot ]
		then
			echo "copy succeed"
			rm $outputRoot
		fi
	done
done
