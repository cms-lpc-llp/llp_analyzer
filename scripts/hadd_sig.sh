ver=V1p17
ver2=/v1/v114
#sample=ggH_HToSSTodddd_MH-125
#VBFH_HToSSTo4b_MH-125
for sample in \
ggH_HToSSTo4Tau_MH-125 \
ggH_HToSSTodddd_MH-125
do
	for year in \
	MC_Fall18_FullGenParticles
	do
		dir1=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/${year}/${ver2}/normalized/
		outDir=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/${year}/${ver2}/normalized/
		outputRoot=${sample}.root

		rm -f $outputRoot
		eval `scram runtime -sh`
		#hadd $outputRoot /mnt/hadoop/${dir1}/Z*root /mnt/hadoop/${dir1}/W*root

		hadd $outputRoot /mnt/hadoop/${dir1}/${sample}*ctau-*root
		#if [ ${sample} == 'WHToSS_SToBB' ]
		#then
		#	hadd $outputRoot /mnt/hadoop/${dir1}/W*HToSS_SToBB*pl100*root
		#else
		#	hadd $outputRoot /mnt/hadoop/${dir1}/${sample}*pl100*root
		#fi
		if [ -f $outputRoot ]
		then
		        echo "hadd done"
		fi

		eval `scram unsetenv -sh`
		LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

		gfal-copy -f --checksum-mode=both $outputRoot gsiftp://transfer.ultralight.org/${outDir}/$outputRoot

		if [ -f /mnt/hadoop/${outDir}/$outputRoot ]
		then
			echo "copy succeed"
			rm $outputRoot
		fi
	done
done
