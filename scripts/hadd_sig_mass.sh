ver=V1p17
ver2=/v1/v108
for year in \
MC_Summer16 \
MC_Fall17 \
MC_Fall18
do
	for m in \
	7 \
	15 \
	40 \
	55
	do
		sample=ggH_HToSSTo4Tau_MH-125_MS-${m}
		dir1=/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/${year}/${ver2}/normalized/
		outDir=/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/${year}/${ver2}/normalized/
		outputRoot=${sample}.root

		rm -f $outputRoot
		eval `scram runtime -sh`
		#hadd $outputRoot /mnt/hadoop/${dir1}/Z*root /mnt/hadoop/${dir1}/W*root

		hadd $outputRoot ${dir1}/${sample}*ctau-100*root

		if [ -f $outputRoot ]
		then
		        echo "hadd done"
		fi

		eval `scram unsetenv -sh`
		LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

		gfal-copy -f --checksum-mode=both $outputRoot gsiftp://transfer-lb.ultralight.org/${outDir}/$outputRoot

		if [ -f ${outDir}/$outputRoot ]
		then
			echo "copy succeed"
			rm $outputRoot
		fi
	done
done
