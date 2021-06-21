ver=V1p17
ver2=/v2/v106
year=MC_Fall18
for decay in \
SToTauTau \
SToBB \
STodd
do
	for mass in \
	ms7_pl100000 \
	ms7_pl10000 \
	ms7_pl1000 \
	ms7_pl100 \
	ms15_pl100000 \
	ms15_pl10000 \
	ms15_pl1000 \
	ms15_pl100 \
	ms40_pl100000 \
	ms40_pl10000 \
	ms40_pl1000 \
	ms40_pl100 \
	ms55_pl100000 \
	ms55_pl10000 \
	ms55_pl1000 \
	ms55_pl100
	do
		sample=${decay}_${mass}
	
		dir1=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/${year}/${ver2}/normalized/
		outDir=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/${year}/${ver2}/normalized/
		outputRoot=WHToSS_${sample}_137000pb_weighted.root
		
		rm -f $outputRoot
		eval `scram runtime -sh`
		#hadd $outputRoot /mnt/hadoop/${dir1}/Z*root /mnt/hadoop/${dir1}/W*root
		
		hadd $outputRoot /mnt/hadoop/${dir1}/Wplus*${sample}_*root /mnt/hadoop/${dir1}/Wminus*${sample}_*root
		
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
