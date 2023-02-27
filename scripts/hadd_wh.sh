ver=V1p17
ver2=/v2/v162
year=MC_Fall18
#SToTauTau \
#SToBB
for decay in \
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
	
		dir1=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/${year}/${ver2}/normalized/
		outputRoot=WHToSS_${sample}_137000pb_weighted.root
		
		rm -f $outputRoot
		eval `scram runtime -sh`
		#hadd $outputRoot /mnt/hadoop/${dir1}/Z*root /mnt/hadoop/${dir1}/W*root
		
		hadd $outputRoot ${dir1}/Wplus*${sample}_*root ${dir1}/Wminus*${sample}_*root
		
		if [ -f $outputRoot ]
		then
		        echo "hadd done"
		fi
		
		eval `scram unsetenv -sh`
		LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
		
		cp ${outputRoot} ${dir1}/${outputRoot}
		if [ -f ${dir1}/$outputRoot ]
		then
			echo "copy succeed"
			rm $outputRoot
		fi
	done
done
