ver=V1p17
ver2=/v1/v169/
#VBFH_HToSSTo4Tau_MH-125 \
#VBFH_HToSSTo4b_MH-125 \
for sample in \
ggH_HToSSTo4Tau_MH-125 \
ggH_HToSSTobbbb_MH-125 \
ggH_HToSSTodddd_MH-125
do
	if [[ ${sample} == "VBFH_HToSSTo4b_MH-125" || ${sample} == "ggH_HToSSTobbbb_MH-125" ]]
	then
	        mass=( 15 40 55 )
	else
	        mass=( 7 15 40 55 )
	fi
	for m in "${mass[@]}"
	do
		for ctau in \
		100 \
		1000 \
		10000 \
		100000
		do
		
			eval `scram runtime -sh`
			rm -f $outputRoot
			file_name=${sample}_MS-${m}_ctau-${ctau}_*0pb_weighted.root
			outputRoot=${sample}_MS-${m}_ctau-${ctau}_137000pb_weighted.root
			path=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/
			hadd ${outputRoot}  ${path}/${ver}/MC_Fall*${ver2}/normalized/${file_name} ${path}/${ver}/MC_Summer*${ver2}/normalized/${file_name}
			outDir=${path}/${ver}/MC_all/${ver2}/normalized/
		
			if [ -f $outputRoot ]
			then
			        echo "hadd done"
			fi
			
			eval `scram unsetenv -sh`
			LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
			#gfal-mkdir -p gsiftp://transfer-lb.ultralight.org//${outDir}
			#gfal-copy -f --checksum-mode=both $outputRoot gsiftp://transfer-lb.ultralight.org/${outDir}/$outputRoot
			mkdir -p ${outDir}	
			cp ${outputRoot} ${outDir}/${outputRoot}
			if [ -f ${outDir}/$outputRoot ]
			then
				echo "copy succeed"
				rm $outputRoot
			fi
		done
	done
done
