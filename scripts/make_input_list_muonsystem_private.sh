#!/bin/bash

for year in \
MC_Fall18
do
	for mode in \
	ZHToSS \
	WplusHToSS \
	WminusHToSS \
	VBFHToSS \
	ttH_HToSS
	do
		for decay in \
		SToBB \
		STodd \
		SToTauTau
		do
			if [ ${decay} == "SToBB" ]
			then
				version=displacedJetMuonNtuple/V1p17/${year}/v1/sixie/
				root_dir=/storage/cms/store/group/phys_exotica/delayedjets/${version}/
        			list_dir=$CMSSW_BASE/src/llp_analyzer/lists/${version}
			else	
				version=displacedJetMuonNtuple/V1p17/${year}/v2/sixie/
				root_dir=/storage/cms/store/group/phys_exotica/delayedjets/${version}/
			fi
			list_dir=$CMSSW_BASE/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/${year}/v2/sixie/
			echo $list_dir
			mkdir -p $list_dir
			for m in \
			7 \
			15 \
			40 \
			55
			do
				for ctau in \
				100 \
				1000 \
				10000 \
				100000
				do				
					sample=${mode}_${decay}_ms${m}_pl${ctau}
	        			echo "${list_dir}${sample}.txt"
	        			rm -f ${list_dir}${sample}.txt
					sample=${sample%.txt}
	        			#find ${root_dir}${sample%_batch*}/*${sample##*ev150000_}*/ -name "*.root" -size +1000c >> ${list_dir}${sample}.txt
					find ${root_dir}${sample}/ -name "*.root" -size +1000c >> ${list_dir}${sample}.txt
					sed -i '/failed/d' ${list_dir}${sample}.txt
	        			echo "input list created for $sample"
					cat ${list_dir}${sample}.txt | wc
				done
			done
		done
	done
done

