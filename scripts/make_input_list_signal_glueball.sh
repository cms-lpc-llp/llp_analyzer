#!/bin/bash
dir=/mnt/hadoop/store/group/phys_exotica/delayedjets/llpntuple/V1p4/MC_Summer16/v1/christiw
list_dir=/data/christiw/LLP/CMSSW_9_4_4/src/cms_lpc_llp/llp_analyzer/lists/llpntuple/V1p4/MC_Summer16/v1/christiw/
ev=100000

for model in ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR
do
	for mh in 125 300 500 1000 2000
	do
        	mx=$((($mh-50)/2))
        	if [ $mh -eq 125 ]
        	then
                	mx=50
        	fi
        	for pl in 500 1000 10000
        	do
			sample=${model}_mh${mh}_mx${mx}_pl${pl}_ev${ev}
			if [ -f ${list_dir}${sample}.txt ]
			then
				echo "[WARNING]: file exist for ${sample}"
			else

				find ${dir}/${sample}/ -name "*.root" >> ${list_dir}${sample}.txt
				#sed -i '/0424/!d' ${list_dir}${sample}.txt
				sed -i '/failed/d' ${list_dir}${sample}.txt

				echo "input list created for $sample"
				echo " ^ ^ "
				echo "  v  "
			fi
		done
	done
done
