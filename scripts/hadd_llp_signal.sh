#!/bin/bash
dir=/store/group/phys_exotica/delayedjets/llp_analyzer/V1p0/MC_Summer16/v1/signals/wH

# hadd each bin
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
                        echo "$sample"
    			if [ -f /mnt/hadoop/${dir}/${sample} ]
			then
				gfal-copy gsiftp://transfer.ultralight.org//${dir}/${sample}            	
		done
        done
done
## normalize
#
if [ -f NormalizeNtuple ]
then
	./NormalizeNtuple lists/filestonormalize/llp_signal.txt 1
fi

# move all samples to hadoop

for sample in \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
WJetsToLNu_Pt-400To600_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
WJetsToLNu_Pt-600ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
do
	if [ -f ${sample}_1pb_weighted.root ]
	then
		gfal-copy ${sample}_1pb_weighted.root ${dir}/normalized/${sample}_1pb_weighted.root
		echo "output file moved to hadoop"
		#rm ${sample}_1pb_weighted.root
	fi
done




