#!/bin/sh

mkdir -p log
mkdir -p submit

cd ../
RazorAnalyzerDir=`pwd`
cd -
inputDir=/store/group/phys_exotica/delayedjets/llp_analyzer/V1p6/MC_Summer16/v1/bkg/wH/
outputDir=${inputDir}normalized
mode=bkg
analyzer=llp_vH
job_script=${RazorAnalyzerDir}/scripts_condor/normalize.sh

for sample in \
ZJetsToNuNu_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8 \
QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1 \
ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1 \
ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1 \
ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1 \
TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
do
	echo "Sample " ${sample}

		
	jdl_file=submit/${analyzer}_normalize_${mode}_${sample}.jdl
	outRoot="/mnt/hadoop/${outputDir}/${sample}_1pb_weighted.root"
	minimumsize=1000
        actualsize=0
        if [ -f ${outRoot} ]
        then
                actualsize=$(wc -c <"${outRoot}")
        fi
        if [ $actualsize -ge $minimumsize ]
	then
		finished=yes
        else
                echo "job ${analyzer}_${sample}  failed, now being resubmitted"
                rm -f log/${analyzer}_normalize_${mode}_${sample}*
                condor_submit $jdl_file
        fi
done

