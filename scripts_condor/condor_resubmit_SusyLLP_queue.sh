#!/bin/sh

mkdir -p log
mkdir -p submit

cd ../
RazorAnalyzerDir=`pwd`
cd -
#eval `scram runtime -sh`

job_script=${RazorAnalyzerDir}/scripts_condor/runAnalyzerJob_SusyLLP_queue.sh
filesPerJob=15
#filesPerJob=10

for sample in \
WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \

#TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#TTJets_SingleLeptFromTbar_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#ZJetsToNuNu_HT-100To200_13TeV-madgraph \
#ZJetsToNuNu_HT-200To400_13TeV-madgraph \
#ZJetsToNuNu_HT-400To600_13TeV-madgraph \
#ZJetsToNuNu_HT-600To800_13TeV-madgraph \
#ZJetsToNuNu_HT-800To1200_13TeV-madgraph \
#ZJetsToNuNu_HT-1200To2500_13TeV-madgraph \
#ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph \
#QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \

#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver1-v1_v4_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver2-v1_v4_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016C-HighMET-07Aug17-v1_v4_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016D-HighMET-07Aug17-v1_v4_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016E-HighMET-07Aug17-v1_v4_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016F-HighMET-07Aug17-v1_v4_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016G-HighMET-07Aug17-v1_v4_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016H-HighMET-07Aug17-v1_v4_v1 \


do
	echo "Sample " ${sample}
	outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p16/v1/MC_Summer16/${sample} 
	inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p16/MC_Summer16/${sample}.txt
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v2/Data2016/${sample} 
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/${sample}.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
        maxjob=`python -c "print int($nfiles.0/$filesPerJob)+1"`
        mod=`python -c "print int($nfiles.0%$filesPerJob)"`
	echo "maxjob " ${maxjob}
	echo "Mod " ${mod}
        if [ ${mod} -eq 0 ]
        then
                maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
	fi
	echo "maxjob " ${maxjob}
	cntjob=`echo ${maxjob}-1 | bc`
	#echo ${cntjob}

	analyzer=SusyLLP

	jdl_file=submit/${analyzer}_${sample}_${maxjob}.jdl
	for jobnumber in `seq 0 1 ${cntjob}`
	#for jobnumber in `seq 0 1 ${maxjob}`
	do
		
		outRoot=/mnt/hadoop/${outputDir}/${sample}_Job${jobnumber}_of_${maxjob}.root
		
		minimumsize=10
                actualsize=0
                if [ -f ${outRoot} ]
                then
                        actualsize=$(wc -c <${outRoot})
                fi
                if [ $actualsize -ge $minimumsize ]
		then
			finished=yes
                else
                        echo "job ${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob} failed, now being resubmitted"
                        rm -f log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}*
			new_jdl_file=submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
			sed 's/$(ProcId)/'${jobnumber}'/g' ${jdl_file} > ${new_jdl_file}
			sed -i 's/Queue '${maxjob}'/Queue 1/g' ${new_jdl_file}
                        condor_submit ${new_jdl_file} -batch-name ${sample}_Job${jobnumber}
                fi
	done
done

