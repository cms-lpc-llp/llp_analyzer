#!/bin/sh

mkdir -p log
mkdir -p submit

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_CaltechT2.sh
filesPerJob=25

for sample in \
QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 


do
	echo "Sample " ${sample}
	inputfilelist=/src/RazorAnalyzer/lists/Run2/razorNtuplerV4p1/MC_Summer16_reMINIAOD/${sample}.cern.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
	maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
	analyzer=DelayedPhotonAnalyzer
	#rm submit/${sample}_Job*.jdl
	#rm log/${sample}_Job*

	for jobnumber in `seq 0 1 ${maxjob}`
	do
		
		jdl_file=submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
                #noFail=`grep YYYY log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}*.out`
		outRoot="/mnt/hadoop/store/group/phys_susy/razor/Run2Analysis/DelayedPhotonAnalysis/2016/orderByPt/jobs/${sample}_Job${jobnumber}_Of_${maxjob}.root"
		
		minimumsize=100000
                actualsize=0
                if [ -f ${outRoot} ]
                then
                        actualsize=$(wc -c <"${outRoot}")
                fi
                if [ $actualsize -ge $minimumsize ]
		then
			finished=yes
			#echo "job ${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob} finished already "
                #elif [ -z "${noFail}" ]
                #then
                #        echo "job ${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob} being processed now, be patient"
                else
                        echo "job ${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob} failed, now being resubmitted"
                        rm log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}*
                        condor_submit submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
                fi
	done
done

