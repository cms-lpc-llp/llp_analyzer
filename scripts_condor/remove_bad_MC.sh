#!/bin/sh

mkdir -p log
mkdir -p submit

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_CaltechT2.sh
filesPerJob=20

for sample in \
GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
GMSB_L100TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau20000cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau4000cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau100cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau20000cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau4000cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau50cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau5cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau0p01cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau0p1cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau100cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau20000cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau50cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau5cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau0p01cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau0p1cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau100cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau50cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau5cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau0p01cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau0p1cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau100cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau50cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau5cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau0p1cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau0p01cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau0p1cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau800cm_13TeV-pythia8


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
		
		outRoot="/mnt/hadoop/store/group/phys_susy/razor/Run2Analysis/DelayedPhotonAnalysis/2016/jobs/${sample}_Job${jobnumber}_Of_${maxjob}.root"
		
		minimumsize=100000
                actualsize=0
                if [ -f ${outRoot} ]
                then
                        actualsize=$(wc -c <"${outRoot}")
                fi
                if [ $minimumsize -ge $actualsize ]
		then
			echo "job ${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob} failed, removing bad output file "
                fi
	done
done

