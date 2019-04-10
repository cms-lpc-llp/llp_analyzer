#!/bin/sh

mkdir -p log
mkdir -p submit

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_CaltechT2.sh
filesPerJob=1

for sample in \
GMSB_L100TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau600cm_13TeV-pythia8 \
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

