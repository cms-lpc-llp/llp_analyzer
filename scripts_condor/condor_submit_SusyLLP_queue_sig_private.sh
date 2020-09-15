#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runAnalyzerJob_SusyLLP_queue.sh
filesPerJob=5
for year in \
Fall18
do
	if [ ${year} == "Summer16" ]
        then
                tune=TuneCUETP8M1
        else
                tune=TuneCP5
        fi
        for sample in \
	GluGluH2_H2ToSSTobbbb_MH-600_MS-150_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	GluGluH2_H2ToSSTobbbb_MH-600_MS-150_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	GluGluH2_H2ToSSTobbbb_MH-600_MS-150_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	GluGluH2_H2ToSSTobbbb_MH-600_MS-150_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	GluGluH2_H2ToSSTobbbb_MH-600_MS-150_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \
	GluGluH2_H2ToSSTobbbb_MH-600_MS-50_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	GluGluH2_H2ToSSTobbbb_MH-600_MS-50_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	GluGluH2_H2ToSSTobbbb_MH-600_MS-50_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	GluGluH2_H2ToSSTobbbb_MH-600_MS-50_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	GluGluH2_H2ToSSTobbbb_MH-600_MS-50_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \

	#GluGluH2_H2ToSSTobbbb_MH-1000_MS-150_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1000_MS-150_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1000_MS-150_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1000_MS-150_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1000_MS-150_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1000_MS-400_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1000_MS-400_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1000_MS-400_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1000_MS-400_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1000_MS-400_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1500_MS-200_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1500_MS-200_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1500_MS-200_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1500_MS-200_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1500_MS-200_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1500_MS-500_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1500_MS-500_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1500_MS-500_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1500_MS-500_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-1500_MS-500_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-2000_MS-250_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-2000_MS-250_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-2000_MS-250_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-2000_MS-250_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-2000_MS-250_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-2000_MS-600_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-2000_MS-600_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-2000_MS-600_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-2000_MS-600_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-2000_MS-600_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-200_MS-25_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-200_MS-25_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-200_MS-25_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-200_MS-25_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-200_MS-25_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-200_MS-50_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-200_MS-50_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-200_MS-50_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-200_MS-50_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-200_MS-50_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-400_MS-100_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-400_MS-100_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-400_MS-100_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-400_MS-100_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-400_MS-100_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-400_MS-50_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-400_MS-50_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-400_MS-50_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-400_MS-50_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
	#GluGluH2_H2ToSSTobbbb_MH-400_MS-50_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \

	do
		echo "Sample " ${sample}
		version=/V1p17
		output=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/${version}/v14/MC_${year}/${sample}
		#${sample%_HToSS*}_HToSSTobbbb_MH-125_MS-${mx}_ctau-${ctau_cm}_${tune}_13TeV-powheg-pythia8
		echo ${output}
		inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/${version}/MC_${year}/${sample}.txt
		nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
		maxjob=`python -c "print int($nfiles.0/$filesPerJob)+1"`
		mod=`python -c "print int($nfiles.0%$filesPerJob)"`
		if [ ${mod} -eq 0 ]
		then
		        maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
		fi
		analyzer=SusyLLP
		if [ ${year} == "Fall18" ]
		then
		        echo ${year}
			analyzerTag=CT2018_17SeptEarlyReReco
			process=ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Fall18_2018
		elif [ ${year} == "Fall17" ]
		then
		        echo ${year}
		        analyzerTag=CT2017_17Nov2017Rereco
			process=ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Fall17_2017
		elif [ ${year} == 'Summer16' ]
		then
		        echo ${year}
		        analyzerTag=CT2016_07Aug2017Rereco
			process=ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Summer16_2016
		else
		        echo "ERROR: NEED TO SET CORRECT YEAR"
		fi
		rm -f submit/${analyzer}_${year}_${sample}_Job*.jdl
		rm -f log/${analyzer}_${year}_${sample}_Job*
		
		echo "number of jobs: " ${maxjob}
		jdl_file=submit/${analyzer}_${year}_${sample}_${maxjob}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
		#echo "Arguments = ${analyzer} ${inputfilelist} no 151 ${filesPerJob} \$(ProcId) ${maxjob} ${output} ${analyzerTag} ${CMSSW_BASE} ${HOME}/" >> ${jdl_file}
		echo "Arguments = ${analyzer} ${inputfilelist} no 150 ${filesPerJob} \$(ProcId) ${maxjob} ${output} ${analyzerTag} ${CMSSW_BASE} ${HOME} ${process}" >> ${jdl_file}
		
		# option should always be 1, when running condor
		echo "Log = log/${analyzer}_${year}_${sample}_Job\$(ProcId)_Of_${maxjob}_PC.log" >> ${jdl_file}
		echo "Output = log/${analyzer}_${year}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
		echo "Error = log/${analyzer}_${year}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}
		
		#echo "Requirements=TARGET.OpSysAndVer==\"CentOS7\"" >> ${jdl_file}
		echo "Requirements=(TARGET.OpSysAndVer==\"CentOS7\" && regexp(\"blade.*\", TARGET.Machine))" >> ${jdl_file}
		
		echo "RequestMemory = 2000" >> ${jdl_file}
		echo "RequestCpus = 1" >> ${jdl_file}
		echo "RequestDisk = 4" >> ${jdl_file}
		
		echo "+RunAsOwner = True" >> ${jdl_file}
		echo "+InteractiveUser = true" >> ${jdl_file}
		#echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7\"" >> ${jdl_file}
		echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7-m202006\"" >> ${jdl_file}
		echo '+SingularityBindCVMFS = True' >> ${jdl_file}
		echo "run_as_owner = True" >> ${jdl_file}
		echo "x509userproxy = ${HOME}/x509_proxy" >> ${jdl_file}
		echo "should_transfer_files = YES" >> ${jdl_file}
		echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
		echo "Queue ${maxjob}" >> ${jdl_file}
		echo "condor_submit ${jdl_file}"
		condor_submit ${jdl_file} --batch-name ${analyzer}_${year}_${sample}
	done
done
