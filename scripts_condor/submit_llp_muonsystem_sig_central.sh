#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_llp_vH.sh
filesPerJob=20
#Fall17 \
#Fall18
for year in \
Fall18
do
	if [ ${year} == "Summer16" ]
        then
                tune=TuneCUETP8M1
        else
                tune=TuneCP5
        fi
	#VBFH_HToSSTo4Tau_MH-125_${tune}_13TeV-powheg-pythia8
        #VBFH_HToSSTo4b_MH-125_${tune}_13TeV-powheg-pythia8 \
        #ggH_HToSSTobbbb_MH-125_${tune}_13TeV-powheg-pythia8
				#ggH_HToSSTo4Tau_MH-125_${tune}_13TeV-powheg-pythia8 \
        #ggH_HToSSTodddd_MH-125_${tune}_13TeV-powheg-pythia8
        for sample in \
	ggH_HToSSTobbbb_MH-125_${tune}_13TeV-powheg-pythia8
	do
		echo "Sample " ${sample}
		version=/V1p17/MC_${year}/v1/
		output=/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${version}/v116/${sample}
		#${sample%_HToSS*}_HToSSTobbbb_MH-125_MS-${mx}_ctau-${ctau_cm}_${tune}_13TeV-powheg-pythia8
		echo ${output}
		inputfilelist=/src/llp_analyzer/lists/displacedJetMuonNtuple/${version}/sixie/${sample}.txt
		nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
		maxjob=`python -c "print int($nfiles.0/$filesPerJob)+1"`
		mod=`python -c "print int($nfiles.0%$filesPerJob)"`
		if [ ${mod} -eq 0 ]
		then
		        maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
		fi
		analyzer=llp_MuonSystem_cluster
		if [[ ${year} == "Fall18" || ${year} == "Fall18_FullGenParticles" ]]
		then
		        echo ${year}
			analyzerTag=Razor2018_17SeptEarlyReReco
		elif [ ${year} == "Fall17" ]
		then
		        echo ${year}
		        analyzerTag=Razor2017_17Nov2017Rereco
		elif [ ${year} == 'Summer16' ]
		then
		        echo ${year}
		        analyzerTag=Razor2016_07Aug2017Rereco
		else
		        echo "ERROR: NEED TO SET CORRECT YEAR"
		fi

		rm -f submit/${analyzer}_${year}_${sample}_Job*.jdl
		rm -f log/${analyzer}_${year}_${sample}_Job*

		echo "number of jobs: " ${maxjob}
		jdl_file=submit/${analyzer}_${year}_${sample}_${maxjob}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
		echo "Arguments = ${analyzer} ${inputfilelist} no 11 ${filesPerJob} \$(ProcId) ${maxjob} ${output} ${analyzerTag} ${CMSSW_BASE} ${HOME}/" >> ${jdl_file}

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
