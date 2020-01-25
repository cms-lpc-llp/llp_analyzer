#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_llp_vH.sh
filesPerJob=50

for sample in \
Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_Run2018A-HighMET-17Sep2018 \
Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018 \
Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018 \
Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_Run2018D-HighMET-PromptReco
do
	echo "Sample " ${sample}
	#output=/storage/user/christiw/displacedJetMuonAnalyzer/V1p7/MC_Summer16/v3/bkg/wH/${sample}
	year=Data2018
	version=/V1p12/${year}/vBDT/
	output=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${version}/v1/${sample}
	echo ${output}
	inputfilelist=/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p12/${year}/${sample}.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
        maxjob=`python -c "print int($nfiles.0/$filesPerJob)+1"`
	mod=`python -c "print int($nfiles.0%$filesPerJob)"`
        if [ ${mod} -eq 0 ]
        then
                maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
        fi
	analyzer=llp_MuonSystem_bdt
	if [ ${year} == "Data2018" ]
        then
                echo "Data18 condition" 
                analyzerTag=Razor2018_17SeptEarlyReReco
        elif [ ${year} == "Data2017" ]
        then

                echo "Data17 condition" 
                analyzerTag=Razor2017_17Nov2017Rereco
        elif [ ${year} == 'Data2016' ]
        then
                echo "Data16 condition" 
                analyzerTag=Razor2016_07Aug2017Rereco
        else
                echo "ERROR: NEED TO SET CORRECT YEAR"
        fi

	rm -f submit/${analyzer}_${sample}*
	rm -f log/${analyzer}_${sample}*

	echo "job " ${maxjob}
        jdl_file=submit/${analyzer}_${sample}_${maxjob}.jdl
        echo "Universe = vanilla" > ${jdl_file}
        echo "Executable = ${job_script}" >> ${jdl_file}
        echo "Arguments = ${analyzer} ${inputfilelist} yes 1 ${filesPerJob} \$(ProcId) ${maxjob} ${output} ${analyzerTag} ${CMSSW_BASE} ${HOME}/" >> ${jdl_file}

        # option should always be 1, when running condor
        echo "Log = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_PC.log" >> ${jdl_file}
        echo "Output = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
        echo "Error = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}

        #echo "Requirements=TARGET.OpSysAndVer==\"CentOS7\"" >> ${jdl_file}
        echo "Requirements=(TARGET.OpSysAndVer==\"CentOS7\" && regexp(\"blade.*\", TARGET.Machine))" >> ${jdl_file}

        echo "RequestMemory = 2000" >> ${jdl_file}
        echo "RequestCpus = 1" >> ${jdl_file}
        echo "RequestDisk = 4" >> ${jdl_file}

        echo "+RunAsOwner = True" >> ${jdl_file}
        echo "+InteractiveUser = true" >> ${jdl_file}
        echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7\"" >> ${jdl_file}
        echo '+SingularityBindCVMFS = True' >> ${jdl_file}
        echo "run_as_owner = True" >> ${jdl_file}
        echo "x509userproxy = ${HOME}/x509_proxy" >> ${jdl_file}
        echo "should_transfer_files = YES" >> ${jdl_file}
        echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
        echo "Queue ${maxjob}" >> ${jdl_file}
        echo "condor_submit ${jdl_file}"
        condor_submit ${jdl_file} -batch-name ${sample}



done
