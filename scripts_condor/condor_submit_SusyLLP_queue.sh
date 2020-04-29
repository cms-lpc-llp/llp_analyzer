#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runAnalyzerJob_SusyLLP_queue.sh
#job_script=${RazorAnalyzerDir}/scripts_condor/runAnalyzerJob_SusyLLP.sh
filesPerJob=15  # qcd/znunu/ttjets bkg, remember to change the arguments 
#filesPerJob=10  #bkg, remember to change the arguments 
#filesPerJob=1   #sig, remember to change the arguments
#filesPerJob=50  #x1n2 sig, remember to change the arguments

for sample in \
TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
TTJets_SingleLeptFromTbar_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
ZJetsToNuNu_HT-100To200_13TeV-madgraph \
ZJetsToNuNu_HT-200To400_13TeV-madgraph \
ZJetsToNuNu_HT-400To600_13TeV-madgraph \
ZJetsToNuNu_HT-600To800_13TeV-madgraph \
ZJetsToNuNu_HT-800To1200_13TeV-madgraph \
ZJetsToNuNu_HT-1200To2500_13TeV-madgraph \
ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph \
QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \

#v1p17 data 2016
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver1-v1_v4_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver2-v1_v4_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016C-HighMET-07Aug17-v1_v4_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016D-HighMET-07Aug17-v1_v4_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016E-HighMET-07Aug17-v1_v4_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016F-HighMET-07Aug17-v1_v4_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016G-HighMET-07Aug17-v1_v4_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016H-HighMET-07Aug17-v1_v4_v1 \

#v1p16 muon ntuple, mc summer 16
#WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \

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

#data timing ntuple, met 2016
#Run2016C-07Aug17 \
#Run2016D-07Aug17 \
#Run2016E-07Aug17 \
#Run2016F-07Aug17 \
#Run2016G-07Aug17 \
#Run2016H-07Aug17 \
#Run2016B-07Aug17 \

do
	echo "Sample " ${sample}
	#outputDir=/store/group/phys_exotica/jmao/susy_llp/llp_analyzer/${sample} 
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v2/Data2016/${sample} 
	outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p16/v1/MC_Summer16/${sample} 
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetTimingNtuple/V1p13/Data_2016/${sample}.txt
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p16/MC_Fall17/${sample}.txt
	inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p16/MC_Summer16/${sample}.txt
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
	analyzer=SusyLLP
	analyzerTag=Razor2016_80X
	rm -f submit/${analyzer}_${sample}_Job*.jdl
	rm -f log/${analyzer}_${sample}_Job*

	jdl_file=submit/${analyzer}_${sample}_${maxjob}.jdl
	echo "Universe = vanilla" > ${jdl_file}
	echo "Executable = ${job_script}" >> ${jdl_file}
	#bkg no 161
	#signal no 151
	#data yes 161 
	echo "Arguments = ${analyzer} ${inputfilelist} yes 161 ${filesPerJob} \$(ProcId) ${maxjob} ${outputDir} ${analyzerTag} ${CMSSW_BASE} ${HOME}" >> ${jdl_file}

	# option should always be 1, when running condor
	echo "Log = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_PC.log" >> ${jdl_file}
	echo "Output = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
	echo "Error = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}
	
	#echo "+JobBatchName  = LLP_SusyAna_${sample}" >> ${jdl_file}

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
	echo "x509userproxy = /storage/user/jmao/x509_proxy" >> ${jdl_file}
	echo "should_transfer_files = YES" >> ${jdl_file}
	echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
	echo "Queue ${maxjob}" >> ${jdl_file}
	echo "condor_submit ${jdl_file}"
	condor_submit ${jdl_file} -batch-name ${sample}
done
