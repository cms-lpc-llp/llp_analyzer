#!/bin/bash
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================

#./RazorRun lists/Run2/razorNtuplerV1p15/MC/50ns/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8.cern.txt RazorQCDStudy 0
##########################
# Baseline Razor Analysis
##########################

outputDir=/afs/cern.ch/work/j/jlawhorn/qcd_extrapolation/
isData=true
#listDir=`pwd`/../lists/Run2/razorNtuplerV3p3/
#listDir=`pwd`/../lists/Run2/razorNtuplerV3p2/
if [[ "${isData}" == "false" ]]; then
    listDir+="MC/"
else
    listDir+="data/"
fi

echo ${outputDir}
echo ${isData}
echo ${listDir}
for sample in HTMHT_2016B_PRv2 JetHT_2016B_PRv2 DoubleEG_2016B_PRv2 DoubleMuon_2016B_PRv2
#for sample in DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
#for sample in ZJetsToNuNu_HT-100To200_13TeV-madgraph ZJetsToNuNu_HT-200To400_13TeV-madgraph ZJetsToNuNu_HT-400To600_13TeV-madgraph ZJetsToNuNu_HT-600ToInf_13TeV-madgraph TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8  QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 
do

    echo "Sample " $sample
    inputfilelist="${listDir}${sample}.cern.txt"
    echo ${inputfilelist}
    filesPerJob=10
    nfiles=`cat $inputfilelist | wc | awk '{print $1}' `
    maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
    echo $nfiles $maxjob
    for jobnumber in `seq 0 1 $maxjob`
    do
	if [[ ! -e ${outputDir}${sample}.Job${jobnumber}Of${maxjob}.root ]]; then
	    echo "job " $jobnumber " out of " $maxjob
	    echo bsub -q 1nh -J RazorQCDStudy_${jobnumber} `pwd`/runRazorQCDJob_CERN_EOS.sh ${CMSSW_BASE} RazorQCDStudy ${inputfilelist} ${isData} 0 ${filesPerJob} ${jobnumber} ${sample}.Job${jobnumber}Of${maxjob}.root ${outputDir}
	    bsub -q 1nh -J RazorQCDStudy_${jobnumber} `pwd`/runRazorQCDJob_CERN_EOS.sh ${CMSSW_BASE} RazorQCDStudy ${inputfilelist} ${isData} 0 ${filesPerJob} ${jobnumber} ${sample}.Job${jobnumber}Of${maxjob}.root ${outputDir}
	    sleep 0.1
	fi
    done
done
