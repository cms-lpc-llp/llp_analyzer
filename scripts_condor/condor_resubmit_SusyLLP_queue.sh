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
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018A-17Sep2018-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018-v1_v5_v2 \
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018-v1_v5_v4 \

# JetHT 2018ABC aod
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018A-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018-v1_v5_v4 \

# MuonEG 2018ABC aod
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018A-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018-v1_v5_v1 \

# MuonEG 2018D aod
#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_AOD_Run2018D-PromptReco-v2_v5_v2 \

# SingleMuon 2018ABC aod
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018-v1_v5_v1 \

# SingleMuon 2018D aod
#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_AOD_Run2018D-PromptReco-v2_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_AOD_Run2018D-PromptReco-v2_v5_v4 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_AOD_Run2018D-PromptReco-v2_v5_v5 \


# JetHT 2016 aod
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016F-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016C-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016G-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016E-07Aug17-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver2-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver1-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17-v1_v5_v5 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver2-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016D-07Aug17-v1_v5_v2 \

# SinglePhoton 2016 aod
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016D-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016C-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016E-07Aug17-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver1-v1_v5_v3 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016G-07Aug17-v1_v5_v3 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016F-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver2-v1_v5_v4 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016C-07Aug17-v1_v5_v3 \

# SingleElectron 2016 aod
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016G-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver2-v2_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016E-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016F-07Aug17-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver1-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016C-07Aug17-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17-v1_v5_v3 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17-v1_v5_v4 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016G-07Aug17-v1_v5_v2 \

# SingleMuon 2016 aod
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016E-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016D-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016F-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016C-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver2-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016G-07Aug17-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver1-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver1-v1_v5_v1 \

# MuonEG 2016 aod
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016F-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver1-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver2-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016C-07Aug17-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016C-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016E-07Aug17-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016D-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016G-07Aug17-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016F-07Aug17-v1_v5_v3 \


# JetHT 2017 aod
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017C-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017D-17Nov2017-v1_v5_v3 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017-v1_v5_v3 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v3 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v4 \


# MuonEG 2018ABC aod
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018A-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018-v1_v5_v1 \

# MuonEG 2018D aod
#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_AOD_Run2018D-PromptReco-v2_v5_v2 \

# MuonEG 2017 aod
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017F-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017D-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017C-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017-v1_v5_v1 \

#WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8 \
#ZJetsToNuNu_HT-100To200_13TeV-madgraph \
#ZJetsToNuNu_HT-200To400_13TeV-madgraph \
#ZJetsToNuNu_HT-400To600_13TeV-madgraph \
#ZJetsToNuNu_HT-600To800_13TeV-madgraph \
#ZJetsToNuNu_HT-800To1200_13TeV-madgraph \
#ZJetsToNuNu_HT-1200To2500_13TeV-madgraph \
#ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph \
#QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8 \
#QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8 \
#QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8 \
#QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8 \
#QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8 \
#QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8 \
#QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8 \
#QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8 \
#TTJets_DiLept_genMET-80_TuneCP5_13TeV-madgraphMLM-pythia8 \
#TTJets_SingleLeptFromT_genMET-80_TuneCP5_13TeV-madgraphMLM-pythia8 \
#TTJets_SingleLeptFromTbar_genMET-80_TuneCP5_13TeV-madgraphMLM-pythia8 \

#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018A-HighMET-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018-v1_v5_v1 \

#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_Run2018D-HighMET-PromptReco-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_Run2018E-HighMET-PromptReco-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_Run2018D-HighMET-PromptReco-v2_v5_v1 \

#ZJetsToNuNu_HT-100To200_13TeV-madgraph \
#ZJetsToNuNu_HT-200To400_13TeV-madgraph \
#ZJetsToNuNu_HT-400To600_13TeV-madgraph \
#ZJetsToNuNu_HT-600To800_13TeV-madgraph \
#ZJetsToNuNu_HT-800To1200_13TeV-madgraph \
#ZJetsToNuNu_HT-1200To2500_13TeV-madgraph \
#ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph \

#WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8 \
#QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8 \
#QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8 \
#QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8 \
#QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8 \
#QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8 \
#QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8 \
#QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8 \
#QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8 \
#QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8 \
#TTJets_SingleLeptFromT_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8 \
#TTJets_SingleLeptFromTbar_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8 \
#TTJets_DiLept_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8 \


#Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017D-HighMET-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017B-HighMET-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017F-HighMET-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017E-HighMET-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017C-HighMET-17Nov2017-v1_v5_v1 \

#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016D-HighMET-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016F-HighMET-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver1-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver2-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016C-HighMET-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016E-HighMET-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016G-HighMET-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016H-HighMET-07Aug17-v1_v5_v1 \

#TTJets_SingleLeptFromT_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#TTJets_SingleLeptFromTbar_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#TTJets_DiLept_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#ZJetsToNuNu_HT-100To200_13TeV-madgraph \
#ZJetsToNuNu_HT-200To400_13TeV-madgraph \
#ZJetsToNuNu_HT-400To600_13TeV-madgraph \
#ZJetsToNuNu_HT-600To800_13TeV-madgraph \
#ZJetsToNuNu_HT-800To1200_13TeV-madgraph \
#ZJetsToNuNu_HT-1200To2500_13TeV-madgraph \
#ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph \

#WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#ZJetsToNuNu_HT-100To200_13TeV-madgraph \
#ZJetsToNuNu_HT-200To400_13TeV-madgraph \
#ZJetsToNuNu_HT-400To600_13TeV-madgraph \
#ZJetsToNuNu_HT-600To800_13TeV-madgraph \
#ZJetsToNuNu_HT-800To1200_13TeV-madgraph \
#ZJetsToNuNu_HT-1200To2500_13TeV-madgraph \
#ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph \
#TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#TTJets_SingleLeptFromTbar_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \


do
	echo "Sample " ${sample}
	#outputDir=/store/group/phys_exotica/jmao/susy_llp/llp_analyzer/${sample} 
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetTimingNtuple/V1p13/Data_2016/${sample}.txt
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p16/MC_Fall17/${sample}.txt
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p16/v6/MC_Summer16/${sample} 
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p16/MC_Summer16/${sample}.txt
	#2016 data
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v5/Data2016/${sample} 
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/${sample}.txt
	#2016 mc
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v4/MC_Summer16/${sample} 
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/MC_Summer16/${sample}.txt
	#2017 mc
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v5/MC_Fall17/${sample} 
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/MC_Fall17/${sample}.txt
	#2017 data
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v5/Data2017/${sample} 
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/${sample}.txt
	#2018 data
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v5/Data2018/${sample} 
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018/${sample}.txt
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018D/${sample}.txt
	#2018 mc
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v5/MC_Fall18/${sample} 
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/MC_Fall18/${sample}.txt

	#2017 aod data MuonEG
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v7/Data2017_AOD/MuonEG/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/MuonEG/${sample}.txt
	#2017 aod data SinglePhoton
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v6/Data2017_AOD/SinglePhoton/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/SinglePhoton/${sample}.txt
	#2017 aod data SingleMuon
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/SingleMuon/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/SingleMuon/${sample}.txt
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/Zmumu/${sample}
	#2017 aod data SingleElectron
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/SingleElectron/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/SingleElectron/${sample}.txt
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/Zee/${sample}
	#2017 aod data JetHT
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/JetHT/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/JetHT/${sample}.txt


	#2016 aod data MuonEG
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2016_AOD/MuonEG/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/MuonEG/${sample}.txt
	#2016 aod data SinglePhoton
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2016_AOD/SinglePhoton/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/SinglePhoton/${sample}.txt
	#2016 aod data SingleMuon
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2016_AOD/SingleMuon/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/SingleMuon/${sample}.txt
	#2016 aod data SingleElectron
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2016_AOD/SingleElectron/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/SingleElectron/${sample}.txt
	#2016 aod data JetHT
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2016_AOD/JetHT/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/JetHT/${sample}.txt

	#2018ABC aod data MuonEG
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v5/Data2018_AOD/MuonEG/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/MuonEG/${sample}.txt
	#2018D aod data MuonEG
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018D_AOD/MuonEG/${sample}.txt
	#2018ABC aod data SingleMuon
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2018_AOD/SingleMuon/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/SingleMuon/${sample}.txt
	#2018D aod data SingleMuon
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018D_AOD/SingleMuon/${sample}.txt
	#2018ABC aod data JetHT
	outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2018_AOD/JetHT/${sample}
	inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/JetHT/${sample}.txt

	
	#isZLL=1
	#2016 aod data SingleMuon ZLL
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2016_AOD/Zmumu/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/SingleMuon/${sample}.txt
	#2016 aod data SingleElectron  ZLL
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2016_AOD/Zee/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/SingleElectron/${sample}.txt

	#2017 aod data SingleMuon ZLL
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/Zmumu/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/SingleMuon/${sample}.txt
	# 2017 aod data SingleElectron ZLL
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/Zee/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/SingleElectron/${sample}.txt

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

	if echo ${inputfilelist} | grep -q 'MuonEG'
	then
		option=171
		op="MuonEG"
	elif echo ${inputfilelist} | grep -q 'SinglePhoton'
	then 
		option=181
		op="SinglePhoton"
	elif echo ${inputfilelist} | grep -q 'JetHT'
	then 
		option=121
		op="JetHT"
	elif echo ${inputfilelist} | grep -q 'SingleMuon'
	then 
		if [ ${isZLL} == "1" ]
		then
			option=191
			op="Zmumu"
		else
			option=131
			op="SingleMuon"
		fi
	elif echo ${inputfilelist} | grep -q 'SingleElectron'
	then 
		if [ ${isZLL} == "1" ]
		then
			option=191
			op="Zee"
		else
			option=141
			op="SingleElectron"
		fi
	else
		option=161
		op="BKG"
	fi
	echo ${op}
	echo ${option}

	jdl_file=submit/${analyzer}_${sample}_${maxjob}_${op}.jdl
	#jdl_file=submit/${analyzer}_${sample}_${maxjob}.jdl
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
			new_jdl_file=submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}_${op}.jdl
			#new_jdl_file=submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
			sed 's/$(ProcId)/'${jobnumber}'/g' ${jdl_file} > ${new_jdl_file}
			sed -i 's/Queue '${maxjob}'/Queue 1/g' ${new_jdl_file}
                        condor_submit ${new_jdl_file} -batch-name ${sample}_Job${jobnumber}
                fi
	done
done

