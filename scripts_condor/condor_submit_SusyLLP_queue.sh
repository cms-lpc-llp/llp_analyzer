#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runAnalyzerJob_SusyLLP_queue.sh
#job_script=${RazorAnalyzerDir}/scripts_condor/runAnalyzerJob_SusyLLP.sh
#filesPerJob=5  # qcd/znunu/ttjets bkg, remember to change the arguments 
filesPerJob=15  # qcd/znunu/ttjets bkg, remember to change the arguments 
#filesPerJob=10  #bkg, remember to change the arguments 
#filesPerJob=1   #sig, remember to change the arguments
#filesPerJob=50  #x1n2 sig, remember to change the arguments

#year=2016
#year=2017
year=2018

for sample in \
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018-v1_v5_v1 \

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

# JetHT 2017 aod
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017C-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017D-17Nov2017-v1_v5_v3 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017-v1_v5_v3 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v2 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v3 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v4 \


# MuonEG 2017 aod
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017F-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017D-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017C-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017-v1_v5_v1 \

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


# mc fall 18 
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

#missing qcd 500-700 one in 2018


# data 2018D
#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_Run2018D-HighMET-PromptReco-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_Run2018E-HighMET-PromptReco-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_Run2018D-HighMET-PromptReco-v2_v5_v1 \

# data 2018
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018A-HighMET-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018-v1_v5_v1 \

# mc fall 17
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

# data 2017
#Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017D-HighMET-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017B-HighMET-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017F-HighMET-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017E-HighMET-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017C-HighMET-17Nov2017-v1_v5_v1 \

# data 2016
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016D-HighMET-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016F-HighMET-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver2-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016C-HighMET-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016E-HighMET-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016G-HighMET-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016H-HighMET-07Aug17-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver1-v1_v5_v1 \

# mc summer 16
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

#n3n2-n1-hbb-hbb_mh200_pl1000_ev100000 \

#above v1p17 data 2016

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
	isZLL=0
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

	#2016 aod data MuonEG
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2016_AOD/MuonEG/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/MuonEG/${sample}.txt
	#2016 aod data SinglePhoton
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2016_AOD/SinglePhoton/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/SinglePhoton/${sample}.txt
	#2016 aod data SingleMuon
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2016_AOD/SingleMuon/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/SingleMuon/${sample}.txt
	#2016 aod data SingleElectron
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2016_AOD/SingleElectron/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/SingleElectron/${sample}.txt
	#2016 aod data JetHT
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2016_AOD/JetHT/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/JetHT/${sample}.txt

	#2017 aod data MuonEG
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v7/Data2017_AOD/MuonEG/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/MuonEG/${sample}.txt
	#2017 aod data SinglePhoton
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v6/Data2017_AOD/SinglePhoton/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/SinglePhoton/${sample}.txt
	#2017 aod data SingleMuon
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/SingleMuon/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/SingleMuon/${sample}.txt
	#2017 aod data SingleElectron
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/SingleElectron/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/SingleElectron/${sample}.txt
	#2017 aod data JetHT
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/JetHT/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/JetHT/${sample}.txt
	
	#2018ABC aod data MuonEG
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2018_AOD/MuonEG/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/MuonEG/${sample}.txt
	#2018D aod data MuonEG
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018D_AOD/MuonEG/${sample}.txt
	#2018ABC aod data SingleMuon
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2018_AOD/SingleMuon/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/SingleMuon/${sample}.txt
	#2018D aod data SingleMuon
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018D_AOD/SingleMuon/${sample}.txt
	#2018ABC aod data JetHT
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2018_AOD/JetHT/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/JetHT/${sample}.txt

	isZLL=1
	#2016 aod data SingleMuon ZLL
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2016_AOD/Zmumu/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/SingleMuon/${sample}.txt
	#2016 aod data SingleElectron  ZLL
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2016_AOD/Zee/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/SingleElectron/${sample}.txt

	#2017 aod data SingleMuon ZLL
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/Zmumu/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/SingleMuon/${sample}.txt
	# 2017 aod data SingleElectron ZLL
	#outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v8/Data2017_AOD/Zee/${sample}
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/SingleElectron/${sample}.txt

	#2018ABC aod data SingleMuon ZLL
	outputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v9/Data2018_AOD/Zmumu/${sample}
	inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/SingleMuon/${sample}.txt
	#2018D aod data SingleMuon ZLL
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018D_AOD/SingleMuon/${sample}.txt

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
	if [ ${year} == "2016" ]
	then
		analyzerTag=CT2016_07Aug2017Rereco
		if echo ${sample} | grep -q "Data"
		then
			isData=yes
			process=ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Summer16_2016
		elif echo ${sample} | grep -q "QCD"
		then
			process=QCD_HT50toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16_2016
			isData=no
		elif echo ${sample} | grep -q "WJetsToLNu"
		then
			process=WJetsToLNu_HT-70ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16_2016
			isData=no
		elif echo ${sample} | grep -q "ZJetsToNuNu"
		then
			process=ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Summer16_2016
			isData=no
		elif echo ${sample} | grep -q "TTJets"
		then
			process=TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16_2016
			isData=no
		fi
	elif [ ${year} == "2017" ]
	then
		analyzerTag=CT2017_17Nov2017Rereco
		if echo ${sample} | grep -q "Data"
		then
			isData=yes
			process=ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Fall17_2017
		elif echo ${sample} | grep -q "QCD"
		then
			process=QCD_HT50toInf_TuneCP5_13TeV-madgraph-pythia8_Fall17_2017
			isData=no
		elif echo ${sample} | grep -q "WJetsToLNu"
		then
			process=WJetsToLNu_HT-70ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_2017
			isData=no
		elif echo ${sample} | grep -q "ZJetsToNuNu"
		then
			process=ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Fall17_2017
			isData=no
		elif echo ${sample} | grep -q "TTJets"
		then
			process=TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_2017
			isData=no
		fi
	elif [ ${year} == "2018" ]
	then
		analyzerTag=CT2018_17SeptEarlyReReco
		if echo ${sample} | grep -q "Data"
		then
			isData=yes
			process=ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Fall18_2018
		elif echo ${sample} | grep -q "QCD"
		then
			process=QCD_HT50toInf_TuneCP5_13TeV-madgraphMLM-pythia8_Fall18_2018
			isData=no
		elif echo ${sample} | grep -q "WJetsToLNu"
		then
			process=WJetsToLNu_HT-70ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_Fall18_2018
			isData=no
		elif echo ${sample} | grep -q "ZJetsToNuNu"
		then
			process=ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Fall18_2018
			isData=no
		elif echo ${sample} | grep -q "TTJets"
		then
			process=TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_Fall18_2018
			isData=no
		fi
	fi
	echo ${year}
	echo ${analyzerTag}
	echo ${process}
	rm -f submit/${analyzer}_${sample}_Job*.jdl
	rm -f log/${analyzer}_${sample}_Job*

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
	#jdl_file=submit/${analyzer}_${sample}_${maxjob}_${option}.jdl
	#if [ ${isZLL} == "1" ]
	#then
	#	jdl_file=submit/${analyzer}_${sample}_${maxjob}_${op}.jdl
	#else
	#	jdl_file=submit/${analyzer}_${sample}_${maxjob}_${option}.jdl
	#fi
	echo "Universe = vanilla" > ${jdl_file}
	echo "Executable = ${job_script}" >> ${jdl_file}
	#bkg no 161
	#signal no 151
	#data yes 161 
	#MR EMU 171 PHO 181 
	echo "Arguments = ${analyzer} ${inputfilelist} ${isData} ${option} ${filesPerJob} \$(ProcId) ${maxjob} ${outputDir} ${analyzerTag} ${CMSSW_BASE} ${HOME} ${process}" >> ${jdl_file}
	#echo "Arguments = ${analyzer} ${inputfilelist} ${isData} 171 ${filesPerJob} \$(ProcId) ${maxjob} ${outputDir} ${analyzerTag} ${CMSSW_BASE} ${HOME} ${process}" >> ${jdl_file}
	#echo "Arguments = ${analyzer} ${inputfilelist} ${isData} 161 ${filesPerJob} \$(ProcId) ${maxjob} ${outputDir} ${analyzerTag} ${CMSSW_BASE} ${HOME} ${process}" >> ${jdl_file}

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
	#echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7\"" >> ${jdl_file}
	echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7-m202006\"" >> ${jdl_file}
	echo '+SingularityBindCVMFS = True' >> ${jdl_file}
	echo "run_as_owner = True" >> ${jdl_file}
	echo "x509userproxy = /storage/user/jmao/x509_proxy" >> ${jdl_file}
	echo "should_transfer_files = YES" >> ${jdl_file}
	echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
	echo "Queue ${maxjob}" >> ${jdl_file}
	echo "condor_submit ${jdl_file}"
	condor_submit ${jdl_file} -batch-name ${op}"_"${sample}
done
