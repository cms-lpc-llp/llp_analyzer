## Inclusive razor analysis definitions
import ROOT as rt
from array import array
import sys
import copy
import csv

#local imports
from framework import Config

#####################################
### WEIGHTS NEEDED
#####################################

razorWeightOpts = {
        "Razor2015":[],
        "Razor2016":[
                     #"reapplyNPUWeights", #remove pileup weight and multiply by new PU weight
                     #"nbjets", #apply an extra btag correction
                     #'toppt', #apply top pt weight
                     #'reapplyLepWeights', #remove lepton efficiency weights and multiply by new weight
                     #'reapplyTrigWeights', #remove lepton trigger weights and multiply by new weight
                     #'removeTrigWeights', #divide event weight by trigger weight 
                     #'removePileupWeights', #divide event weight by pileup weight
                     ], 
        }
razorWeightOpts["Razor2016G_SUSYUnblind_80X"] = razorWeightOpts["Razor2016"]
razorWeightOpts["Razor2016_MoriondRereco"] = razorWeightOpts["Razor2016"]
razorWeightOpts["Razor2016_80X"] = razorWeightOpts["Razor2016"]
razorWeightOpts["Razor2016_ICHEP_80X"] = razorWeightOpts["Razor2016"]
razorExtraWeightOpts = {
        "Razor2016G_SUSYUnblind_80X":{
            #'TTJets':['nisr'], # weight ttbar sample according to number of ISR jets
            #'TTJets1L':['nisr'], 
            #'TTJets2L':['nisr'],
                 }
        }
razorExtraWeightOpts["Razor2016_MoriondRereco"] = razorExtraWeightOpts["Razor2016G_SUSYUnblind_80X"]
razorExtraWeightOpts["Razor2015"] = {}
razorExtraWeightOpts["Razor2016_80X"] = razorExtraWeightOpts["Razor2016G_SUSYUnblind_80X"]
razorExtraWeightOpts["Razor2016_ICHEP_80X"] = razorExtraWeightOpts["Razor2016G_SUSYUnblind_80X"]
razorWeightHists = {
        "Razor2015":{},
        "Razor2016":{ 
            #"pileup": #("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data/PileupWeights/PileupReweight2016G_partial.root", "PileupReweight"),
            #"muoneff":("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightMuonSelectionEffDenominatorReco_2016G_Golden.root", "ScaleFactor_TightMuonSelectionEffDenominatorReco"),
            #"eleeff":("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightElectronSelectionEffDenominatorReco_2016G_Golden.root", "ScaleFactor_TightElectronSelectionEffDenominatorReco"),
            #"muontrig":("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleMuonTriggerEfficiency_2016G_Golden.root", "hEffEtaPt"),
            #"eletrig":("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleElectronTriggerEfficiency_2016G_Golden.root", "hEffEtaPt"),
                }
        }
razorWeightHists["Razor2016G_SUSYUnblind_80X"] = {}
razorWeightHists["Razor2016_MoriondRereco"] = {}
razorWeightHists["Razor2016_80X"] = {}
razorWeightHists["Razor2016_ICHEP_80X"] = {}

#####################################
### NTUPLES
#####################################

razorNtuples = { 
        "SingleLepton":{},
        "Dilepton":{},
        "SingleLeptonInv":{},
        "DileptonInv":{},
        "VetoLepton":{},
        "VetoTau":{},
        "SignalHadronic":{},
        "SignalMuon":{},
        "SignalElectron":{},
        "SignalLepton":{},
        "Photon":{},
        "SusySync":{},
        }

razorSamples = {
        "TTJetsSingleLepton":["Other", "DYJets", "SingleTop", "WJets", "TTJets"],
        "WJetsSingleLepton":["Other", "DYJets", "SingleTop", "TTJets", "WJets"],
        "TTJetsDilepton":["Other", "DYJets", "SingleTop", "WJets", "TTJets"],
        "VetoLepton":["Other", "ZInv", "QCD", "DYJets", "SingleTop", "WJets", "TTJets"],
        "VetoTau":["Other", "ZInv", "QCD", "DYJets", "SingleTop", "WJets", "TTJets"],
        "WJetsSingleLeptonInv":["Other","DYJets", "SingleTop", "TTJets", "WJetsInv"],
        "DYJetsDileptonInv":["Other", "SingleTop", "WJets", "TTJets", "DYJetsInv"],
        "SignalHadronic":["Other", "QCD", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets1L","TTJets2L"],
        "SignalLeptonic":["Other", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets1L", "TTJets2L"],
        "Photon":["Other", "QCD", "GJetsFrag", "GJetsInv"],
        "SusySync":["SingleTop", "WJets", "TTJets"],
        }
razorSamplesReduced = {
        "TTJetsSingleLepton":["Other", "WJets", "TTJets"],
        "WJetsSingleLepton":["Other", "TTJets", "WJets"],
        "TTJetsDilepton":["Other", "WJets", "TTJets"],
        "VetoLepton":["Other", "ZInv", "QCD", "WJets", "TTJets"],
        "VetoTau":["Other", "ZInv", "QCD", "WJets", "TTJets"],
        "WJetsSingleLeptonInv":["Other", "TTJets", "WJetsInv"],
        "DYJetsDileptonInv":["Other", "WJets", "TTJets", "DYJetsInv"],
        "SignalHadronic":["Other", "QCD", "ZInv", "WJets", "TTJets"],
        "SignalLeptonic":["Other", "ZInv", "WJets", "TTJets"],
        "Photon":["Other", "QCD", "GJetsFrag", "GJetsInv"],
        "SusySync":["SingleTop", "WJets", "TTJets"],
        }

### 2015 ntuples
dir1L2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2015/OneLeptonFull_1p23_2015Final/RazorSkim/"
dir1LInv2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2015/OneLeptonAddToMETFull_1p23_2015Final/RazorSkim/"
dir2LInv2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2015/DileptonFullAddToMET_1p23_2015Final/RazorNJets80Skim/"
dirVetoL2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2015/VetoLeptonFull_1p23_2015Final/RazorNJets80Skim/"
dirVetoTau2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2015/VetoTauFull_1p23_2015Final/"
dirSignalMC2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2015/V1p23_Background_20160108/"
dirSignalData2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/2015/V1p23_ForMoriond20160119/RazorSkim/"

prefix1L2015 = "RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim"
prefix1LInv2015 = "RunTwoRazorControlRegions_OneLeptonAddToMetFull_SingleLeptonSkim" 
prefix2LInv2015 = "RunTwoRazorControlRegions_DileptonAddToMetFull_DileptonSkim"
prefixVetoL2015 = "RunTwoRazorControlRegions_VetoLeptonFull"
prefixVetoTau2015 = "RunTwoRazorControlRegions_VetoTauFull_RazorSkim"

### 1-lepton control region
razorNtuples["SingleLepton"]["Razor2015"] = {
        "TTJets"   : dir1L2015+"/"+prefix1L2015+"_TTJets_1pb_weighted_RazorSkim.root",
        "WJets"    : dir1L2015+"/"+prefix1L2015+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
        "SingleTop": dir1L2015+"/"+prefix1L2015+"_SingleTop_1pb_weighted_RazorSkim.root",
        "DYJets"   : dir1L2015+"/"+prefix1L2015+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
        "Other"    : dir1L2015+"/"+prefix1L2015+"_Other_1pb_weighted_RazorSkim.root",
        "ZInv"     : dir1L2015+"/"+prefix1L2015+"_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",
        "QCD"      : dir1L2015+"/"+prefix1L2015+"_QCD_HTBinned_1pb_weighted_RazorSkim.root",       
        "Data"     : dir1L2015+"/"+prefix1L2015+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
        }
 
### 1-lepton invisible control region
razorNtuples["SingleLeptonInv"]["Razor2015"] = {
        "TTJets"   : dir1LInv2015+"/"+prefix1LInv2015+"_TTJets_1pb_weighted_RazorSkim.root",
        "WJetsInv" : dir1LInv2015+"/"+prefix1LInv2015+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
        "SingleTop": dir1LInv2015+"/"+prefix1LInv2015+"_SingleTop_1pb_weighted_RazorSkim.root",
        "DYJets"   : dir1LInv2015+"/"+prefix1LInv2015+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
        "Other"    : dir1LInv2015+"/"+prefix1LInv2015+"_Other_1pb_weighted_RazorSkim.root",
        "ZInv"     : dir1LInv2015+"/"+prefix1LInv2015+"_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",
        "QCD"      : dir1LInv2015+"/"+prefix1LInv2015+"_QCD_HTBinned_1pb_weighted_RazorSkim.root",
        "Data"     : dir1LInv2015+"/"+prefix1LInv2015+"_SingleLepton_Run2015D_RazorSkim_GoodLumiGolden_NoDuplicates.root"
        }

### 2-lepton invisible control region
razorNtuples["DileptonInv"]["Razor2015"] = {
        "TTJets"    : dir2LInv2015+"/"+prefix2LInv2015+"_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_RazorSkim.root",
        "WJets"     : dir2LInv2015+"/"+prefix2LInv2015+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
        "SingleTop" : dir2LInv2015+"/"+prefix2LInv2015+"_SingleTop_1pb_weighted_RazorSkim.root",
        "DYJetsInv" : dir2LInv2015+"/"+prefix2LInv2015+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
        "Other"     : dir2LInv2015+"/"+prefix2LInv2015+"_Other_1pb_weighted_RazorSkim.root",
        "Data"      : dir2LInv2015+"/"+prefix2LInv2015+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
        }

### Veto lepton control region
razorNtuples["VetoLepton"]["Razor2015"] = {
        "TTJets"   : dirVetoL2015+"/"+prefixVetoL2015+"_TTJets_1pb_weighted_RazorSkim.root",
        "WJets"    : dirVetoL2015+"/"+prefixVetoL2015+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
        "SingleTop": dirVetoL2015+"/"+prefixVetoL2015+"_SingleTop_1pb_weighted_RazorSkim.root",
        "Data"     : dirVetoL2015+"/"+prefixVetoL2015+"_HTMHT_Run2015D_GoodLumiGolden_RazorSkim.root"
        }

### Veto tau control region
razorNtuples["VetoTau"]["Razor2015"] = {
        "TTJets"   : dirVetoTau2015+"/"+prefixVetoTau2015+"_TTJets_1pb_weighted.root",
        "WJets"    : dirVetoTau2015+"/"+prefixVetoTau2015+"_WJetsToLNu_HTBinned_1pb_weighted.root",
        "SingleTop": dirVetoTau2015+"/"+prefixVetoTau2015+"_SingleTop_1pb_weighted.root",
        "DYJets"   : dirVetoTau2015+"/"+prefixVetoTau2015+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted.root",
        "Other"    : dirVetoTau2015+"/"+prefixVetoTau2015+"_Other_1pb_weighted.root",
        "Data"     : dirVetoTau2015+"/"+prefixVetoTau2015+"_HTMHT_Run2015D_GoodLumiGolden.root"
        }

### Signal region
razorNtuples["SignalHadronic"]["Razor2015"] = {
        "TTJets1L" : dirSignalMC2015+"/FullRazorInclusive_TTJets1L_1pb_weighted.root",
        "TTJets2L" : dirSignalMC2015+"/FullRazorInclusive_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root",
        "WJets"    : dirSignalMC2015+"/FullRazorInclusive_WJetsToLNu_HTBinned_1pb_weighted.root",
        "SingleTop": dirSignalMC2015+"/FullRazorInclusive_SingleTop_1pb_weighted.root",
        "Other"    : dirSignalMC2015+"/FullRazorInclusive_Other_1pb_weighted.root",
        "DYJets"   : dirSignalMC2015+"/FullRazorInclusive_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted.root",
        "ZInv"     : dirSignalMC2015+"/FullRazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted.root",
        #data-driven QCD prediction for MultiJet
        "QCD"      : dirSignalData2015+"/RazorInclusive_HTMHT_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter.root", 
        "Data"     : dirSignalData2015+"/RazorInclusive_HTMHT_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter.root"
        }

razorNtuples["SignalMuon"]["Razor2015"] = razorNtuples["SignalHadronic"]["Razor2015"].copy()

razorNtuples["SignalMuon"]["Razor2015"]["Data"] = dirSignalData2015+"RazorInclusive_SingleMuon_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter.root"

razorNtuples["SignalElectron"]["Razor2015"] = razorNtuples["SignalHadronic"]["Razor2015"].copy()

razorNtuples["SignalElectron"]["Razor2015"]["Data"] = dirSignalData2015+"RazorInclusive_SingleElectron_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter.root"

razorNtuples["SignalLepton"]["Razor2015"] = razorNtuples["SignalHadronic"]["Razor2015"].copy()

razorNtuples["SignalLepton"]["Razor2015"]["Data"] = dirSignalData2015+"RazorInclusive_SingleLepton_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter.root"

### 2016 ntuples
dirCR2016 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/"
dirSR2016 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/"
versionMC2016 = "V3p8_19Jan2017_DM"
versionData2016 = "V3p8_19Jan2017_DM"

sampleTags2016 = { "Razor2016":"",
               "Razor2016_80X":"_Razor2016_80X",
               "Razor2016G_SUSYUnblind_80X":"_Razor2016G_SUSYUnblind_80X",
               "Razor2016_MoriondRereco":"_Razor2016_MoriondRereco",
               "Razor2016_ICHEP_80X":"_Razor2016_ICHEP_80X" }

prefixes2016 = { tag:{} for tag in sampleTags2016 }
for tag, suffix in sampleTags2016.iteritems():
    prefixes2016[tag]["SingleLepton"] = "RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim"+suffix
    prefixes2016[tag]["Dilepton"] = "RunTwoRazorControlRegions_DileptonFull_DileptonSkim"+suffix
    prefixes2016[tag]["SingleLeptonInv"] = "RunTwoRazorControlRegions_OneLeptonAddToMetFull_SingleLeptonSkim"+suffix
    prefixes2016[tag]["DileptonInv"] = "RunTwoRazorControlRegions_DileptonAddToMetFull_DileptonSkim"+suffix
    prefixes2016[tag]["VetoLepton"] = "RunTwoRazorControlRegions_VetoLeptonFull"+suffix
    prefixes2016[tag]["VetoTau"] = "RunTwoRazorControlRegions_VetoTauFull_RazorSkim"+suffix
    prefixes2016[tag]["Photon"] = "RunTwoRazorControlRegions_PhotonFull"+suffix
    prefixes2016[tag]["Signal"] = "FullRazorInclusive"+suffix
    prefixes2016[tag]["SusySync2016"] = "RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim"+suffix
skimstr = ""

#on EOS
#dir1L2016 = dirCR2016+'/'+versionMC2016+'/OneLeptonFull'
#dir1LInv2016 = dirCR2016+'/'+versionMC2016+'/OneLeptonAddToMET'
#dir2LInv2016 = dirCR2016+'/'+versionMC2016+'/DileptonAddToMET'
#dir2L2016 = dirCR2016+'/'+versionMC2016+'/DileptonFull'
#dirVetoL2016 = dirCR2016+'/'+versionMC2016+'/VetoLeptonRazorSkim'
#dirVetoL2016 = dirCR2016+'/V3p6_25October2016_CustomType1MET_TestTightVeto/VetoLepton'
#dirVetoTau2016 = dirCR2016+'/'+versionMC2016+'/VetoTauRazorSkim'
#dirPhoton2016 = dirCR2016+'/'+versionMC2016+'/PhotonAddToMET'
#dirSignal2016 = dirSR2016+'/'+versionMC2016+'/Signal'
dirSusySync2016 = "eos/cms/store/group/phys_susy/razor/Run2Analysis/SusySync/2016/V3p6_25October2016_CustomType1MET/OneLeptonFull/"

#local directories
dir1L2016 = 'Backgrounds/1L'
dir2L2016 = 'Backgrounds/2L'
dir1LInv2016 = 'Backgrounds/1LInv'
dir2LInv2016 = 'Backgrounds/2LInv'
dirVetoL2016 = 'Backgrounds/VetoL'
dirVetoTau2016 = 'Backgrounds/VetoTau'
dirSignal2016 = 'Backgrounds/Signal'
dirPhoton2016 = 'Backgrounds/Photon'    

for tag in sampleTags2016:
    razorNtuples["SingleLepton"][tag] = {
            "TTJets"   : dir1L2016+"/"+prefixes2016[tag]["SingleLepton"]+"_TTJets_1pb_weighted"+skimstr+".root",
            "WJets"    : dir1L2016+"/"+prefixes2016[tag]["SingleLepton"]+"_WJets_1pb_weighted"+skimstr+".root",
            "SingleTop": dir1L2016+"/"+prefixes2016[tag]["SingleLepton"]+"_SingleTop_1pb_weighted"+skimstr+".root",
            "DYJets"   : dir1L2016+"/"+prefixes2016[tag]["SingleLepton"]+"_DYJets_1pb_weighted"+skimstr+".root",
            "Other"    : dir1L2016+"/"+prefixes2016[tag]["SingleLepton"]+"_Other_1pb_weighted"+skimstr+".root",
            "Data"     : dir1L2016+"/"+prefixes2016[tag]["SingleLepton"]+"_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root"
            }
    razorNtuples["Dilepton"][tag] = {
            "TTJets"   : dir2L2016+"/"+prefixes2016[tag]["Dilepton"]+"_TTJets_1pb_weighted"+skimstr+".root",
            "WJets"    : dir2L2016+"/"+prefixes2016[tag]["Dilepton"]+"_WJets_1pb_weighted"+skimstr+".root",
            "SingleTop": dir2L2016+"/"+prefixes2016[tag]["Dilepton"]+"_SingleTop_1pb_weighted"+skimstr+".root",
            "DYJets"   : dir2L2016+"/"+prefixes2016[tag]["Dilepton"]+"_DYJets_1pb_weighted"+skimstr+".root",
            "Other"    : dir2L2016+"/"+prefixes2016[tag]["Dilepton"]+"_Other_1pb_weighted"+skimstr+".root",
            "Data"     : dir2L2016+"/"+prefixes2016[tag]["Dilepton"]+"_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root"
            }
    razorNtuples["SingleLeptonInv"][tag] = {
            "TTJets"   : dir1LInv2016+"/"+prefixes2016[tag]["SingleLeptonInv"]+"_TTJets_1pb_weighted"+skimstr+".root",
            "WJetsInv"    : dir1LInv2016+"/"+prefixes2016[tag]["SingleLeptonInv"]+"_WJets_1pb_weighted"+skimstr+".root",
            "SingleTop": dir1LInv2016+"/"+prefixes2016[tag]["SingleLeptonInv"]+"_SingleTop_1pb_weighted"+skimstr+".root",
            "DYJets"   : dir1LInv2016+"/"+prefixes2016[tag]["SingleLeptonInv"]+"_DYJets_1pb_weighted"+skimstr+".root",
            "Other"    : dir1LInv2016+"/"+prefixes2016[tag]["SingleLeptonInv"]+"_Other_1pb_weighted"+skimstr+".root",
            "Data"     : dir1LInv2016+"/"+prefixes2016[tag]["SingleLeptonInv"]+"_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root"
            }
    razorNtuples["DileptonInv"][tag] = {
            "TTJets"   : dir2LInv2016+"/"+prefixes2016[tag]["DileptonInv"]+"_TTJets_1pb_weighted"+skimstr+".root",
            "WJets"    : dir2LInv2016+"/"+prefixes2016[tag]["DileptonInv"]+"_WJets_1pb_weighted"+skimstr+".root",
            "SingleTop": dir2LInv2016+"/"+prefixes2016[tag]["DileptonInv"]+"_SingleTop_1pb_weighted"+skimstr+".root",
            "DYJetsInv"   : dir2LInv2016+"/"+prefixes2016[tag]["DileptonInv"]+"_DYJets_1pb_weighted"+skimstr+".root",
            "Other"    : dir2LInv2016+"/"+prefixes2016[tag]["DileptonInv"]+"_Other_1pb_weighted"+skimstr+".root",
            "Data"     : dir2LInv2016+"/"+prefixes2016[tag]["DileptonInv"]+"_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root"
            }
    razorNtuples["VetoLepton"][tag] = {
            "TTJets"   : dirVetoL2016+"/"+prefixes2016[tag]["VetoLepton"]+"_TTJets_1pb_weighted"+skimstr+".root",
            "WJets"    : dirVetoL2016+"/"+prefixes2016[tag]["VetoLepton"]+"_WJets_1pb_weighted"+skimstr+".root",
            "SingleTop": dirVetoL2016+"/"+prefixes2016[tag]["VetoLepton"]+"_SingleTop_1pb_weighted"+skimstr+".root",
            "DYJets"   : dirVetoL2016+"/"+prefixes2016[tag]["VetoLepton"]+"_DYJets_1pb_weighted"+skimstr+".root",
            "Other"    : dirVetoL2016+"/"+prefixes2016[tag]["VetoLepton"]+"_Other_1pb_weighted"+skimstr+".root",
            "ZInv"     : dirVetoL2016+"/"+prefixes2016[tag]["VetoLepton"]+"_ZInv_1pb_weighted"+skimstr+".root",
            "QCD"      : dirVetoL2016+"/"+prefixes2016[tag]["VetoLepton"]+"_QCD_1pb_weighted"+skimstr+".root",       
            "Data"     : dirVetoL2016+"/"+prefixes2016[tag]["VetoLepton"]+"_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root"
            }
    razorNtuples["VetoTau"][tag] = {
            "TTJets"   : dirVetoTau2016+"/"+prefixes2016[tag]["VetoTau"]+"_TTJets_1pb_weighted"+skimstr+".root",
            "WJets"    : dirVetoTau2016+"/"+prefixes2016[tag]["VetoTau"]+"_WJets_1pb_weighted"+skimstr+".root",
            "SingleTop": dirVetoTau2016+"/"+prefixes2016[tag]["VetoTau"]+"_SingleTop_1pb_weighted"+skimstr+".root",
            "DYJets"   : dirVetoTau2016+"/"+prefixes2016[tag]["VetoTau"]+"_DYJets_1pb_weighted"+skimstr+".root",
            "Other"    : dirVetoTau2016+"/"+prefixes2016[tag]["VetoTau"]+"_Other_1pb_weighted"+skimstr+".root",
            "ZInv"     : dirVetoTau2016+"/"+prefixes2016[tag]["VetoTau"]+"_ZInv_1pb_weighted"+skimstr+".root",
            "QCD"      : dirVetoTau2016+"/"+prefixes2016[tag]["VetoTau"]+"_QCD_1pb_weighted"+skimstr+".root",       
            "Data"     : dirVetoTau2016+"/"+prefixes2016[tag]["VetoTau"]+"_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root"
            }
    razorNtuples["Photon"][tag] = {
            "GJetsInv"    : dirPhoton2016+"/"+prefixes2016[tag]["Photon"]+"_GJets_1pb_weighted"+skimstr+".root",
            "GJetsFrag": dirPhoton2016+"/"+prefixes2016[tag]["Photon"]+"_QCD_1pb_weighted"+skimstr+".root",
            "Other"    : dirPhoton2016+"/"+prefixes2016[tag]["Photon"]+"_Other_1pb_weighted"+skimstr+".root",
            #QCD predicted using data driven method
            "QCD"      : dirPhoton2016+"/"+prefixes2016[tag]["Photon"]+"_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root",       
            "Data"     : dirPhoton2016+"/"+prefixes2016[tag]["Photon"]+"_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root",       
            }
    razorNtuples["SignalHadronic"][tag] = {
            "TTJets1L"   : dirSignal2016+"/"+prefixes2016[tag]["Signal"]+"_TTJets1L_1pb_weighted"+skimstr+".root",
            "TTJets2L"   : dirSignal2016+"/"+prefixes2016[tag]["Signal"]+"_TTJets2L_1pb_weighted"+skimstr+".root",
            "WJets"    : dirSignal2016+"/"+prefixes2016[tag]["Signal"]+"_WJets_1pb_weighted"+skimstr+".root",
            "SingleTop": dirSignal2016+"/"+prefixes2016[tag]["Signal"]+"_SingleTop_1pb_weighted"+skimstr+".root",
            "DYJets"   : dirSignal2016+"/"+prefixes2016[tag]["Signal"]+"_DYJets_1pb_weighted"+skimstr+".root",
            "Other"    : dirSignal2016+"/"+prefixes2016[tag]["Signal"]+"_Other_1pb_weighted"+skimstr+".root",
            "ZInv"     : dirSignal2016+"/"+prefixes2016[tag]["Signal"]+"_ZInv_1pb_weighted"+skimstr+".root",
            #QCD predicted using data driven method
            "QCD"      : dirSignal2016+"/"+prefixes2016[tag]["Signal"]+"_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root",
            "Data"     : dirSignal2016+"/"+prefixes2016[tag]["Signal"]+"_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root"
            }
    razorNtuples["SignalLepton"][tag] = razorNtuples["SignalHadronic"][tag].copy()
    razorNtuples["SignalMuon"][tag] = razorNtuples["SignalHadronic"][tag].copy()
    razorNtuples["SignalElectron"][tag] = razorNtuples["SignalHadronic"][tag].copy()
    razorNtuples["SusySync"][tag] = {
            "TTJets"    : dirSusySync2016+"/"+prefixes2016[tag]["SusySync2016"]+"_TTJets_1pb_weighted_HTMETSkim.root",
            "WJets"     : dirSusySync2016+"/"+prefixes2016[tag]["SusySync2016"]+"_WJets_1pb_weighted_HTMETSkim.root",
            "SingleTop" : dirSusySync2016+"/"+prefixes2016[tag]["SusySync2016"]+"_SingleTop_1pb_weighted_HTMETSkim.root",
            "Data"      : dirSusySync2016+"/"+prefixes2016[tag]["SusySync2016"]+"_HTMHT_2016BCD_PRv2_NoDuplicates_GoodLumiGolden.root"
            }

    #update version number for data, in case different
    for name,files in razorNtuples.iteritems():
        if "Razor2016" in files:
            files[tag]["Data"] = files[tag]["Data"].replace(
                    '/'+versionMC2016+'/','/'+versionData2016+'/')


#####################################
### FITS
#####################################

razorFitDirs = { 
        "Razor2016G_SUSYUnblind_80X":"/afs/cern.ch/work/j/jlawhorn/public/Razor_Moriond2017/CMSSW_7_1_5/src/RazorAnalyzer/fits_2017_01_09/ReReco2016_02Jan/"
        }
razorFitFiles = { tag:{} for tag in razorFitDirs }
for tag,path in razorFitDirs.iteritems():
    for box in ["MultiJet","DiJet","LeptonMultiJet","LeptonJet"]:
        razorFitFiles[tag][box] = "%s/%s/Sideband/toys_Bayes_noStat_%s.root"%(path,box,box)


#####################################
### TRIGGER
#####################################

triggerNames = {}
for ttype in ["Razor", "SingleLepton", "Dilepton", "Photon", "SusySync"]: triggerNames[ttype] = { "Data":{}, "MC":{} }

# 2015 triggers
# Data
triggerNames["Razor"]["Data"]["Razor2015"] = [
        'HLT_RsqMR240_Rsq0p09_MR200',
        'HLT_RsqMR240_Rsq0p09_MR200_4jet',
        'HLT_RsqMR270_Rsq0p09_MR200',
        'HLT_RsqMR270_Rsq0p09_MR200_4jet',
        'HLT_RsqMR260_Rsq0p09_MR200',
        'HLT_RsqMR260_Rsq0p09_MR200_4jet',
        'HLT_RsqMR300_Rsq0p09_MR200',
        'HLT_RsqMR300_Rsq0p09_MR200_4jet',
        'HLT_Rsq0p25',
        'HLT_Rsq0p30',
        'HLT_Rsq0p36',
        ]
triggerNames["SingleLepton"]["Data"]["Razor2015"] = [
        'HLT_Mu50',
        'HLT_IsoMu20',
        'HLT_IsoMu27',
        'HLT_IsoTkMu20',
        'HLT_IsoTkMu27',
        'HLT_Ele23_WPLoose_Gsf',
        'HLT_Ele27_WPLoose_Gsf',
        'HLT_Ele27_eta2p1_WPLoose_Gsf',
        'HLT_Ele27_eta2p1_WPTight_Gsf',
        'HLT_Ele32_eta2p1_WPLoose_Gsf',
        'HLT_Ele32_eta2p1_WPTight_Gsf',
        'HLT_Ele105_CaloIdVT_GsfTrkIdT',
        'HLT_Ele115_CaloIdVT_GsfTrkIdT',
        ]
triggerNames["Dilepton"]["Data"]["Razor2015"] = [
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
        'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
        'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
        'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL',
        'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL',
        ]
# MC
triggerNames["Razor"]["MC"]["Razor2015"] = copy.copy(triggerNames["Razor"]["Data"]["Razor2015"])
triggerNames["SingleLepton"]["MC"]["Razor2015"] = [
        'HLT_Mu50',
        'HLT_IsoMu20',
        'HLT_IsoMu27',
        'HLT_IsoTkMu20',
        'HLT_IsoTkMu27',
        'HLT_Ele23_WP75_Gsf',
        'HLT_Ele27_WP85_Gsf',
        'HLT_Ele27_eta2p1_WP75_Gsf',
        'HLT_Ele32_eta2p1_WP75_Gsf',
        'HLT_Ele105_CaloIdVT_GsfTrkIdT',
        'HLT_Ele115_CaloIdVT_GsfTrkIdT',
        'HLT_Ele22_eta2p1_WP75_Gsf',
        ]
triggerNames["Dilepton"]["MC"]["Razor2015"] = copy.copy(triggerNames["Dilepton"]["Data"]["Razor2015"])

# 2016 trigger
# Data
triggerNames["Razor"]["Data"]["Razor2016"] = copy.copy(triggerNames["Razor"]["Data"]["Razor2015"])
triggerNames["SingleLepton"]["Data"]["Razor2016"] = [
        "HLT_IsoMu20", 
        "HLT_IsoTkMu20", 
        "HLT_IsoMu22", 
        "HLT_IsoTkMu22", 
        "HLT_IsoMu24", 
        "HLT_IsoTkMu24", 
        "HLT_IsoMu27", 
        "HLT_IsoTkMu27", 
        "HLT_Mu50",
        #"HLT_TkMu50",
        "HLT_Ele23_WPLoose_Gsf",
        "HLT_Ele27_WPLoose_Gsf",
        "HLT_Ele27_WPTight_Gsf",
        "HLT_Ele27_eta2p1_WPLoose_Gsf",
        "HLT_Ele27_eta2p1_WPTight_Gsf",
        "HLT_Ele32_eta2p1_WPLoose_Gsf",
        "HLT_Ele32_eta2p1_WPTight_Gsf",
        "HLT_Ele105_CaloIdVT_GsfTrkIdT",
        "HLT_Ele115_CaloIdVT_GsfTrkIdT",
        ]
triggerNames["Dilepton"]["Data"]["Razor2016"] = copy.copy(triggerNames["Dilepton"]["Data"]["Razor2015"])
triggerNames["Photon"]["Data"]["Razor2016"] = [
        'HLT_Photon165_HE10',
        ]
triggerNames["SusySync"]["Data"]["Razor2016"] = [
        'HLT_PFMET100_PFMHT100_IDTight',
        'HLT_PFMETNoMu100_PFMHTNoMu100_IDTight',
        ]
#MC (samples have broken HLT menu)
triggerNames["Razor"]["MC"]["Razor2016"] = [ 'PASS' ]
triggerNames["SingleLepton"]["MC"]["Razor2016"] = [ 'PASS' ]
triggerNames["Dilepton"]["MC"]["Razor2016"] = [ 'PASS' ]
triggerNames["Photon"]["MC"]["Razor2016"] = [ 'PASS' ]
triggerNames["SusySync"]["MC"]["Razor2016"] = [ 'PASS' ]

#Use 2016 triggers for all 2016 eras
for tag in ["Razor2016_80X", "Razor2016_ICHEP_80X", "Razor2016G_SUSYUnblind_80X", "Razor2016_MoriondRereco"]:
    for proc in triggerNames:
        for datamc in ["Data","MC"]:
            triggerNames[proc][datamc][tag] = triggerNames[proc][datamc]["Razor2016"]

class TriggerUtils(object):
    """Stores and retrieves trigger names"""
    def __init__(self, trigType, tag, isData=False):
        self.trigType = trigType
        self.tag = tag
        if isData:
            self.isDataTag = "Data"
        else:
            self.isDataTag = "MC"

        # get trigger path name/number correspondence
        self.triggerDict = {}
        if tag == 'Razor2015':
            trigFile = 'data/RazorHLTPathnames_2015.dat'
        else:
            trigFile = 'data/RazorHLTPathnames_2016.dat'
        with open(trigFile) as f:
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                self.triggerDict[row[-1]] = int(row[0]) 

        # get trigger names
        if self.trigType in triggerNames and self.isDataTag in triggerNames[self.trigType] and self.tag in triggerNames[self.trigType][self.isDataTag]:
            self.names = triggerNames[self.trigType][self.isDataTag][self.tag]
        else:
            self.badInit(trigType,tag)

        # get trigger numbers
        self.getTriggerNums()

    def badInit(self, trigType, tag):
        print "Error in triggerUtils: trigger type/tag combination '",trigType,tag,"' not recognized!"
        self.names = None

    def getTriggerNums(self):
        self.nums = []
        for name in self.names:
            if name == 'PASS': # indicates that no trigger cuts should be applied
                self.nums = [-1] 
                return
            elif name in self.triggerDict:
                self.nums.append( self.triggerDict[name] )
            else:
                print "Warning in triggerUtils.getTriggerNums: trigger name",name,"not found!"

    def appendTriggerCuts(self, cutsString):
        """Append a string of the form "(HLTDecision[t1] || HLTDecision[t2] || ... || HLTDecision[tN]) && " 
           to the provided cut string, where t1...tN are the desired trigger numbers"""
        if -1 in self.nums: # trigger requirement not applied
            return cutsString
        return '('+(' || '.join(['HLTDecision['+str(n)+']' for n in self.nums]))+") && "+cutsString

#####################################
### CUTS
#####################################

razorCuts = {}
razorExtraCuts = {} #per-process cuts

def appendBoxCuts(cuts, boxNums):
    """Append a string of the form "(box == b1 || box == b2 || ... || box == bN) && " to the provided cut string, where b1...bN are the desired box numbers"""
    return '('+(' || '.join(['box == '+str(n) for n in boxNums])) + ") && " + cuts

recommendedNoiseFilters = ["Flag_HBHENoiseFilter","Flag_HBHEIsoNoiseFilter","Flag_goodVertices",
        "Flag_eeBadScFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_CSCTightHaloFilter",
        "Flag_badChargedCandidateFilter","Flag_badMuonFilter"
        ]
def appendNoiseFilters(cuts, tree=None):
    ret = copy.copy(cuts)
    for bit in recommendedNoiseFilters:
        if (tree is None) or hasattr(tree, bit):
            ret += " && " + bit
    return ret

# 1-lepton control region
razorCuts["SingleLepton"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && ((abs(lep1Type) == 11 && lep1.Pt() > 30) || (abs(lep1Type) == 13 && lep1.Pt() > 25)) && MET > 30 && lep1MT > 30 && lep1MT < 100 && MR > 150 && Rsq > 0.4"

# 1-lepton invisible control region
razorCuts["SingleLeptonInv"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && ((abs(lep1Type) == 11 && lep1.Pt() > 30) || (abs(lep1Type) == 13 && lep1.Pt() > 25)) && MET > 30 && lep1MT > 30 && lep1MT < 100 && NJets80_NoW >= 2 && MR_NoW > 150 && Rsq_NoW > 0.4 && ( weight < 0.01 || weight == 1)"

# TTJets 2-lepton control region
razorCuts["TTJetsDilepton"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && (abs(lep2Type) == 11 || abs(lep2Type) == 13) && lep1PassTight && lep2PassTight && lep1.Pt() > 30 && lep2.Pt() > 30 && mll > 20 && ((abs(lep1Type) != abs(lep2Type)) || (mll < 76 || mll > 106)) && NBJetsMedium > 0 && MET > 40 && MR > 150 && Rsq > 0.4 && Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter"

# DYJets 2-lepton invisible control region
razorCuts["DYJetsDileptonInv"] = "((abs(lep1Type) == 11 && abs(lep2Type) == 11) || (abs(lep1Type) == 13 && abs(lep2Type) == 13)) && lep1.Pt() > 30 && lep2.Pt() > 20 && recoZmass > 80 && recoZmass < 110 && NBJetsMedium == 0 && NJets80_NoZ >= 2 && MR_NoZ > 150 && Rsq_NoZ > 0.4"

# Veto lepton control region
razorCuts["VetoLepton"]  = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassVeto && lep1.Pt() > 5 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 150 && Rsq > 0.4"

razorCuts["VetoElectron"] = "(abs(lep1Type) == 11) && lep1PassVeto && lep1.Pt() > 5 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 150 && Rsq > 0.4"

razorCuts["VetoMuon"] = "(abs(lep1Type) == 13) && lep1PassVeto && lep1.Pt() > 5 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 150 && Rsq > 0.4"

razorCuts["VetoTau"] = "(abs(lep1Type) == 15) && lep1PassLoose && lep1.Pt() > 20 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 150 && Rsq > 0.4"

# Photon control region
razorCuts["Photon"] = "pho1.Pt() > 185 && pho1_chargediso < 2.5 && ((abs(pho1.Eta()) < 1.5 && pho1_sigmaietaieta < 0.01031) || (abs(pho1.Eta()) >= 1.5 && pho1_sigmaietaieta < 0.03013)) && NJets80_NoPho >= 2 && MR_NoPho > 150 && Rsq_NoPho > 0.4 && Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter"
razorExtraCuts["Photon"] = {
        "GJetsInv":"!(abs(pho1_motherID) > 5 && pho1_motherID != 21 && pho1_motherID != 2212) && minDRGenPhotonToParton >= 0.4",
        "GJetsFrag":"!(abs(pho1_motherID) > 5 && pho1_motherID != 21 && pho1_motherID != 2212) && minDRGenPhotonToParton < 0.4",
        }

# For October 2016 SUSY synchronization exercise
razorCuts["SusySync"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && lep1.Pt() > 25 && abs(lep1.Eta()) < 2.4 && MET > 250 && HT > 250 && NJets40 >= 2 && Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter"

# Hadronic boxes with inverted lepton veto
cutsMultiJetForVeto = "(box == 11 || box == 12 || box == 21) && MR > 150.000000 && Rsq > 0.40000 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"
cutsDiJetForVeto = "(box == 14) && MR > 150.000000 && Rsq > 0.40000 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"

razorCuts["MultiJetForVetoLepton"] = cutsMultiJetForVeto+" && ( abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13 ) && leadingGenLeptonPt > 5"

razorCuts["MultiJetForVetoTau"] = cutsMultiJetForVeto+" && abs(leadingGenLeptonType) == 15 && leadingGenLeptonPt > 20"

razorCuts["DiJetForVetoLepton"] = cutsDiJetForVeto+" && ( abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13 ) && leadingGenLeptonPt > 5"

razorCuts["DiJetForVetoTau"] = cutsDiJetForVeto+" && abs(leadingGenLeptonType) == 15 && leadingGenLeptonPt > 20"

# Signal boxes
razorBoxes = {
        "MuEle" : [0], 
        "MuMu" : [1],
        "EleEle" : [2],
        "MuSixJet" : [3],
        "MuFourJet" : [4],
        "MuJet" : [5],
        "EleSixJet" : [6],
        "EleFourJet" : [7],
        "EleJet" : [8],
        "LeptonSixJet" : [3,6],
        "LeptonFourJet" : [4,7],
        "LeptonJet" : [5,8],
        "LooseLeptonSixJet" : [9],
        "LooseLeptonFourJet" : [10],
        "LooseLeptonDiJet" : [13],
        "SixJet" : [11],
        "FourJet" : [12],
        "DiJet" : [14],	  
        "MuMultiJet" : [3,4,18],
        "EleMultiJet" : [6,7,19],
        "LeptonMultiJet" : [3,4,18,6,7,19],
        "LooseLeptonMultiJet" : [9,10,20],
        "MultiJet" : [11,12,21],
        "FourToSixJet" : [11,12,21],
        "SevenJet" : [11,12,21],
        }

hadronicRazorBoxes = ["DiJet", "FourJet", "SixJet", "MultiJet", "FourToSixJet", "SevenJet"]
looseLeptonRazorBoxes = ["LooseLeptonDiJet", "LooseLeptonFourJet", "LooseLeptonSixJet", "LooseLeptonMultiJet"]
leptonicRazorBoxes = ["MuJet", "MuFourJet", "MuSixJet", "MuMultiJet","EleJet", "EleFourJet", "EleSixJet", "EleMultiJet", "LeptonJet", "LeptonFourJet", "LeptonSixJet", "LeptonMultiJet"]
electronRazorBoxes = ["EleJet", "EleFourJet", "EleSixJet", "EleMultiJet"]
muonRazorBoxes = ["MuJet", "MuFourJet", "MuSixJet", "MuMultiJet"]
leptonRazorBoxes = ["LeptonJet", "LeptonFourJet", "LeptonSixJet", "LeptonMultiJet"]
dileptonRazorBoxes = ["MuEle", "MuMu", "EleEle"]

dileptonSignalRegionCuts = "MR > 150.000000 && MR < 4000 && Rsq > 0.40000 && Rsq < 2.5"
leptonicSignalRegionCuts = "MR > 150.000000 && MR < 4000 && Rsq > 0.40000 && Rsq < 2.5 && mT > 120"
looseLeptonSignalRegionCuts = "MR > 150.000000 && MR < 4000 && Rsq > 0.40000 && Rsq < 2.5 && mTLoose > 100 && nJets80 >= 2"
hadronicSignalRegionCuts = "MR > 150.000000 && MR < 4000 && Rsq > 0.40000 && Rsq < 2.5 && nJets80 >= 2"
#dileptonSignalRegionCuts = "MR > 150.000000 && MR < 4000 && Rsq > 0.40000 && Rsq < 1.5 && abs(dPhiRazor) < 2.8"
#leptonicSignalRegionCuts = "MR > 150.000000 && MR < 4000 && Rsq > 0.40000 && Rsq < 1.5 && mT > 120"
#looseLeptonSignalRegionCuts = "MR > 150.000000 && MR < 4000 && Rsq > 0.40000 && Rsq < 1.5 && mTLoose > 100 && nJets80 >= 2"
#hadronicSignalRegionCuts = "MR > 150.000000 && MR < 4000 && Rsq > 0.40000 && Rsq < 1.5 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"

for box in razorBoxes:
    if box in hadronicRazorBoxes: 
        razorCuts[box] = appendBoxCuts(hadronicSignalRegionCuts, razorBoxes[box])
    elif box in looseLeptonRazorBoxes:
        razorCuts[box] = appendBoxCuts(looseLeptonSignalRegionCuts, razorBoxes[box])
    elif box in leptonicRazorBoxes:
        razorCuts[box] = appendBoxCuts(leptonicSignalRegionCuts, razorBoxes[box])
    elif box in dileptonRazorBoxes:
        razorCuts[box] = appendBoxCuts(dileptonSignalRegionCuts, razorBoxes[box])
razorCuts['FourToSixJet'] += ' && nSelectedJets < 7'
razorCuts['SevenJet'] += ' && nSelectedJets >= 7'

#####################################
### BINNING
#####################################

razorBinning = {}

### Single Lepton Control Region
razorBinning["SingleLepton"] = {
        "MR" : [150, 200, 300, 400, 600, 3500],
        "Rsq" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "NJets40" : [1,2,4,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT" : range(0, 2000, 200),
        "MET" : range(0, 1000, 100),
        "lep1.Pt()" : [20,30,40,100,1000],
        "abs(lep1.Eta())" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }
razorBinning["SingleLeptonInv"] = {
        "MR_NoW" : [150, 200, 300, 400, 600, 3500],
        "Rsq_NoW" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "NJets_NoW" : [1,2,4,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT_NoW" : range(0, 2000, 200),
        "MET_NoW" : range(0, 1000, 100),
        }

### Veto lepton control region
razorBinning["VetoLepton"] =  {
        "MR" : [150, 200, 300, 400, 600, 3500],
        "Rsq" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "NJets40" : [1,2,4,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "lep1.Pt()" : [5, 10, 15, 20.,30.,40.,100,1000],
        "abs(lep1.Eta())" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }
razorBinning["SignalRegionForVetoLepton"] =  {
        "MR" : [150, 200, 300, 400, 600, 3500],
        "Rsq" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "nSelectedJets" : [1,2,4,20],
        "leadingGenLeptonPt" : [5, 10, 15, 20.,30.,40.,100,1000],
        "abs(leadingGenLeptonEta)" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }

### Veto Tau Control Region
razorBinning["VetoTau"] = {
        "MR" : [150, 200, 300, 400, 600, 3500],
        "Rsq" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "NJets40" : [1,2,4,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "lep1.Pt()" : [20,30,40,100,1000],
        "abs(lep1.Eta())" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }
razorBinning["SignalRegionForVetoTau"] = {
        "MR" : [150, 200, 300, 400, 600, 3500],
        "Rsq" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "nSelectedJets" : [1,2,4,20],
        "leadingGenLeptonPt" : [20,30,40,100,1000],
        "abs(leadingGenLeptonEta)" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }

### DYJets Dilepton Invisible Control Region
razorBinning["TTJetsDilepton"] = {
        "MR" : [150, 200, 300, 400, 600, 3500],
        "Rsq" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "NJets40" : [1,2,3,4,5,6,20],
        "lep1.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep1.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "MET" : range(250, 1000, 50),
        "NBJetsMedium" : [0,1,2,3,4],
        }
razorBinning["TTJetsDileptonReduced"] = {
        "MR" : [150, 200, 300, 400, 600, 3500],
        "Rsq" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "NJets40" : [4,5,6,20],
        "lep1.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep1.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "MET" : range(250, 1000, 50),
        "NBJetsMedium" : [0,1,2,3,4],
        }

### DYJets Dilepton Invisible Control Region
razorBinning["DYJetsDileptonInv"] = {
        "MR_NoZ" : [150, 200, 300, 400, 600, 3500],
        "Rsq_NoZ" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "NJets_NoZ" : [1,2,4,20],
        "1" : [0.5, 1.5],
        "lep1.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep1.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "lep2.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep2.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT_NoZ" : range(0, 2000, 200),
        "MET_NoZ" : range(0, 1000, 100),
        }
razorBinning["DYJetsDileptonInvReduced"] = {
        "MR_NoZ" : [150, 200, 300, 400, 600, 3500],
        "Rsq_NoZ" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "NJets_NoZ" : [1,2,4,20],
        "1" : [0.5, 1.5],
        "lep1.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep1.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "lep2.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep2.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT_NoZ" : range(0, 2000, 200),
        "MET_NoZ" : range(0, 1000, 100),
        }

### Photon control region
razorBinning["Photon"] = {
        "MR_NoPho" : [150, 200, 300, 400, 600, 3500],
        "Rsq_NoPho" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "NJets_NoPho" : [0,4,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT_NoPho" : range(0, 2000, 200),
        "MET_NoPho" : range(0, 1000, 100),
        "1" : [0.5, 1.5],
        "pho1.Pt()" : range(185, 1500, 50),
        "pho1.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "jet1.Pt()" : range(40, 1500, 50),
        "jet1.Eta()" : [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        "jet2.Pt()" : range(40, 1500, 50),
        "jet2.Eta()" : [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        }

### For October 2016 SUSY Synchronization exercise 
razorBinning["SusySync"] = {
        "MR" : [150, 200, 300, 400, 600, 3500],
        "Rsq" : [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        "NJets40" : [2,3,4,5,6,7,8,9,10],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT" : range(250, 2050, 200),
        "MET" : range(250, 1000, 50),
        "lep1.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep1.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "NPV" : range(60),
        }

### Signal region
signalConfig = "config/run2_2017_01_07_Run2016G_SUSYUnblind_Sep23ReReco.config"
cfg = Config.Config(signalConfig)
for box in ["MultiJet", "DiJet", "LeptonMultiJet", "LeptonJet"]:
    razorBinning[box] = {
            "MR" : cfg.getBinning(box)[0],
            "Rsq" : cfg.getBinning(box)[1],
            }

### Add 2D razor binning to each region
for region,binning in razorBinning.iteritems():
    if "MR" in binning and "Rsq" in binning:
        binning[("MR","Rsq")] = []
    if "MR_NoW" in binning and "Rsq_NoW" in binning:
        binning[("MR_NoW","Rsq_NoW")] = []
    if "MR_NoZ" in binning and "Rsq_NoZ" in binning:
        binning[("MR_NoZ","Rsq_NoZ")] = []
    if "MR_NoPho" in binning and "Rsq_NoPho" in binning:
        binning[("MR_NoPho","Rsq_NoPho")] = []
    if "HT" in binning and "MET" in binning:
        binning[("HT","MET")] = []
    if "HT" in binning and "met" in binning:
        binning[("HT","met")] = []

### Non-grid binning 
xbinsSignal = { "MultiJet":{}, "MuMultiJet":{}, "EleMultiJet":{}, "LeptonMultiJet":{},
                "DiJet":{}, "MuJet":{}, "EleJet":{}, "FourToSixJet":{}, "SevenJet":{}}
colsSignal = { "MultiJet":{}, "MuMultiJet":{}, "EleMultiJet":{}, "LeptonMultiJet":{},
                "DiJet":{}, "MuJet":{}, "EleJet":{}, "FourToSixJet":{}, "SevenJet":{}}

xbinsSignal["MultiJet"]["0B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["MultiJet"]["1B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["MultiJet"]["2B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["MultiJet"]["3B"] = [150, 200, 300, 400, 600, 3500]

xbinsSignal["FourToSixJet"] = xbinsSignal["MultiJet"]

xbinsSignal["MuMultiJet"]["0B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["MuMultiJet"]["1B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["MuMultiJet"]["2B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["MuMultiJet"]["3B"] = [150, 200, 300, 400, 600, 3500]

xbinsSignal["EleMultiJet"]["0B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["EleMultiJet"]["1B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["EleMultiJet"]["2B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["EleMultiJet"]["3B"] = [150, 200, 300, 400, 600, 3500]

xbinsSignal["DiJet"] = xbinsSignal["MultiJet"]
xbinsSignal["MuJet"] = xbinsSignal["MuMultiJet"]
xbinsSignal["EleJet"] = xbinsSignal["EleMultiJet"]
xbinsSignal["LeptonMultiJet"] = xbinsSignal["MuMultiJet"]
xbinsSignal["LeptonJet"] = xbinsSignal["LeptonMultiJet"]

colsSignal["MultiJet"]["0B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["MultiJet"]["1B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["MultiJet"]["2B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["MultiJet"]["3B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["MuMultiJet"]["0B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["MuMultiJet"]["1B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["MuMultiJet"]["2B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["MuMultiJet"]["3B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["EleMultiJet"]["0B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["EleMultiJet"]["1B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["EleMultiJet"]["2B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["EleMultiJet"]["3B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]

colsSignal["FourToSixJet"] = colsSignal["MultiJet"]

colsSignal["DiJet"] = colsSignal["MultiJet"]
colsSignal["MuJet"] = colsSignal["MuMultiJet"]
colsSignal["EleJet"] = colsSignal["EleMultiJet"]
colsSignal["LeptonMultiJet"] = colsSignal["MuMultiJet"]
colsSignal["LeptonJet"] = colsSignal["LeptonMultiJet"]

xbinsSignal["SevenJet"]["0B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["SevenJet"]["1B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["SevenJet"]["2B"] = [150, 200, 300, 400, 600, 3500]
xbinsSignal["SevenJet"]["3B"] = [150, 200, 300, 400, 600, 3500]

colsSignal["SevenJet"]["0B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["SevenJet"]["1B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["SevenJet"]["2B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]
colsSignal["SevenJet"]["3B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
        ]

#binning for scale factor histograms
xbinsSignal["TTJetsSingleLepton"] = {}
colsSignal["TTJetsSingleLepton"] = {}
xbinsSignal["TTJetsSingleLepton"]["0B"] = [150, 200, 300, 400, 600, 3500]
colsSignal["TTJetsSingleLepton"]["0B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
    ]

xbinsSignal["WJetsSingleLepton"] = {}
colsSignal["WJetsSingleLepton"] = {}
xbinsSignal["WJetsSingleLepton"]["0B"] = [150, 200, 300, 400, 600, 3500]
colsSignal["WJetsSingleLepton"]["0B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
    ]

xbinsSignal["WJetsSingleLeptonInv"] = {}
colsSignal["WJetsSingleLeptonInv"] = {}
xbinsSignal["WJetsSingleLeptonInv"]["0B"] = [150, 200, 300, 400, 600, 3500]
colsSignal["WJetsSingleLeptonInv"]["0B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
    ]

xbinsSignal["GJetsInv"] = {}
colsSignal["GJetsInv"] = {}
xbinsSignal["GJetsInv"]["0B"] = [150, 200, 300, 400, 600, 3500]
colsSignal["GJetsInv"]["0B"] = [
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.575,0.65,0.75,0.85,0.95,2.5],
        [0.40,0.5,0.6,0.60,0.7,0.95,2.5],
    ]

class Analysis:
    """Class to hold analysis cuts, binning, input files, and trigger info"""
    def __init__(self, region, tag, njetsMin=-1, njetsMax=-1, nbMin=-1, nbMax=-1):
        self.region = region
        self.tag = tag
        self.njetsMin = njetsMin
        self.njetsMax = njetsMax
        self.nbMin = nbMin
        self.nbMax = nbMax

        if tag == "Razor2015":
            self.lumi = 2300
        elif tag == "Razor2016" or tag == "Razor2016_80X":
            self.lumi = 26400
        elif tag == "Razor2016_ICHEP_80X":
            self.lumi = 12900
        elif tag == "Razor2016G_SUSYUnblind_80X":
            self.lumi = 4394
        elif tag == "Razor2016_MoriondRereco":
            self.lumi = 36800
        else:
            sys.exit("Error: tag"+tag+"is not supported!")

        #load information for analysis
        self.getAnalysis()
        #weight options
        self.getWeightOpts()

    def getAnalysis(self):
        btag = str( max(0, self.nbMin) )+"B"
        self.extraCuts = {}
        self.extraWeightOpts = razorExtraWeightOpts[self.tag]
        if self.region == "SingleLepton":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLepton"][self.tag]
            self.samples = razorSamples["TTJetsSingleLepton"]
            self.samplesReduced = razorSamplesReduced["TTJetsSingleLepton"]
            self.cuts = razorCuts["SingleLepton"]
            self.binning = razorBinning["SingleLepton"]
            self.unrollBins = (xbinsSignal["TTJetsSingleLepton"]["0B"], colsSignal["TTJetsSingleLepton"]["0B"])

        elif self.region == "SingleLeptonInv":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoW"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLeptonInv"][self.tag]
            self.samples = razorSamples["WJetsSingleLeptonInv"]
            self.samplesReduced = razorSamplesReduced["WJetsSingleLeptonInv"]
            self.cuts = razorCuts["SingleLeptonInv"]
            self.binning = razorBinning["SingleLeptonInv"]
            self.unrollBins = (xbinsSignal["WJetsSingleLeptonInv"]["0B"], colsSignal["WJetsSingleLeptonInv"]["0B"])

        elif self.region == "TTJetsSingleLepton":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLepton"][self.tag]
            self.samples = razorSamples["TTJetsSingleLepton"]
            self.samplesReduced = razorSamplesReduced["TTJetsSingleLepton"]
            self.cuts = razorCuts["SingleLepton"]
            self.binning = razorBinning["SingleLepton"]
            self.unrollBins = (xbinsSignal["TTJetsSingleLepton"]["0B"], colsSignal["TTJetsSingleLepton"]["0B"])

        elif self.region == "WJetsSingleLepton":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLepton"][self.tag]
            self.samples = razorSamples["WJetsSingleLepton"]
            self.samplesReduced = razorSamplesReduced["WJetsSingleLepton"]
            self.cuts = razorCuts["SingleLepton"]
            self.binning = razorBinning["SingleLepton"]
            self.unrollBins = (xbinsSignal["WJetsSingleLepton"]["0B"], colsSignal["WJetsSingleLepton"]["0B"])

        elif self.region == "WJetsSingleLeptonInv":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoW"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLeptonInv"][self.tag]
            self.samples = razorSamples["WJetsSingleLeptonInv"]
            self.samplesReduced = razorSamplesReduced["WJetsSingleLeptonInv"]
            self.cuts = razorCuts["SingleLeptonInv"]
            self.binning = razorBinning["SingleLeptonInv"]
            self.unrollBins = (xbinsSignal["WJetsSingleLeptonInv"]["0B"], colsSignal["WJetsSingleLeptonInv"]["0B"])

        elif self.region == "TTJetsDilepton":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["Dilepton"][self.tag]
            self.samples = razorSamples["TTJetsDilepton"]
            self.samplesReduced = razorSamplesReduced["TTJetsDilepton"]
            self.cuts = razorCuts["TTJetsDilepton"]
            self.binning = razorBinning["TTJetsDilepton"]
            self.unrollBins = (None,None)

        elif self.region == "TTJetsDileptonMultiJet":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["Dilepton"][self.tag]
            self.samples = razorSamples["TTJetsDilepton"]
            self.samplesReduced = razorSamplesReduced["TTJetsDilepton"]
            self.cuts = razorCuts["TTJetsDilepton"]
            self.binning = razorBinning["TTJetsDileptonReduced"]
            self.unrollBins = (None,None)

        elif self.region == "DYJetsDileptonInv":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoZ"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["DileptonInv"][self.tag]
            self.samples = razorSamples["DYJetsDileptonInv"]
            self.samplesReduced = razorSamplesReduced["DYJetsDileptonInv"]
            self.cuts = razorCuts["DYJetsDileptonInv"]
            self.binning = razorBinning["DYJetsDileptonInv"]
            self.unrollBins = (None,None)

        elif self.region == "DYJetsDileptonInvMultiJet":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoZ"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["DileptonInv"][self.tag]
            self.samples = razorSamples["DYJetsDileptonInv"]
            self.samplesReduced = razorSamplesReduced["DYJetsDileptonInv"]
            self.cuts = razorCuts["DYJetsDileptonInv"]
            self.binning = razorBinning["DYJetsDileptonInvReduced"]
            self.unrollBins = (None,None)

        elif self.region == "GJetsInv":
            self.trigType = "Photon"
            self.jetVar = "NJets_NoPho"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["Photon"][self.tag]
            self.samples = razorSamples["Photon"]
            self.samplesReduced = razorSamplesReduced["Photon"]
            self.cuts = razorCuts["Photon"]
            self.binning = razorBinning["Photon"]
            self.unrollBins = (xbinsSignal["GJetsInv"]["0B"],colsSignal["GJetsInv"]["0B"])
            self.extraCuts = razorExtraCuts["Photon"]

        elif self.region == "SusySync":
            self.trigType = "SusySync"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SusySync"][self.tag]
            self.samples = razorSamples["SusySync"]
            self.samplesReduced = razorSamplesReduced["SusySync"]
            self.cuts = razorCuts["SusySync"]
            self.binning = razorBinning["SusySync"]
            self.unrollBins = (None,None)

        elif self.region == "VetoLeptonControlRegion":
            self.trigType = "Razor"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["VetoLepton"][self.tag]
            self.samples = razorSamples["VetoLepton"]
            self.samplesReduced = razorSamplesReduced["VetoLepton"]
            self.cuts = razorCuts["VetoLepton"]
            self.binning = razorBinning["VetoLepton"]
            self.unrollBins = (None,None)

        elif self.region == "MultiJetForVetoLeptonControlRegion":
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["MultiJetForVetoLepton"]
            self.binning = razorBinning["SignalRegionForVetoLepton"]
            self.unrollBins = (None,None)
        
        elif self.region == "DiJetForVetoLeptonControlRegion":
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["DiJetForVetoLepton"]
            self.binning = razorBinning["SignalRegionForVetoLepton"]
            self.unrollBins = (None,None)

        elif self.region == "VetoTauControlRegion":
            self.trigType = "Razor"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["VetoTau"][self.tag]
            self.samples = razorSamples["VetoTau"]
            self.samplesReduced = razorSamplesReduced["VetoTau"]
            self.cuts = razorCuts["VetoTau"]
            self.binning = razorBinning["VetoTau"]
            self.unrollBins = (None,None)

        elif self.region == "MultiJetForVetoTauControlRegion":
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["MultiJetForVetoTau"]
            self.binning = razorBinning["SignalRegionForVetoTau"]
            self.unrollBins = (None,None)

        elif self.region == "DiJetForVetoTauControlRegion":
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["DiJetForVetoTau"]
            self.binning = razorBinning["SignalRegionForVetoTau"]
            self.unrollBins = (None,None)

        elif self.region in electronRazorBoxes:
            self.trigType = "SingleLepton"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalElectron"][self.tag]
            self.samples = razorSamples["SignalLeptonic"]
            self.samplesReduced = razorSamplesReduced["SignalLeptonic"]
            self.cuts = razorCuts[self.region]
            self.binning = razorBinning["SignalLeptonic"]
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        elif self.region in muonRazorBoxes:
            self.trigType = "SingleLepton"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalMuon"][self.tag]
            self.samples = razorSamples["SignalLeptonic"]
            self.samplesReduced = razorSamplesReduced["SignalLeptonic"]
            self.cuts = razorCuts[self.region]
            self.binning = razorBinning["SignalLeptonic"]
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        elif self.region in leptonRazorBoxes:
            self.trigType = "SingleLepton"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalLepton"][self.tag]
            self.samples = razorSamples["SignalLeptonic"]
            self.samplesReduced = razorSamplesReduced["SignalLeptonic"]
            self.cuts = razorCuts[self.region]
            self.binning = razorBinning[self.region]
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        elif self.region in hadronicRazorBoxes:
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts[self.region]
            self.binning = razorBinning[self.region]
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        else: 
            sys.exit("Error: analysis region "+self.region+" is not implemented!")

        #add jet and bjet and trigger requirements to cuts
        if self.region in razorBoxes: self.cuts = appendNoiseFilters(self.cuts)
        self.addJetCuts()
        self.addTriggerCuts()

    def addJetCuts(self):
        if self.njetsMin >= 0:
            self.cuts = self.cuts + ' && %s >= %d'%(self.jetVar,self.njetsMin)
        if self.njetsMax >= 0:
            self.cuts = self.cuts + ' && %s <= %d'%(self.jetVar,self.njetsMax)
        if self.nbMin >= 0:
            self.cuts = self.cuts + ' && %s >= %d'%(self.bjetVar,self.nbMin)
        if self.nbMax >= 0:
            self.cuts = self.cuts + ' && %s <= %d'%(self.bjetVar,self.nbMax)

    def addTriggerCuts(self):
        self.trigInfoMC = TriggerUtils(self.trigType, self.tag, isData=False)
        self.trigInfoData = TriggerUtils(self.trigType, self.tag, isData=True)
        self.cutsMC = self.trigInfoMC.appendTriggerCuts( self.cuts )+' && TMath::Finite(weight) && (!TMath::IsNaN(weight))'
        self.cutsData = self.trigInfoData.appendTriggerCuts( self.cuts )
        if self.region == 'GJetsInv':
            self.cutsData += ' && pho1HLTFilter[36]'

    def getWeightOpts(self):
        self.weightOpts = copy.copy(razorWeightOpts[self.tag])
        self.dataWeightOpts = []
        self.weightFiles = {}
        self.weightHists = {}
        for wType,wInfo in razorWeightHists[self.tag].iteritems():
            self.weightFiles[wType] = rt.TFile.Open(wInfo[0])
            self.weightHists[wType] = self.weightFiles[wType].Get(wInfo[1])

