## Inclusive razor analysis definitions
import argparse
import ROOT as rt
from array import array
import sys
import copy
import csv

#local imports
from framework import Config

def make_parser():
    """
    Gets a basic command line parser with common arguments needed for
    razor analysis scripts
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument("--tag", help="Analysis tag, e.g. Razor2015", default="Razor2016_MoriondRereco")
    parser.add_argument('--no-fill', help="dry run -- do not fill histograms", action="store_true", 
            dest='noFill')
    parser.add_argument("--no-save", dest="noSave", action="store_true", help="Do not save SFs or histograms")
    parser.add_argument('--no-boost-cuts', dest='noBoostCuts', action='store_true')
    return parser

#####################################
### WEIGHTS NEEDED
#####################################

razorWeightOpts = {
        "Razor2015":[],
        "Razor2016":[], 
        }
razorWeightOpts["Razor2016G_SUSYUnblind_80X"] = razorWeightOpts["Razor2016"]
razorWeightOpts["Razor2016_MoriondRereco"] = razorWeightOpts["Razor2016"]+['boost']
razorWeightOpts["Razor2016_80X"] = razorWeightOpts["Razor2016"]
razorWeightOpts["Razor2016_ICHEP_80X"] = razorWeightOpts["Razor2016"]
razorExtraWeightOpts = {
        "Razor2016G_SUSYUnblind_80X":{
            #'TTJets':['toppt'], 
            #'TTJets1L':['toppt'], 
            #'TTJets2L':['toppt'],
            'Signal':['nisr'],
                 }
        }
razorExtraWeightOpts["Razor2016_MoriondRereco"] = razorExtraWeightOpts["Razor2016G_SUSYUnblind_80X"]
razorExtraWeightOpts["Razor2015"] = {}
razorExtraWeightOpts["Razor2016_80X"] = razorExtraWeightOpts["Razor2016G_SUSYUnblind_80X"]
razorExtraWeightOpts["Razor2016_ICHEP_80X"] = razorExtraWeightOpts["Razor2016G_SUSYUnblind_80X"]
razorWeightHists = {
        "Razor2015":{},
        "Razor2016":{ 
            "qcdfunction{}".format(box):(
            "macros/BackgroundStudies/QCD/qcd_best_fit_2d_{}.root".format(box),
            "qcd{}".format(box)) for box in ['dijet', 'multijet', 'sevenjet'] }
        }
razorWeightHists["Razor2016"].update({
    "qcdaltfunction{}".format(box):(
        "macros/BackgroundStudies/QCD/qcd_best_fit_2d_{}_quadratic.root".format(box), "qcd{}".format(box)) for box in ['dijet', 'multijet', 'sevenjet'] })
razorWeightHists["Razor2016"].update({
    "qcdbtags{}".format(box):(
        "data/ScaleFactors/RazorMADD2015/RazorQCDBTagScaleFactors_Razor2016_MoriondRereco.root",
        "{}btags".format(box)) for box in ['dijet', 'multijet', 'sevenjet']})
razorWeightHists["Razor2016"].update({
    "qcdbtags{}sys".format(box):(
        "data/ScaleFactors/RazorMADD2015/RazorQCDBTagScaleFactors_Razor2016_MoriondRereco.root",
        "{}btagssys".format(box)) for box in ['dijet', 'multijet', 'sevenjet']})
razorWeightHists["Razor2016G_SUSYUnblind_80X"] = {}
razorWeightHists["Razor2016_MoriondRereco"] = copy.copy(razorWeightHists['Razor2016'])
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
dirCR2016 = "/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/"
dirSR2016 = "/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/"
versionMC2016 = "V3p15_27Nov2017"
versionMCSR2016 = "V3p15_27Nov2017"
versionData2016 = "V3p15_27Nov2017"

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
dir1L2016 = dirCR2016+'/'+versionMC2016+'/OneLeptonFull'
dir1LInv2016 = dirCR2016+'/'+versionMC2016+'/OneLeptonAddToMET'
dir2L2016 = dirCR2016+'/'+versionMC2016+'/DileptonFull'
dir2LInv2016 = dirCR2016+'/'+versionMC2016+'/DileptonAddToMET'
dirVetoL2016 = dirCR2016+'/'+versionMC2016+'/VetoLepton'
dirVetoTau2016 = dirCR2016+'/'+versionMC2016+'/VetoTau'
dirPhoton2016 = dirCR2016+'/'+versionMC2016+'/PhotonAddToMET'
dirSignal2016 = dirSR2016+'/'+versionMCSR2016+'/Signal'
dirSusySync2016 = "eos/cms/store/group/phys_susy/razor/Run2Analysis/SusySync/2016/V3p6_25October2016_CustomType1MET/OneLeptonFull/"

#local directories
#dir1L2016 = 'Backgrounds/1L'
#dir2L2016 = 'Backgrounds/2L'
#dir1LInv2016 = 'Backgrounds/1LInv'
#dir2LInv2016 = 'Backgrounds/2LInv'
#dirVetoL2016 = 'Backgrounds/VetoL'
#dirVetoTau2016 = 'Backgrounds/VetoTau'
#dirSignal2016 = 'Backgrounds/Signal'
#dirPhoton2016 = 'Backgrounds/Photon'    

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
            #"QCD"      : dirSignal2016+"/"+prefixes2016[tag]["Signal"]+"_QCD_1pb_weighted"+skimstr+".root",
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
            if 'QCD' in files[tag] and 'Data' in files[tag]['QCD']:
                files[tag]['QCD'] = files[tag]['QCD'].replace(
                        '/'+versionMC2016+'/', '/'+versionData2016+'/')

    #Signal ntuples
    razorSignalDirs = {
            "Razor2016_MoriondRereco": "/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_27Nov2017/SignalFastsim/"
            }


#####################################
### FITS
#####################################

razorFitDirs = { 
        "Razor2016G_SUSYUnblind_80X":"/afs/cern.ch/work/j/jlawhorn/public/Razor_Moriond2017/CMSSW_7_1_5/src/RazorAnalyzer/fits_2017_01_09/ReReco2016_02Jan/",
        "Razor2016_MoriondRereco":"/afs/cern.ch/work/j/jlawhorn/public/Razor_Moriond2017/newFitFxn/CMSSW_7_1_5/src/RazorAnalyzer/Plots/Razor2016_MoriondRereco"
        }
razorFitFiles = { tag:{} for tag in razorFitDirs }
for tag,path in razorFitDirs.iteritems():
    for box in ["MultiJet","MultiJet_0b","MultiJet_1b","MultiJet_2b","DiJet","DiJet_0b","DiJet_1b","DiJet_2b","LeptonMultiJet","LeptonMultiJet_0b","LeptonMultiJet_1b","LeptonMultiJet_2b","LeptonJet","LeptonJet_0b","LeptonJet_1b","LeptonJet_2b"]:
        razorFitFiles[tag][box] = "%s/toys_Bayes_noStat_%s_combinedBtags.root"%(path,box)
        #razorFitFiles[tag][box] = "%s/Fits_%s/toys_Bayes_noStat_%s_combinedBtags.root"%(path,box,box)


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
triggerNames["Razor"]["Data"]["Razor2016"] = [
        'HLT_PFHT900',
        'HLT_RsqMR270_Rsq0p09_MR200',
        'HLT_RsqMR270_Rsq0p09_MR200_4jet',
        ]
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

recommendedNoiseFilters = [
        "Flag_HBHENoiseFilter","Flag_HBHEIsoNoiseFilter",
        "Flag_goodVertices", "Flag_eeBadScFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_CSCTightHaloFilter",
        "Flag_badChargedCandidateFilter","Flag_badMuonFilter",
        "!Flag_badGlobalMuonFilter", "!Flag_duplicateMuonFilter"
        ]
def appendNoiseFilters(cuts, tree=None):
    ret = copy.copy(cuts)
    for bit in recommendedNoiseFilters:
        if (tree is None) or hasattr(tree, bit):
            ret += " && " + bit
    return ret

# 1-lepton control region
razorCuts["SingleLepton"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && ((abs(lep1Type) == 11 && lep1.Pt() > 30) || (abs(lep1Type) == 13 && lep1.Pt() > 25)) && MET > 30 && lep1MT > 30 && lep1MT < 100 && MR > 300 && Rsq > 0.15"

# 1-lepton invisible control region
razorCuts["SingleLeptonInv"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && ((abs(lep1Type) == 11 && lep1.Pt() > 30) || (abs(lep1Type) == 13 && lep1.Pt() > 25)) && MET > 30 && lep1MT > 30 && lep1MT < 100 && NJets80_NoW >= 2 && MR_NoW > 300 && Rsq_NoW > 0.15 && ( weight < 0.01 || weight == 1)"

# TTJets 2-lepton control region
razorCuts["TTJetsDilepton"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && (abs(lep2Type) == 11 || abs(lep2Type) == 13) && lep1PassTight && lep2PassTight && lep1.Pt() > 30 && lep2.Pt() > 30 && mll > 20 && ((abs(lep1Type) != abs(lep2Type)) || (mll < 76 || mll > 106)) && NBJetsMedium > 0 && MET > 40 && MR > 300 && Rsq > 0.15 && Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter"

# DYJets 2-lepton invisible control region
razorCuts["DYJetsDileptonInv"] = "((abs(lep1Type) == 11 && abs(lep2Type) == 11) || (abs(lep1Type) == 13 && abs(lep2Type) == 13)) && lep1.Pt() > 30 && lep2.Pt() > 20 && recoZmass > 80 && recoZmass < 110 && NBJetsMedium == 0 && NJets80_NoZ >= 2 && MR_NoZ > 400 && Rsq_NoZ > 0.25"

# Veto lepton control region
razorCuts["VetoLepton"]  = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassVeto && lep1.Pt() > 5 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 400 && Rsq > 0.25"

razorCuts["VetoElectron"] = "(abs(lep1Type) == 11) && lep1PassVeto && lep1.Pt() > 5 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 400 && Rsq > 0.25"

razorCuts["VetoMuon"] = "(abs(lep1Type) == 13) && lep1PassVeto && lep1.Pt() > 5 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 400 && Rsq > 0.25"

razorCuts["VetoTau"] = "(abs(lep1Type) == 15) && lep1PassLoose && lep1.Pt() > 20 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 400 && Rsq > 0.25"

# Photon control region
razorCuts["Photon"] = "pho1.Pt() > 185 && pho1_chargediso < 2.5 && ((abs(pho1.Eta()) < 1.5 && pho1_sigmaietaieta < 0.01031) || (abs(pho1.Eta()) >= 1.5 && pho1_sigmaietaieta < 0.03013)) && NJets80_NoPho >= 2 && MR_NoPho > 400 && Rsq_NoPho > 0.25 && Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter"
razorExtraCuts["Photon"] = {
        "GJetsInv":"!(abs(pho1_motherID) > 5 && pho1_motherID != 21 && pho1_motherID != 2212) && minDRGenPhotonToParton >= 0.4",
        "GJetsFrag":"!(abs(pho1_motherID) > 5 && pho1_motherID != 21 && pho1_motherID != 2212) && minDRGenPhotonToParton < 0.4",
        }

# For October 2016 SUSY synchronization exercise
razorCuts["SusySync"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && lep1.Pt() > 25 && abs(lep1.Eta()) < 2.4 && MET > 250 && HT > 250 && NJets40 >= 2 && Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter"

# Hadronic boxes with inverted lepton veto
cutsMultiJetForVeto = "(box == 11 || box == 12 || box == 21) && MR > 400.000000 && Rsq > 0.250000 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"
cutsDiJetForVeto = "(box == 14) && MR > 400.000000 && Rsq > 0.250000 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"

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
        "LeptonJet_0b" : [5,8],
        "LeptonJet_1b" : [5,8],
        "LeptonJet_2b" : [5,8],
        "LooseLeptonSixJet" : [9],
        "LooseLeptonFourJet" : [10],
        "LooseLeptonDiJet" : [13],
        "SixJet" : [11],
        "FourJet" : [12],
        "DiJet" : [14],	  
        "DiJet_0b" : [14],	  
        "DiJet_1b" : [14],	  
        "DiJet_2b" : [14],	  
        "MuMultiJet" : [3,4,18],
        "EleMultiJet" : [6,7,19],
        "LeptonMultiJet" : [3,4,18,6,7,19],
        "LeptonSevenJet" : [3,4,18,6,7,19],
        "LeptonMultiJet_0b" : [3,4,18,6,7,19],
        "LeptonMultiJet_1b" : [3,4,18,6,7,19],
        "LeptonMultiJet_2b" : [3,4,18,6,7,19],
        "LooseLeptonMultiJet" : [9,10,20],
        "MultiJet" : [11,12,21],
        "MultiJet_0b" : [11,12,21],
        "MultiJet_1b" : [11,12,21],
        "MultiJet_2b" : [11,12,21],
        "FourToSixJet" : [11,12,21],
        "SevenJet" : [11,12,21],
        }

hadronicRazorBoxes = ["DiJet", "FourJet", "SixJet", "MultiJet", "FourToSixJet", "SevenJet","DiJet_0b","DiJet_1b","DiJet_2b","MultiJet_0b","MultiJet_1b","MultiJet_2b"]
looseLeptonRazorBoxes = ["LooseLeptonDiJet", "LooseLeptonFourJet", "LooseLeptonSixJet", "LooseLeptonMultiJet"]
leptonicRazorBoxes = ["MuJet", "MuFourJet", "MuSixJet", "MuMultiJet","EleJet", "EleFourJet", "EleSixJet", "EleMultiJet", "LeptonJet", "LeptonFourJet", "LeptonSixJet", "LeptonMultiJet","LeptonSevenJet","LeptonJet_0b","LeptonJet_1b","LeptonJet_2b","LeptonMultiJet_0b","LeptonMultiJet_1b","LeptonMultiJet_2b"]
electronRazorBoxes = ["EleJet", "EleFourJet", "EleSixJet", "EleMultiJet"]
muonRazorBoxes = ["MuJet", "MuFourJet", "MuSixJet", "MuMultiJet"]
leptonRazorBoxes = ["LeptonJet", "LeptonFourJet", "LeptonSixJet", "LeptonMultiJet","LeptonSevenJet","LeptonJet_0b","LeptonJet_1b","LeptonJet_2b","LeptonMultiJet_0b","LeptonMultiJet_1b","LeptonMultiJet_2b"]
dileptonRazorBoxes = ["MuEle", "MuMu", "EleEle"]

dileptonSignalRegionCuts = "MR > 400 && MR < 4000 && Rsq > 0.15 && Rsq < 1.5 && abs(dPhiRazor) < 2.8"
leptonicSignalRegionCuts = "MR > 550 && MR < 4000 && Rsq > 0.20 && Rsq < 1.5 && mT > 120 && nSelectedJets >= 2"
looseLeptonSignalRegionCuts = "MR > 500 && MR < 4000 && Rsq > 0.25 && Rsq < 1.5 && mTLoose > 100 && nJets80 >= 2"
hadronicSignalRegionCuts = "MR > 650 && MR < 4000 && Rsq > 0.30 && Rsq < 1.5 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"

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
razorCuts['MultiJet'] += ' && nSelectedJets < 7'
razorCuts['LeptonMultiJet'] += ' && nSelectedJets < 7'
razorCuts['SevenJet'] += ' && nSelectedJets >= 7'
razorCuts['LeptonSevenJet'] += ' && nSelectedJets >= 7'
for b in ['MultiJet', 'SevenJet']:
    razorCuts[b] = razorCuts[b].replace('Rsq > 0.30',
            '(Rsq > 0.30 || (Rsq > 0.2 && MR > 1600))')

zeroBjetBoxes = ["MultiJet_0b", "LeptonJet_0b","LeptonMultiJet_0b","DiJet_0b"]
oneBjetBoxes  = ["MultiJet_1b", "LeptonJet_1b","LeptonMultiJet_1b","DiJet_1b"]
twoBjetBoxes  = ["MultiJet_2b", "LeptonJet_2b","LeptonMultiJet_2b","DiJet_2b"]

zeroBjetCut = "nBTaggedJets==0"
oneBjetCut  = "nBTaggedJets==1"
twoBjetCut  = "nBTaggedJets>=2"

for box in razorBoxes:
    if box in zeroBjetBoxes:
        razorCuts[box] += ' && nBTaggedJets==0'
    elif box in oneBjetBoxes:
        razorCuts[box] += ' && nBTaggedJets==1'
    elif box in twoBjetBoxes:
        razorCuts[box] += ' && nBTaggedJets>=2'

razorCuts["SingleLeptonForSignal"] = "MR > 300 && Rsq > 0.15 && (leadingTightMuPt > 25 || leadingTightElePt > 30) && met > 30 && mT > 30 && mT < 100"
razorCuts["SingleLeptonForSignal"] = appendBoxCuts(
        razorCuts["SingleLeptonForSignal"], 
        razorBoxes['LeptonJet']+razorBoxes['LeptonMultiJet'])

#####################################
### BINNING
#####################################

razorBinning = {}

### Single Lepton Control Region
razorBinning["SingleLepton"] = {
        "MR" : [300, 400, 500, 600, 700, 900, 1200, 4000],
        "Rsq" : [0.15,0.20,0.25,0.30,0.41,0.52,0.64,1.5],
        "NJets40" : [1,2,4,7,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT" : range(0, 2000, 200),
        "MET" : range(0, 500, 50),
        "lep1.Pt()" : [20,30,40,100,1000],
        "abs(lep1.Eta())" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "dPhiRazor":[ 0.2*x for x in range(16) ],
        "lep1MT":range(0,300,10),
        "NJets80" : [0, 1, 2, 3, 4, 5, 6],
        }
razorBinning["SingleLeptonInv"] = {
        "MR_NoW" : [300, 400, 500, 600, 700, 900, 1200, 4000],
        "Rsq_NoW" : [0.15,0.20,0.25,0.30,0.41,0.52,0.64,1.5],
        "NJets_NoW" : [1,2,4,7,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT_NoW" : range(0, 2000, 200),
        "MET_NoW" : range(0, 500, 50),
        }

### Veto lepton control region
razorBinning["VetoLepton"] =  {
        "MR" : [400, 500, 600, 700, 900, 1200, 4000],
        "Rsq" : [0.25,0.30,0.41,0.52,0.64,1.5],
        "NJets40" : [1,2,4,7,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "lep1.Pt()" : [5, 10, 15, 20.,30.,40.,100,1000],
        "abs(lep1.Eta())" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }
razorBinning["SignalRegionForVetoLepton"] =  {
        "MR" : [400, 500, 600, 700, 900, 1200, 4000],
        "Rsq" : [0.25,0.30,0.41,0.52,0.64,1.5],
        "nSelectedJets" : [1,2,4,7,20],
        "leadingGenLeptonPt" : [5, 10, 15, 20.,30.,40.,100,1000],
        "abs(leadingGenLeptonEta)" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }

### Veto Tau Control Region
razorBinning["VetoTau"] = {
        "MR" : [400, 500, 600, 700, 900, 4000],
        "Rsq" : [0.25,0.30,0.41,1.5],
        "NJets40" : [1,2,4,7,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "lep1.Pt()" : [20,30,40,100,1000],
        "abs(lep1.Eta())" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }
razorBinning["SignalRegionForVetoTau"] = {
        "MR" : [400, 500, 600, 700, 900, 4000],
        "Rsq" : [0.25,0.30,0.41,1.5],
        "nSelectedJets" : [1,2,4,7,20],
        "leadingGenLeptonPt" : [20,30,40,100,1000],
        "abs(leadingGenLeptonEta)" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }

### DYJets Dilepton Invisible Control Region
razorBinning["TTJetsDilepton"] = {
        "MR" : [300, 400, 500, 700, 900, 4000],
        "Rsq" : [0.15, 0.20, 0.25, 0.30, 0.41, 1.5],
        "NJets40" : [1,2,3,4,5,6,7,20],
        "lep1.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep1.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "MET" : range(250, 1000, 50),
        "NBJetsMedium" : [0,1,2,3,4],
        }
razorBinning["TTJetsDileptonReduced"] = {
        "MR" : [300, 400, 500, 700, 900, 4000],
        "Rsq" : [0.15, 0.20, 0.30, 1.5],
        "NJets40" : [4,5,6,7,20],
        "lep1.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep1.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "MET" : range(250, 1000, 50),
        "NBJetsMedium" : [0,1,2,3,4],
        }

### DYJets Dilepton Invisible Control Region
razorBinning["DYJetsDileptonInv"] = {
        "MR_NoZ" : [400, 500, 600, 700, 900, 1200, 4000],
        "Rsq_NoZ" : [0.25,0.30,0.41,0.52,0.64,1.5],
        "NJets_NoZ" : [1,2,4,7,20],
        "1" : [0.5, 1.5],
        "lep1.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep1.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "lep2.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep2.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT_NoZ" : range(0, 2000, 200),
        "MET_NoZ" : range(0, 500, 50),
        }
razorBinning["DYJetsDileptonInvReduced"] = {
        "MR_NoZ" : [400, 500, 600, 700, 900, 4000],
        "Rsq_NoZ" : [0.25,0.30,0.41,1.5],
        "NJets_NoZ" : [1,2,4,7,20],
        "1" : [0.5, 1.5],
        "lep1.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep1.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "lep2.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep2.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT_NoZ" : range(0, 2000, 200),
        "MET_NoZ" : range(0, 500, 50),
        }

### Photon control region
razorBinning["Photon"] = {
        "MR_NoPho" : [400, 500, 600, 700, 900, 1200, 4000],
        "Rsq_NoPho" : [0.25, 0.30, 0.41, 0.52, 0.64, 1.5],
        "NJets_NoPho" : [1,2,4,7,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT_NoPho" : range(0, 2000, 200),
        "MET_NoPho" : range(0, 500, 50),
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
        "MR" : [300, 400, 500, 600, 700, 900, 1200, 4000],
        "Rsq" : [0.15,0.20,0.25,0.30,0.41,0.52,0.64,1.5],
        "NJets40" : [2,3,4,5,6,7,8,9,10],
        "NBJetsMedium" : [0,1,2,3,4],
        "HT" : range(250, 2050, 200),
        "MET" : range(250, 1000, 50),
        "lep1.Pt()" : [20.,25.,30.,35.,40.,50.,60.,80.,100,150,200],
        "lep1.Eta()" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "NPV" : range(60),
        }

### Signal region
signalConfig = "config/run2_MADD.config"
cfg = Config.Config(signalConfig)
for box in ["MultiJet","MultiJet_0b","MultiJet_1b","MultiJet_2b","DiJet","DiJet_0b","DiJet_1b","DiJet_2b","LeptonMultiJet","LeptonMultiJet_0b","LeptonMultiJet_1b","LeptonMultiJet_2b","LeptonJet","LeptonJet_0b","LeptonJet_1b","LeptonJet_2b","SevenJet","LeptonSevenJet"]:
    razorBinning[box] = {
            "MR" : cfg.getBinning(box)[0],
            "Rsq" : cfg.getBinning(box)[1],
            }
for box in ['LooseLeptonDiJet', 'LooseLeptonMultiJet']:
    razorBinning[box] = razorBinning[box.replace('LooseLepton','')].copy()

# Config for control regions
controlConfig = "config/razorControlRegions_Razor2016_MoriondRereco.config"
razorBinning["SingleLeptonForSignal"] = {
        "MR" : razorBinning["SingleLepton"]["MR"],
        "Rsq": razorBinning["SingleLepton"]["Rsq"],
        }

# Obligatory super regions
for box in ['DiJet', 'MultiJet', 'SevenJet', 'LeptonMultiJet', 'LeptonSevenJet']:
    razorBinning['Super'+box] = {
            'MR' : [900, 1600, 4000],
            'Rsq': [0.2, 0.3, 1.5]
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
xbinsSignal = { "MultiJet":{}, "MuMultiJet":{}, "EleMultiJet":{}, 
                "LeptonMultiJet":{}, "DiJet":{}, "MuJet":{}, 
                "EleJet":{}, "FourToSixJet":{}, "SevenJet":{},
                "LeptonSevenJet":{} }
colsSignal = { "MultiJet":{}, "MuMultiJet":{}, "EleMultiJet":{}, 
                "LeptonMultiJet":{}, "DiJet":{}, "MuJet":{}, 
                "EleJet":{}, "FourToSixJet":{}, "SevenJet":{},
                "LeptonSevenJet":{} }

xbinsSignalSideband = { "MultiJet":{}, "LeptonMultiJet":{},
        "DiJet":{}, "LeptonJet":{}, "SevenJet":{} }
colsSignalSideband = { "MultiJet":{}, "LeptonMultiJet":{},
        "DiJet":{}, "LeptonJet":{}, "SevenJet":{} }

xbinsSignal["MultiJet"]["0B"] = [ 650, 750, 900, 1200, 1600, 4000 ]
xbinsSignal["MultiJet"]["1B"] = [ 650, 750, 900, 1200, 1600, 4000 ]
xbinsSignal["MultiJet"]["2B"] = [ 650, 750, 900, 1200, 1600, 4000 ]
xbinsSignal["MultiJet"]["3B"] = [ 650, 750, 900, 1200, 1600, 4000 ]

xbinsSignal["FourToSixJet"] = xbinsSignal["MultiJet"]

xbinsSignal["MuMultiJet"]["0B"] = [ 550, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["MuMultiJet"]["1B"] = [ 550, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["MuMultiJet"]["2B"] = [ 550, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["MuMultiJet"]["3B"] = [ 550, 700, 900, 1200, 1600, 4000 ]


xbinsSignal["EleMultiJet"]["0B"] = [ 550, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["EleMultiJet"]["1B"] = [ 550, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["EleMultiJet"]["2B"] = [ 550, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["EleMultiJet"]["3B"] = [ 550, 700, 900, 1200, 1600, 4000 ]

xbinsSignal["DiJet"] = xbinsSignal["MultiJet"]
xbinsSignal["MuJet"] = xbinsSignal["MuMultiJet"]
xbinsSignal["EleJet"] = xbinsSignal["EleMultiJet"]
xbinsSignal["LeptonMultiJet"] = xbinsSignal["MuMultiJet"]

xbinsSignal["LeptonJet"] = xbinsSignal["LeptonMultiJet"]

xbinsSignal["LeptonJet_0b"] = xbinsSignal["LeptonJet"]
xbinsSignal["LeptonJet_1b"] = xbinsSignal["LeptonJet"]
xbinsSignal["LeptonJet_2b"] = xbinsSignal["LeptonJet"]

xbinsSignal["LeptonMultiJet_0b"] = xbinsSignal["LeptonMultiJet"]
xbinsSignal["LeptonMultiJet_1b"] = xbinsSignal["LeptonMultiJet"]
xbinsSignal["LeptonMultiJet_2b"] = xbinsSignal["LeptonMultiJet"]

xbinsSignal["MultiJet_0b"] = xbinsSignal["MultiJet"]
xbinsSignal["MultiJet_1b"] = xbinsSignal["MultiJet"]
xbinsSignal["MultiJet_2b"] = xbinsSignal["MultiJet"]

xbinsSignal["DiJet_0b"] = xbinsSignal["DiJet"]
xbinsSignal["DiJet_1b"] = xbinsSignal["DiJet"]
xbinsSignal["DiJet_2b"] = xbinsSignal["DiJet"]

for box in ['MultiJet']:
    # don't include last bin, it's part of the signal region
    xbinsSignalSideband[box] = { btags:[500, 575]+bins[:-1]
            for btags, bins in xbinsSignal[box].iteritems() }
for box in ['DiJet']:
    xbinsSignalSideband[box] = { btags:[500, 575]+bins 
            for btags, bins in xbinsSignal[box].iteritems() }
for box in ['LeptonMultiJet', 'LeptonJet']:
    xbinsSignalSideband[box] = { btags:[400, 475]+bins 
            for btags, bins in xbinsSignal[box].iteritems() }

colsSignal["MultiJet"]["0B"] = [
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        ]
colsSignal["MultiJet"]["1B"] = [
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        ]
colsSignal["MultiJet"]["2B"] = [
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        ]
colsSignal["MultiJet"]["3B"] = [
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 1.5 ],
        [ 0.30, 0.41, 1.5 ],
        [ 0.20, 0.30, 1.5 ],
        ]
colsSignal["MuMultiJet"]["0B"] = [
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        ]
colsSignal["MuMultiJet"]["1B"] = [
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        ]
colsSignal["MuMultiJet"]["2B"] = [
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        ]
colsSignal["MuMultiJet"]["3B"] = [
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        ]
colsSignal["EleMultiJet"]["0B"] = [
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        ]
colsSignal["EleMultiJet"]["1B"] = [
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        ]
colsSignal["EleMultiJet"]["2B"] = [
        [ 0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 1.5 ],
        [ 0.20, 0.25, 0.30, 0.41, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        ]
colsSignal["EleMultiJet"]["3B"] = [
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        [ 0.20, 0.25, 1.5 ],
        ]

colsSignal["DiJet"]["0B"] = [
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        ]
colsSignal["DiJet"]["1B"] = [
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        ]
colsSignal["DiJet"]["2B"] = [
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.30, 0.41, 0.52, 1.5 ],
        [ 0.30, 0.41, 0.52, 1.5 ],
        ]

colsSignal["FourToSixJet"] = colsSignal["MultiJet"]

colsSignal["MuJet"] = colsSignal["MuMultiJet"]
colsSignal["EleJet"] = colsSignal["EleMultiJet"]
colsSignal["LeptonMultiJet"] = colsSignal["MuMultiJet"]
colsSignal["LeptonJet"] = colsSignal["LeptonMultiJet"]

colsSignal["LeptonJet_0b"] = colsSignal["LeptonJet"]
colsSignal["LeptonJet_1b"] = colsSignal["LeptonJet"]
colsSignal["LeptonJet_2b"] = colsSignal["LeptonJet"]

colsSignal["LeptonMultiJet_0b"] = colsSignal["LeptonMultiJet"]
colsSignal["LeptonMultiJet_1b"] = colsSignal["LeptonMultiJet"]
colsSignal["LeptonMultiJet_2b"] = colsSignal["LeptonMultiJet"]

colsSignal["MultiJet_0b"] = colsSignal["MultiJet"]
colsSignal["MultiJet_1b"] = colsSignal["MultiJet"]
colsSignal["MultiJet_2b"] = colsSignal["MultiJet"]

colsSignal["DiJet_0b"] = colsSignal["DiJet"]
colsSignal["DiJet_1b"] = colsSignal["DiJet"]
colsSignal["DiJet_2b"] = colsSignal["DiJet"]

for box in ['MultiJet', 'DiJet']:
    for btags, bins in colsSignal[box].iteritems():
        colsSignalSideband[box][btags] = [
                [0.25]+bins[0],
                [0.25]+bins[1]
                ]
        colsSignalSideband[box][btags] += [[0.25, 0.30] for _ in range(len(bins))]
        # skip last bin, it's part of the signal region
        if box == 'MultiJet':
            colsSignalSideband[box][btags] = colsSignalSideband[box][btags][:-1]
for box in ['LeptonMultiJet', 'LeptonJet']:
    for btags, bins in colsSignal[box].iteritems():
        colsSignalSideband[box][btags] = [
                [0.15]+bins[0],
                [0.15]+bins[1]
                ]
        colsSignalSideband[box][btags] += [[0.15, 0.20] for _ in range(len(bins))]

xbinsSignal["SevenJet"]["0B"] = [ 650, 900, 1200, 1600, 4000 ]
xbinsSignal["SevenJet"]["1B"] = [ 650, 900, 1200, 1600, 4000 ]
xbinsSignal["SevenJet"]["2B"] = [ 650, 900, 1200, 1600, 4000 ]
xbinsSignal["SevenJet"]["3B"] = [ 650, 1600, 4000 ]
xbinsSignal["LeptonSevenJet"]["0B"] = [ 550, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["LeptonSevenJet"]["1B"] = [ 550, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["LeptonSevenJet"]["2B"] = [ 550, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["LeptonSevenJet"]["3B"] = [ 550, 700, 1600, 4000 ]

colsSignal["SevenJet"]["0B"] = [
        [ 0.30, 0.41, 1.5 ],
        [ 0.30, 0.41, 1.5 ],
        [ 0.30, 0.41, 1.5 ],
        [ 0.20, 1.5 ],
        ]
colsSignal["SevenJet"]["1B"] = [
        [ 0.30, 0.41, 1.5 ],
        [ 0.30, 0.41, 1.5 ],
        [ 0.30, 0.41, 1.5 ],
        [ 0.20, 1.5 ],
        ]
colsSignal["SevenJet"]["2B"] = [
        [ 0.30, 0.41, 1.5 ],
        [ 0.30, 0.41, 1.5 ],
        [ 0.30, 1.5 ],
        [ 0.20, 1.5 ],
        ]
colsSignal["SevenJet"]["3B"] = [
        [ 0.30, 0.41, 1.5 ],
        [ 0.20, 1.5 ],
        ]

colsSignal["LeptonSevenJet"]["0B"] = [
        [ 0.20, 0.30, 1.5 ],
        [ 0.20, 0.30, 1.5 ],
        [ 0.20, 0.30, 1.5 ],
        [ 0.20, 1.5 ],
        [ 0.20, 1.5 ],
        ]
colsSignal["LeptonSevenJet"]["1B"] = [
        [ 0.20, 0.30, 1.5 ],
        [ 0.20, 0.30, 1.5 ],
        [ 0.20, 0.30, 1.5 ],
        [ 0.20, 1.5 ],
        [ 0.20, 1.5 ],
        ]
colsSignal["LeptonSevenJet"]["2B"] = [
        [ 0.20, 0.30, 1.5 ],
        [ 0.20, 0.30, 1.5 ],
        [ 0.20, 0.30, 1.5 ],
        [ 0.20, 1.5 ],
        [ 0.20, 1.5 ],
        ]
colsSignal["LeptonSevenJet"]["3B"] = [
        [ 0.20, 0.30, 1.5 ],
        [ 0.20, 0.30, 1.5 ],
        [ 0.20, 1.5 ],
        ]
xbinsSignalSideband['SevenJet'] = copy.copy(
        xbinsSignalSideband['MultiJet'])
colsSignalSideband['SevenJet'] = copy.copy(
        colsSignalSideband['MultiJet'])
xbinsSignalSideband['LeptonSevenJet'] = copy.copy(
        xbinsSignalSideband['LeptonMultiJet'])
colsSignalSideband['LeptonSevenJet'] = copy.copy(
        colsSignalSideband['LeptonMultiJet'])

#binning for scale factor histograms
xbinsSignal["TTJetsSingleLepton"] = {}
colsSignal["TTJetsSingleLepton"] = {}
xbinsSignal["TTJetsSingleLepton"]["0B"] = [300, 400, 500, 600, 700, 900, 1200, 4000]
colsSignal["TTJetsSingleLepton"]["0B"] = [
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 1.5 ],
    [ 0.15, 0.2, 0.25, 1.5 ],
    [ 0.15, 0.2, 1.5 ],
    ]

xbinsSignal["WJetsSingleLepton"] = {}
colsSignal["WJetsSingleLepton"] = {}
xbinsSignal["WJetsSingleLepton"]["0B"] = [300, 400, 500, 600, 700, 900, 1200, 4000]
colsSignal["WJetsSingleLepton"]["0B"] = [
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 1.5 ],
    [ 0.15, 0.2, 0.25, 1.5 ],
    ]

xbinsSignal["WJetsSingleLeptonInv"] = {}
colsSignal["WJetsSingleLeptonInv"] = {}
xbinsSignal["WJetsSingleLeptonInv"]["0B"] = [300, 400, 500, 600, 700, 900, 1200, 4000]
colsSignal["WJetsSingleLeptonInv"]["0B"] = [
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.30, 0.41, 1.5 ],
    ]

xbinsSignal["GJetsInv"] = {}
colsSignal["GJetsInv"] = {}
xbinsSignal["GJetsInv"]["0B"] = [400, 500, 600, 700, 900, 1200, 4000]
colsSignal["GJetsInv"]["0B"] = [
    [ 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.25, 0.3, 0.41, 0.52, 1.5 ],
    ]

for box in ['LooseLeptonMultiJet','LooseLeptonDiJet']:
    xbinsSignal[box] = xbinsSignal[box.replace('LooseLepton','')]
    colsSignal[box] = colsSignal[box.replace('LooseLepton','')]

for box in ['DiJet', 'MultiJet', 'SevenJet', 'LeptonMultiJet', 'LeptonSevenJet']:
    rsqMin = 0.3
    if 'Lepton' in box:
        rsqMin = 0.2
    xbinsSignal['Super'+box] = [900, 1600, 4000]
    colsSignal['Super'+box] = [
            [rsqMin, 1.5],
            [0.2, 1.5],
            ]


class Analysis:
    """Class to hold analysis cuts, binning, input files, and trigger info"""
    def __init__(self, region, tag, njetsMin=-1, njetsMax=-1, nbMin=-1, nbMax=-1,
            boostCuts=True):
        self.region = region
        self.tag = tag
        self.njetsMin = njetsMin
        self.njetsMax = njetsMax
        self.nbMin = nbMin
        self.nbMax = nbMax
        self.boostCuts = boostCuts

        if tag == "Razor2015":
            self.lumi = 2300
        elif tag == "Razor2016" or tag == "Razor2016_80X":
            self.lumi = 26400
        elif tag == "Razor2016_ICHEP_80X":
            self.lumi = 12900
        elif tag == "Razor2016G_SUSYUnblind_80X":
            self.lumi = 4394
        elif tag == "Razor2016_MoriondRereco":
            self.lumi = 35867
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
            self.binning = copy.copy(razorBinning["SingleLepton"])
            self.unrollBins = (xbinsSignal["TTJetsSingleLepton"]["0B"], colsSignal["TTJetsSingleLepton"]["0B"])

        elif self.region == "SingleLeptonInv":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoW"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLeptonInv"][self.tag]
            self.samples = razorSamples["WJetsSingleLeptonInv"]
            self.samplesReduced = razorSamplesReduced["WJetsSingleLeptonInv"]
            self.cuts = razorCuts["SingleLeptonInv"]
            self.binning = copy.copy(razorBinning["SingleLeptonInv"])
            self.unrollBins = (xbinsSignal["WJetsSingleLeptonInv"]["0B"], colsSignal["WJetsSingleLeptonInv"]["0B"])

        elif self.region == "TTJetsSingleLepton":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLepton"][self.tag]
            self.samples = razorSamples["TTJetsSingleLepton"]
            self.samplesReduced = razorSamplesReduced["TTJetsSingleLepton"]
            self.cuts = razorCuts["SingleLepton"]
            self.binning = copy.copy(razorBinning["SingleLepton"])
            self.unrollBins = (xbinsSignal["TTJetsSingleLepton"]["0B"], colsSignal["TTJetsSingleLepton"]["0B"])

        elif self.region == "WJetsSingleLepton":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLepton"][self.tag]
            self.samples = razorSamples["WJetsSingleLepton"]
            self.samplesReduced = razorSamplesReduced["WJetsSingleLepton"]
            self.cuts = razorCuts["SingleLepton"]
            self.binning = copy.copy(razorBinning["SingleLepton"])
            self.unrollBins = (xbinsSignal["WJetsSingleLepton"]["0B"], colsSignal["WJetsSingleLepton"]["0B"])

        elif self.region == "WJetsSingleLeptonInv":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoW"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLeptonInv"][self.tag]
            self.samples = razorSamples["WJetsSingleLeptonInv"]
            self.samplesReduced = razorSamplesReduced["WJetsSingleLeptonInv"]
            self.cuts = razorCuts["SingleLeptonInv"]
            self.binning = copy.copy(razorBinning["SingleLeptonInv"])
            self.unrollBins = (xbinsSignal["WJetsSingleLeptonInv"]["0B"], colsSignal["WJetsSingleLeptonInv"]["0B"])

        elif self.region == "TTJetsDilepton":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["Dilepton"][self.tag]
            self.samples = razorSamples["TTJetsDilepton"]
            self.samplesReduced = razorSamplesReduced["TTJetsDilepton"]
            self.cuts = razorCuts["TTJetsDilepton"]
            self.binning = copy.copy(razorBinning["TTJetsDilepton"])
            self.unrollBins = (None,None)

        elif self.region == "TTJetsDileptonMultiJet":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["Dilepton"][self.tag]
            self.samples = razorSamples["TTJetsDilepton"]
            self.samplesReduced = razorSamplesReduced["TTJetsDilepton"]
            self.cuts = razorCuts["TTJetsDilepton"]
            self.binning = copy.copy(razorBinning["TTJetsDileptonReduced"])
            self.unrollBins = (None,None)

        elif self.region == "DYJetsDileptonInv":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoZ"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["DileptonInv"][self.tag]
            self.samples = razorSamples["DYJetsDileptonInv"]
            self.samplesReduced = razorSamplesReduced["DYJetsDileptonInv"]
            self.cuts = razorCuts["DYJetsDileptonInv"]
            self.binning = copy.copy(razorBinning["DYJetsDileptonInv"])
            self.unrollBins = (None,None)

        elif self.region == "DYJetsDileptonInvMultiJet":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoZ"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["DileptonInv"][self.tag]
            self.samples = razorSamples["DYJetsDileptonInv"]
            self.samplesReduced = razorSamplesReduced["DYJetsDileptonInv"]
            self.cuts = razorCuts["DYJetsDileptonInv"]
            self.binning = copy.copy(razorBinning["DYJetsDileptonInvReduced"])
            self.unrollBins = (None,None)

        elif self.region == "GJetsInv":
            self.trigType = "Photon"
            self.jetVar = "NJets_NoPho"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["Photon"][self.tag]
            self.samples = razorSamples["Photon"]
            self.samplesReduced = razorSamplesReduced["Photon"]
            self.cuts = razorCuts["Photon"]
            self.binning = copy.copy(razorBinning["Photon"])
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
            self.binning = copy.copy(razorBinning["SusySync"])
            self.unrollBins = (None,None)

        elif self.region == "VetoLeptonControlRegion":
            self.trigType = "Razor"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["VetoLepton"][self.tag]
            self.samples = razorSamples["VetoLepton"]
            self.samplesReduced = razorSamplesReduced["VetoLepton"]
            self.cuts = razorCuts["VetoLepton"]
            self.binning = copy.copy(razorBinning["VetoLepton"])
            self.unrollBins = (None,None)

        elif (self.region == "MultiJetForVetoLeptonControlRegion"
            or self.region == "SevenJetForVetoLeptonControlRegion"):
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["MultiJetForVetoLepton"]
            self.binning = copy.copy(razorBinning["SignalRegionForVetoLepton"])
            self.unrollBins = (None,None)
        
        elif self.region == "DiJetForVetoLeptonControlRegion":
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["DiJetForVetoLepton"]
            self.binning = copy.copy(razorBinning["SignalRegionForVetoLepton"])
            self.unrollBins = (None,None)

        elif self.region == "VetoTauControlRegion":
            self.trigType = "Razor"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["VetoTau"][self.tag]
            self.samples = razorSamples["VetoTau"]
            self.samplesReduced = razorSamplesReduced["VetoTau"]
            self.cuts = razorCuts["VetoTau"]
            self.binning = copy.copy(razorBinning["VetoTau"])
            self.unrollBins = (None,None)

        elif (self.region == "MultiJetForVetoTauControlRegion"
            or self.region == "SevenJetForVetoTauControlRegion"):
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["MultiJetForVetoTau"]
            self.binning = copy.copy(razorBinning["SignalRegionForVetoTau"])
            self.unrollBins = (None,None)

        elif self.region == "DiJetForVetoTauControlRegion":
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["DiJetForVetoTau"]
            self.binning = copy.copy(razorBinning["SignalRegionForVetoTau"])
            self.unrollBins = (None,None)

        elif self.region in electronRazorBoxes:
            self.trigType = "SingleLepton"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalElectron"][self.tag]
            self.samples = razorSamples["SignalLeptonic"]
            self.samplesReduced = razorSamplesReduced["SignalLeptonic"]
            self.cuts = razorCuts[self.region]
            self.binning = copy.copy(razorBinning["SignalLeptonic"])
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        elif self.region in muonRazorBoxes:
            self.trigType = "SingleLepton"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalMuon"][self.tag]
            self.samples = razorSamples["SignalLeptonic"]
            self.samplesReduced = razorSamplesReduced["SignalLeptonic"]
            self.cuts = razorCuts[self.region]
            self.binning = copy.copy(razorBinning["SignalLeptonic"])
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        elif self.region in leptonRazorBoxes:
            self.trigType = "SingleLepton"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalLepton"][self.tag]
            self.samples = razorSamples["SignalLeptonic"]
            self.samplesReduced = razorSamplesReduced["SignalLeptonic"]
            self.cuts = razorCuts[self.region]
            self.binning = copy.copy(razorBinning[self.region])
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])
            self.boostCuts = False
            print "Note: disabling W/top tag cuts for box {}".format(self.region)

        elif self.region in looseLeptonRazorBoxes:
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts[self.region]
            self.binning = copy.copy(razorBinning[self.region])
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        elif self.region in hadronicRazorBoxes:
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts[self.region]
            self.binning = copy.copy(razorBinning[self.region])
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        elif self.region == "TTJetsSingleLeptonForSignal":
            self.trigType = "SingleLepton"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalLepton"][self.tag]
            self.samples = razorSamples["SignalLeptonic"]
            self.samplesReduced = razorSamplesReduced["SignalLeptonic"]
            self.cuts = razorCuts["SingleLeptonForSignal"]
            self.binning = copy.copy(razorBinning["SingleLeptonForSignal"])
            self.unrollBins = (xbinsSignal["TTJetsSingleLepton"]["0B"], 
                    colsSignal["TTJetsSingleLepton"]["0B"])

        elif self.region == "WJetsSingleLeptonForSignal":
            self.trigType = "SingleLepton"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalLepton"][self.tag]
            self.samples = razorSamples["SignalLeptonic"]
            self.samplesReduced = razorSamplesReduced["SignalLeptonic"]
            self.cuts = razorCuts["SingleLeptonForSignal"]
            self.binning = copy.copy(razorBinning["SingleLeptonForSignal"])
            self.unrollBins = (xbinsSignal["WJetsSingleLepton"]["0B"], 
                    colsSignal["WJetsSingleLepton"]["0B"])

        else: 
            sys.exit("Error: analysis region "+self.region+" is not implemented!")

        #add jet and bjet and trigger requirements to cuts
        #self.cuts = appendNoiseFilters(self.cuts) # MC does not have all noise filters yet -- just add filter cuts to data
        if self.boostCuts:
            self.addBoostCuts()
        self.addJetCuts()
        self.addTriggerCuts()

    def addBoostCuts(self):
        self.cuts += ' && nWTags == 0 && nTopTags == 0'

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
        self.cutsData = appendNoiseFilters( self.trigInfoData.appendTriggerCuts( self.cuts ) )
        if self.region == 'GJetsInv':
            self.cutsData += ' && pho1HLTFilter[36]'

    def getWeightOpts(self):
        self.weightOpts = copy.copy(razorWeightOpts[self.tag])
        if not self.boostCuts and 'boost' in self.weightOpts:
            self.weightOpts.remove('boost')
        self.dataWeightOpts = []
        self.weightFiles = {}
        self.weightHists = {}
        for wType,wInfo in razorWeightHists[self.tag].iteritems():
            self.weightFiles[wType] = rt.TFile.Open(wInfo[0])
            self.weightHists[wType] = self.weightFiles[wType].Get(wInfo[1])

