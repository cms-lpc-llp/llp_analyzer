#Script to create efficiency and acceptance histograms for DY+Jets, W+Jets, and G+Jets 
#Runs on the output of the RazorPhotonStudy analyzer
import sys, os
import string
import ROOT as rt
import array
import copy

def makeEfficiencyHistos(trees, particleType, ptBins, etaBins, tight=False):
    """Creates a (pt, eta) histogram of the efficiency for selecting the given particle type in the given samples"""
    #convert the ptBins and etaBins into arrays, for use with TH2F
    ptBinsArray = array.array('d', ptBins)
    etaBinsArray = array.array('d', etaBins)
    #create numerator and denominator histograms for efficiency calculation
    denominator = rt.TH2F(particleType+"Denominator", particleType+" Denominator", len(ptBins)-1, ptBinsArray, len(etaBins)-1, etaBinsArray)
    numerator = rt.TH2F(particleType+"Efficiency"+(tight*("Tight")), particleType+" Efficiency"+(tight*("Tight")), len(ptBins)-1, ptBinsArray, len(etaBins)-1, etaBinsArray)
    #fill the trees
    for tree in trees:
        tree.Draw("leadingGen"+particleType+"Eta:leadingGen"+particleType+"Pt>>+"+denominator.GetName(), "weight*(nGen"+particleType+"s == 1)") #denominator: all events with one gen-lepton
        tree.Draw("leadingGen"+particleType+"Eta:leadingGen"+particleType+"Pt>>+"+numerator.GetName(), "weight*(nGen"+particleType+"s == 1 && leadGen"+particleType+"IsFound"+((tight)*("Tight"))+")") #numerator: all events with one gen-lepton that is found (loose selection)

    #divide to get the efficiency in each bin
    efficiency = numerator.Clone()
    for i in range(efficiency.GetNbinsX()+1):
        for j in range(efficiency.GetNbinsY()+1):
            if denominator.GetBinContent(i, j) > 0:
                efficiency.SetBinContent(i, j, numerator.GetBinContent(i, j)*1.0/denominator.GetBinContent(i, j))
                errorLow = efficiency.GetBinContent(i, j) - rt.TEfficiency.ClopperPearson(int(denominator.GetBinContent(i, j)), int(numerator.GetBinContent(i, j)), 0.68269, rt.kFALSE)
                errorHigh = rt.TEfficiency.ClopperPearson(int(denominator.GetBinContent(i, j)), int(numerator.GetBinContent(i, j)), 0.68269, rt.kTRUE) - efficiency.GetBinContent(i, j)
                efficiency.SetBinError(i, j, abs((errorHigh - errorLow)/2))
            else:
                efficiency.SetBinContent(i, j, 0.0)
                efficiency.SetBinError(i, j, 0.0)
    return efficiency

def makeAcceptanceHistos(trees, particleType, ptBins, etaBins):
    """Same as makeEfficiencyHistos, but corrects for lepton acceptance (as a function of Z boson pt and eta) instead of efficiency"""
    acceptanceRange = {"Muon":2.4, "Electron":2.5} #max eta for each particle type
    #convert the ptBins and etaBins into arrays, for use with TH2F
    ptBinsArray = array.array('d', ptBins)
    etaBinsArray = array.array('d', etaBins)
    #create numerator and denominator histograms for efficiency calculation
    denominator = rt.TH2F(particleType+"Denominator", particleType+" Denominator", len(ptBins)-1, ptBinsArray, len(etaBins)-1, etaBinsArray)
    numerator = rt.TH2F(particleType+"Acceptance", particleType+" Acceptance", len(ptBins)-1, ptBinsArray, len(etaBins)-1, etaBinsArray)
    #fill the trees
    for tree in trees:
        tree.Draw("genZeta:genZpt>>+"+denominator.GetName(), "weight*(nGen"+particleType+"s == 2 && leadingGen"+particleType+"Pt > 10 && subleadingGen"+particleType+"Pt > 10)") #denominator: all events with two gen-leptons having pt above 10 GeV
        tree.Draw("genZeta:genZpt>>+"+numerator.GetName(), "weight*(nGen"+particleType+"s == 2 && leadingGen"+particleType+"Pt > 10 && subleadingGen"+particleType+"Pt > 10 && abs(leadingGen"+particleType+"Eta) < "+str(acceptanceRange[particleType])+" && abs(subleadingGen"+particleType+"Eta) < "+str(acceptanceRange[particleType])+")") #numerator: all events with two gen-leptons having pt above 10 GeV and eta within the desired acceptance 

    #divide to get the efficiency in each bin
    efficiency = numerator.Clone()
    for i in range(efficiency.GetNbinsX()+1):
        for j in range(efficiency.GetNbinsY()+1):
            if denominator.GetBinContent(i,j) > 0:
                efficiency.SetBinContent(i,j, numerator.GetBinContent(i,j)*1.0/denominator.GetBinContent(i,j))
                errorLow = efficiency.GetBinContent(i,j) - rt.TEfficiency.ClopperPearson(int(denominator.GetBinContent(i,j)), int(numerator.GetBinContent(i,j)), 0.68269, rt.kFALSE)
                errorHigh = rt.TEfficiency.ClopperPearson(int(denominator.GetBinContent(i,j)), int(numerator.GetBinContent(i,j)), 0.68269, rt.kTRUE) - efficiency.GetBinContent(i,j)
                efficiency.SetBinError(i,j, abs((errorHigh - errorLow)/2))
            else:
                efficiency.SetBinContent(i,j, 0.0)
                efficiency.SetBinError(i,j, 0.0)

    return efficiency

##begin main program

#set root to batch mode
rt.gROOT.SetBatch()

#get the files
prefix = ""
postfix = "Run1Eff_19700pb.root"

wnames = ["WJets"]
gnames = ["GJets"]
dynames = ["DYJets"]

wfiles = [rt.TFile(prefix+name+postfix) for name in wnames]
gfiles = [rt.TFile(prefix+name+postfix) for name in gnames]
dyfiles = [rt.TFile(prefix+name+postfix) for name in dynames]

#get the trees
wtrees = [tfile.Get("RazorInclusive") for tfile in wfiles]
gtrees = [tfile.Get("RazorInclusive") for tfile in gfiles]
dytrees = [tfile.Get("RazorInclusive") for tfile in dyfiles]

#set bins in pt and eta for each particle type
muonPtBins = [10, 20, 30, 40, 60, 80, 100, 300, 500, 1000]
electronPtBins = copy.copy(muonPtBins)
photonPtBins = copy.copy(muonPtBins)

muonEtaBins = [0.0, 0.4, 0.7, 1.0, 1.25, 1.479, 1.6, 1.8, 2.0, 2.2, 2.4]
electronEtaBins = copy.copy(muonEtaBins)
electronEtaBins.append(2.5)
photonEtaBins = copy.copy(electronEtaBins)

acceptancePtBins = [0.0, 10, 20, 40, 60, 80, 100, 150, 200, 300, 400, 600, 800, 3000]
acceptanceEtaBins = [0.0, 0.4, 0.8, 1.2, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 4.0, 6.0]

outfile = rt.TFile("Run1LeptonPhotonEfficiency.root", "recreate")
outfile.cd()
results = []
#efficiency for muons
print("Computing efficiency for muons...")
results.append(makeEfficiencyHistos(wtrees, "Muon", muonPtBins, muonEtaBins))
results.append(makeEfficiencyHistos(wtrees, "Muon", muonPtBins, muonEtaBins, tight=True))
#efficiency for photons
print("Computing efficiency for photons...")
results.append(makeEfficiencyHistos(gtrees, "Photon", photonPtBins, photonEtaBins))
#acceptance for muons
print("Creating acceptance histograms for muons...")
results.append(makeAcceptanceHistos(dytrees, "Muon", acceptancePtBins, acceptanceEtaBins))

for result in results: 
    result.Write()
