#Computes photon selection efficiencies (input: the ntuples produced by the photon ntupler in RazorAnalyzer) as a function of pt and eta
import sys, os
import string
import numpy as np
import itertools
import ROOT as rt
import rootTools

##initialization

#check that the input file was specified
if len(sys.argv) < 2:
    print 'Usage: PhotonValidation.py filename.root'
    exit()

#switch ROOT to batch mode
rt.gROOT.SetBatch()
#turn on fit stats
rt.gStyle.SetOptFit(0111)

#get the path of this script and load the razor library
pathname = os.path.dirname(sys.argv[0]) 
fullpath = os.path.abspath(pathname) #get the full path of this script
outpath = fullpath+"/output"
libpath = fullpath+"/lib"
rt.gSystem.Load(libpath+"/libRazorRun2.so")

print('Will compute photon selection efficiencies and validate the photon ID')

#load the TTree from the input file
inFile = rt.TFile(sys.argv[1])
events = inFile.Get("PhotonNtuple")

#drawPhoton =False
drawPhoton =True 
#plot photon efficiency as a function of pT, eta
c = rt.TCanvas("c", "c", 800, 600)
for level in ['Loose', 'Medium', 'Tight']:
    if not drawPhoton: break
    for ptCut in [0, 25, 40, 100]:
        #photon efficiency
        numerator = rt.TH1F("numerator", level+" Photon efficiency vs #eta, p_{T} > "+str(ptCut)+"; #eta; #epsilon", 50, -2.5, 2.5)
        denominator = rt.TH1F("denominator", level+" Photon efficiency vs #eta; #eta; #epsilon", 50, -2.5, 2.5)
        num = events.Draw("pho_superClusterEta>>numerator", "phoMatchesGen && !pho_hasPixelSeed && phoIs"+level+" && phoPt > "+str(ptCut))
        denom = events.Draw("pho_superClusterEta>>denominator", "phoMatchesGen && !pho_hasPixelSeed && phoPt > "+str(ptCut))
        numerator.Divide(denominator)
        numerator.SetStats(0)
        numerator.SetMarkerStyle(20)
        numerator.Draw("pc")
        c.Print(outpath+"/Photon"+level+"EfficiencyVsEtaForPt"+str(ptCut)+".pdf")
        print("\nOverall efficiency for "+level+" photons with pT > "+str(ptCut)+" is "+str(num*1.0/denom)+"\n")
        #check barrel only
        num = events.Draw("", "phoMatchesGen && !pho_hasPixelSeed && phoIs"+level+" && phoPt > "+str(ptCut)+" && abs(pho_superClusterEta) < 1.479")
        denom = events.Draw("", "phoMatchesGen && !pho_hasPixelSeed && phoPt > "+str(ptCut)+" && abs(pho_superClusterEta) < 1.479")
        print("\nOverall efficiency for "+level+" barrel photons with pT > "+str(ptCut)+" is "+str(num*1.0/denom)+"\n")
        #photon energy resolution
        c.SetLogy()
        EnergyResolution = rt.TH1F("EnergyResolution", level+" photon energy resolution in the barrel, p_{T} > "+str(ptCut)+"; #Delta E/E_{gen}", 50, 0., 1)
        events.Draw("deltaEOverEBest>>EnergyResolution", "phoMatchesGen && !pho_hasPixelSeed && phoIs"+level+" && phoPt > "+str(ptCut))
        EnergyResolution.Draw()
        EnergyResolution.Fit("expo")
        c.Print(outpath+"/Photon"+level+"EnergyResolutionFoPt"+str(ptCut)+".pdf")
        c.SetLogy(rt.kFALSE)
        #photon isolation efficiency
        numerator = rt.TH1F("numerator", level+" Photon isolation efficiency vs #eta, p_{T} > "+str(ptCut)+"; #eta; #epsilon", 50, -2.5, 2.5)
        denominator = rt.TH1F("denominator", level+" Photon isolation efficiency vs #eta; #eta; #epsilon", 50, -2.5, 2.5)
        events.Draw("pho_superClusterEta>>numerator", "phoMatchesGen && !pho_hasPixelSeed && phoIsIsolated"+level+" && phoPt > "+str(ptCut))
        events.Draw("pho_superClusterEta>>denominator", "phoMatchesGen && !pho_hasPixelSeed && phoPt > "+str(ptCut))
        numerator.Divide(denominator)
        numerator.SetStats(0)
        numerator.SetMarkerStyle(20)
        numerator.Draw("pc")
        c.Print(outpath+"/Photon"+level+"IsoEfficiencyVsEtaForPt"+str(ptCut)+".pdf")
    for etaCut in [2.5, 1.479]:
        #photon efficiency
        numerator = rt.TH1F("numerator", level+" Photon efficiency vs p_{T}, |#eta| < "+str(etaCut)+"; p_{T}; #epsilon", 50, 0, 400)
        denominator = rt.TH1F("denominator", level+" Photon efficiency vs p_{T}, |#eta| < "+str(etaCut)+"; p_{T}; #epsilon", 50, 0, 400)
        events.Draw("phoPt>>numerator", "phoMatchesGen && !pho_hasPixelSeed && phoIs"+level+" && abs(pho_superClusterEta) < "+str(etaCut))
        events.Draw("phoPt>>denominator", "phoMatchesGen && !pho_hasPixelSeed && abs(pho_superClusterEta) < "+str(etaCut))
        numerator.Divide(denominator)
        numerator.SetStats(0)
        numerator.SetMarkerStyle(20)
        numerator.Draw("pc")
        c.Print(outpath+"/Photon"+level+"EfficiencyVsPtForEta"+str(etaCut)+".pdf")
        #photon isolation efficiency
        numerator = rt.TH1F("numerator", level+" Photon isolation efficiency vs p_{T}, |#eta| < "+str(etaCut)+"; p_{T}; #epsilon", 50, 0, 400)
        denominator = rt.TH1F("denominator", level+" Photon isolation efficiency vs p_{T}, |#eta| < "+str(etaCut)+"; p_{T}; #epsilon", 50, 0, 400)
        events.Draw("phoPt>>numerator", "phoMatchesGen && !pho_hasPixelSeed && phoIsIsolated"+level+" && abs(pho_superClusterEta) < "+str(etaCut))
        events.Draw("phoPt>>denominator", "phoMatchesGen && !pho_hasPixelSeed && abs(pho_superClusterEta) < "+str(etaCut))
        numerator.Divide(denominator)
        numerator.SetStats(0)
        numerator.SetMarkerStyle(20)
        numerator.Draw("pc")
        c.Print(outpath+"/Photon"+level+"IsoEfficiencyVsPtForEta"+str(etaCut)+".pdf")

##make distributions of photon variables with cut levels drawn in
c.SetLogy()
maxes = [0.15, 0.05, 20.0, 30.0, 50.0]
loose = [0.553, 0.0099, 2.49, 15.43, 9.42]
medium = [0.058, 0.0099, 1.91, 4.66, 4.29]
tight = [0.019, 0.0099, 1.61, 3.98, 3.01]
titles = ['H/E', '#sigma_{i#eta i#eta}', 'Charged hadron isolation', 'Neutral hadron isolation', 'Photon energy isolation']
for index, var in enumerate(['pho_HoverE', 'phoFull5x5SigmaIetaIeta', 'phoChHadIsolation', 'phoNeuHadIsolation', 'phoPhotIsolation-phoPt']):
    hist = rt.TH1F("hist", titles[index], 100, 0.0, maxes[index])
    events.Draw(var+">>hist", "phoMatchesGen && !pho_hasPixelSeed && abs(pho_superClusterEta) < 1.479 && phoPt > 40")
    hist.SetStats(0)
    hist.SetLineWidth(2)
    hist.Draw()
    looseLine = rt.TLine(loose[index], 0, loose[index], hist.GetMaximum())
    looseLine.SetLineColor(rt.kGreen)
    looseLine.SetLineWidth(2)
    looseLine.Draw("lsame")
    mediumLine = rt.TLine(medium[index], 0, medium[index], hist.GetMaximum())
    mediumLine.SetLineColor(rt.kOrange)
    mediumLine.SetLineWidth(2)
    mediumLine.Draw("lsame")
    tightLine = rt.TLine(tight[index], 0, tight[index], hist.GetMaximum())
    tightLine.SetLineColor(rt.kRed)
    tightLine.SetLineWidth(2)
    tightLine.Draw("lsame")
    c.Print(outpath+"/Photon"+var+"barrel.pdf")
##plot photon relative isolation to check footprint removal
relHadIso = rt.TH1F("relHadIso", "Relative charged hadron isolation (barrel)", 100, 0.0, 2.0)
relNeuIso = rt.TH1F("relNeuIso", "Relative neutral hadron isolation (barrel)", 100, 0.0, 2.0)
relPhoIso = rt.TH1F("relPhoIso", "Relative photon energy isolation (barrel)", 100, 0.0, 2.0)
events.Draw("phoChHadIsolation/phoPt>>relHadIso", "phoMatchesGen && !pho_hasPixelSeed && abs(pho_superClusterEta) < 1.479")
relHadIso.SetStats(0)
relHadIso.SetLineWidth(2)
relHadIso.Draw()
c.Print(outpath+"/PhotonRelChHadIsolationBarrel.pdf")
events.Draw("phoNeuHadIsolation/phoPt>>relNeuIso", "phoMatchesGen && !pho_hasPixelSeed && abs(pho_superClusterEta) < 1.479")
relNeuIso.SetStats(0)
relNeuIso.SetLineWidth(2)
relNeuIso.Draw()
c.Print(outpath+"/PhotonRelNeuHadIsolationBarrel.pdf")
events.Draw("(phoPhotIsolation-phoPt)/phoPt>>relPhoIso", "phoMatchesGen && !pho_hasPixelSeed && abs(pho_superClusterEta) < 1.479")
relPhoIso.SetStats(0)
relPhoIso.SetLineWidth(2)
relPhoIso.Draw()
c.Print(outpath+"/PhotonRelPhoHadIsolationBarrel.pdf")
