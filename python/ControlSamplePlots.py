#Plotting script to compare kinematic distributions in Z->nunu with those obtained from the DY+jets, W+jets, and photon+jets control samples
import sys, os
import string
import ROOT as rt

def compareControlDistributions(trees, quantitiesToPlot, titles, xmins, xmaxs, nBins, effHistogram = None, accHistogram = None, reweighBy = None, reweighByHistogram=None, selectByGenMuons=False):
    """Draws histograms of the desired quantity for the Z->invisible samples, the DY samples, the W samples, and the Gamma samples, applies formatting, and prints the plot as a PDF.  Supply reweighByHistogram (it should be a dict of histograms returned by a different invocation of this function) to reweigh the measured distributions according to the given histograms, using the ZJets distribution as reference."""
    if reweighBy is not None and reweighByHistogram is None:
        print("Error in compareControlDistributions: reweighing expression provided without accompanying histograms to reweigh by!")
        return
    elif reweighBy is None and reweighByHistogram is not None:
        print("Error in compareControlDistributions: reweighing histograms provided without accompanying variable name!")
        return

    #names of samples to draw
    sampleNames = ["ZJets", "DYJets", "WJets", "GJets"]
    normFactors = {"ZJets": 1.0, "DYJets": 0.03366/0.2, "WJets": 2363.02/489.503*10.57/32.57, "GJets": 1.0}

    #legend entries
    sampleTitles = {}
    sampleTitles["ZJets"] = "Z -> #nu#nu + Jets"
    sampleTitles["DYJets"] = "DY -> ll + Jets"
    sampleTitles["WJets"] = "W -> l#nu + Jets"
    sampleTitles["GJets"] = "#gamma + Jets"

    #postfixes for each sample
    postfixes = {}
    postfixes["ZJets"] = ""
    if selectByGenMuons: #select events with two gen muons
        postfixes["DYJets"] = "_noGenZ"
        postfixes["WJets"] = "_noGenZ"
    else: #select events with two reco muons
        postfixes["DYJets"] = "_noZ"
        postfixes["WJets"] = "_noW"
    postfixes["GJets"] = "_noPho"

    #efficiency expressions 
    effPtExps = {}
    effEtaExps = {}
    if effHistogram is not None:
        effPtExps["DYJets"] = "leadingMuonPt"
        effEtaExps["DYJets"] = "abs(leadingMuonEta)"
        effPtExps["WJets"] = "leadingMuonPt"
        effEtaExps["WJets"] = "abs(leadingMuonEta)"
        effPtExps["GJets"] = "leadingPhotonPt"
        effEtaExps["GJets"] = "abs(leadingPhotonEta)"
        effPtExps2 = effPtExps["DYJets"].replace("leading", "subleading") #for DY+Jets
        effEtaExps2 = effEtaExps["DYJets"].replace("leading", "subleading") #for DY+Jets

    #acceptance expressions (pt and eta of the Z boson, used to correct for lepton acceptance in DY+Jets)
    accPtExps = {}
    accEtaExps = {}
    if accHistogram is not None:
        accPtExps["DYJets"] = "recoZpt"
        accEtaExps["DYJets"] = "recoZeta"

    #expression to reweigh by 
    reweighExps = {}
    if reweighBy is not None:
        for s in sampleNames: reweighExps[s] = reweighBy+postfixes[s]

    #baseline selection for each sample
#    conditions = {}
#    conditions["ZJets"] = "nLooseMuons == 0 && nLooseElectrons == 0 && numJets80 >= 2 && MR > 300 && Rsq > 0.15"
#    if selectByGenMuons:
#        conditions["DYJets"] = "nGenMuons == 2 && numJets80"+postfixes["DYJets"]+" >= 2 && MR"+postfixes["DYJets"]+" > 300 && Rsq"+postfixes["DYJets"]+" > 0.15 && genZmass > 71 && genZmass < 111"
#    else:
#        conditions["DYJets"] = "nLooseMuons == 2 && numJets80"+postfixes["DYJets"]+" >= 2 && MR"+postfixes["DYJets"]+" > 300 && Rsq"+postfixes["DYJets"]+" > 0.15 && recoZmass > 71 && recoZmass < 111"
#    conditions["WJets"] = "nTightMuons == 1 && nLooseMuons == 1 && numJets80_noW >= 2 && MR_noW > 300 && Rsq_noW > 0.15 && mTLepMet > 30 && mTLepMet < 100"
#    conditions["GJets"] = "numJets80_noPho >= 2 && MR_noPho > 300 && Rsq_noPho > 0.15 && nSelectedPhotons >= 1"
    conditions = {}
    conditions["ZJets"] = "nGenNeutrinos == 2 && MR > 300 && Rsq > 0.15 && numJets80 >= 2"
    conditions["DYJets"] = "nGenMuons == 2 && genZmass > 71 && genZmass < 111"
    conditions["WJets"] = "nGenMuons == 1"
    conditions["GJets"] = "nGenPhotons >= 1"

    #histogram colors
    colors = {}
    colors["ZJets"] = rt.kCyan+3
    colors["DYJets"] = rt.kAzure
    colors["WJets"] = rt.kOrange+10
    colors["GJets"] = rt.kTeal+10

    #histograms, canvas and legend
    histos = {}
    for name in sampleNames: 
        histos[name] = {}
        for i, q in enumerate(quantitiesToPlot): histos[name][q] = rt.TH1F(name+q+((reweighBy is not None)*("Reweigh"))+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr")), titles[i]+"; "+titles[i], nBins[i], xmins[i], xmaxs[i])
    c = rt.TCanvas("c", "c", 800, 600)
    leg = rt.TLegend(0.7, 0.7, 0.9, 0.9)
    legNoGamma = rt.TLegend(0.7, 0.7, 0.9, 0.9)

    #loop over samples and fill histograms
    for treeName in trees:
        for name in sampleNames:
            if name == "GJets" or name == "ZJets": continue #save gamma for later
            if name in treeName: #fill appropriate histogram, with weight
                if reweighByHistogram is None and effHistogram is None and accHistogram is None: #reweigh only by sample cross section
                    print("Filling "+name+" from tree "+treeName+" (weighing by sample cross section)")
                    for q in quantitiesToPlot:
                        trees[treeName].Draw(q+postfixes[name]+">>+"+histos[name][q].GetName(), "weight*("+conditions[name]+")*1.0/"+str(normFactors[name]))
                else: #reweigh by sample cross section and according to reweighByHistogram and/or effHistogram
                    print("Filling "+name+" from tree "+treeName+" (weighing by sample cross section "+((reweighBy is not None)*("and "+str(reweighBy)))+((effHistogram is not None)*(" and selection efficiency"))+")")

                    #formula to fill
                    form = {}
                    for q in quantitiesToPlot:
                        form[q] = rt.TTreeFormula(name+q+"Formula", q+postfixes[name], trees[treeName])
                        form[q].GetNdata()

                    #formulas for reweighing
                    reweighForm = None
                    if reweighByHistogram is not None: 
                        reweighForm = rt.TTreeFormula(name+"ReweighFormula", reweighExps[name], trees[treeName])
                        reweighForm.GetNdata()
                    #formulas for reweighing based on pt and eta
                    effPtForm = None
                    effEtaForm = None
                    effPtForm2 = None
                    effEtaForm2 = None
                    accPtForm = None
                    accEtaForm = None
                    if effHistogram is not None and name != "ZJets":
                        effPtForm = rt.TTreeFormula(name+"PtEfficiencyFormula", effPtExps[name], trees[treeName])
                        effEtaForm = rt.TTreeFormula(name+"EtaEfficiencyFormula", effEtaExps[name], trees[treeName])
                        effPtForm.GetNdata()
                        effEtaForm.GetNdata()
                        if name == "DYJets": 
                            effPtForm2 = rt.TTreeFormula(name+"PtEfficiencyFormula2", effPtExps2, trees[treeName]) #this one corrects for the efficiency of the second lepton in DY+Jets events
                            effEtaForm2 = rt.TTreeFormula(name+"EtaEfficiencyFormula2", effEtaExps2, trees[treeName]) #this one corrects for the efficiency of the second lepton in DY+Jets events
                            effPtForm2.GetNdata()
                            effEtaForm2.GetNdata()
                        #get the correct histogram for the efficiency
                        thisEffHisto = effHistogram[name]
                    if accHistogram is not None and name == "DYJets":
                        accPtForm = rt.TTreeFormula(name+"ZPtFormula", accPtExps[name], trees[treeName])
                        accEtaForm = rt.TTreeFormula(name+"ZEtaFormula", accEtaExps[name], trees[treeName])
                        accPtForm.GetNdata()
                        accEtaForm.GetNdata()
                        thisAccHisto = accHistogram[name]

                    #condition for event selection 
                    conditionsForm = rt.TTreeFormula(name+"ConditionsFormula", conditions[name], trees[treeName])
                    conditionsForm.GetNdata()

                    nEvents = trees[treeName].GetEntries()
                    for n, event in enumerate(trees[treeName]):
                        #count events
                        if n % 1000000 == 0: print("Processing event "+str(n)+" of "+str(nEvents))

                        #check if the event passes
                        passesCondition = conditionsForm.EvalInstance()
                        if not passesCondition: continue
                        
                        #value of the quantity we are filling
                        exprValue = {}
                        for q in quantitiesToPlot: exprValue[q] = form[q].EvalInstance()

                        #initialize weight
                        weight = event.weight
                        weight /= normFactors[name] #normalize to Z->nunu cross section
                        #value of the quantity whose distribution should be made to match that of Z+jets
                        if reweighForm is not None: #compute contribution to the weight from the reweighing histogram
                            reweighExprValue = reweighForm.EvalInstance() 
                            #multiply weight by (value of reweighing histogram in Z+jets)/(value of reweighing histogram in the given sample)
                            reweighRatio = reweighByHistogram["DYJets"].GetBinContent(reweighByHistogram["DYJets"].FindBin(reweighExprValue))*1.0/reweighByHistogram[name].GetBinContent(reweighByHistogram[name].FindBin(reweighExprValue))
                            weight = weight*reweighRatio

                        #get appropriate pt's and etas for reweighing
                        if effHistogram is not None and name != "ZJets":
                            effPtExprValue = effPtForm.EvalInstance()
                            effEtaExprValue = effEtaForm.EvalInstance()
                            if name == "DYJets":
                                effPtExprValue2 = effPtForm2.EvalInstance()
                                effEtaExprValue2 = effEtaForm2.EvalInstance()
                        if accHistogram is not None and name == "DYJets":
                            accPtExprValue = accPtForm.EvalInstance()
                            accEtaExprValue = accEtaForm.EvalInstance()

                        #reweigh according to selection efficiency
                        if effHistogram is not None and name != "ZJets":
                            ptMax = 1000
                            if effPtExprValue > ptMax: effPtExprValue = ptMax - 1 #if pt is higher than the maximum of the histogram, use the last pt bin
                            if name == "DYJets" and effPtExprValue2 > ptMax: effPtExprValue2 = ptMax - 1
                            efficiency = thisEffHisto.GetBinContent(thisEffHisto.FindBin(effPtExprValue, effEtaExprValue))
                            if name == "DYJets": efficiency *= thisEffHisto.GetBinContent(thisEffHisto.FindBin(effPtExprValue2, effEtaExprValue2))

                            weight = weight/efficiency

                        #reweigh according to acceptance
                        if accHistogram is not None and name == "DYJets":
                            ptMax = 3000
                            etaMax = 6.0
                            if accPtExprValue > ptMax: effPtExprValue = ptMax - 1 
                            if abs(accEtaExprValue) > etaMax: 
                                print("Warning: Z boson pt is outside of histogram range: "+str(accEtaExprValue))
                                accEtaExprValue = etaMax - 0.1
                            accEfficiency = thisAccHisto.GetBinContent(thisAccHisto.FindBin(accPtExprValue, accEtaExprValue))
                            if accEfficiency > 0.001:
                                weight = weight/accEfficiency
                            else:
                                print("Warning: trying to divide by zero in acceptance correction calculation!")
                                weight = 0

                        for q in quantitiesToPlot: histos[name][q].Fill(exprValue[q], weight)
    #end of loop over samples

    #format histos
    for name in sampleNames:
        for q in quantitiesToPlot:
            histos[name][q].SetLineWidth(3)
            histos[name][q].SetStats(0)
            histos[name][q].SetLineColor(colors[name])
            histos[name][q].SetTitleSize(0.3,"t"); 
            print("Integral is "+str(histos[name][q].Integral()))

    #add histos to legend
    for name in sampleNames: 
        #add to the legend (the itervalues.next expression just gives you a histogram from the dict
        if name != "ZJets": leg.AddEntry(histos[name].itervalues().next(), sampleTitles[name])
        if name != "GJets" and name != "ZJets": legNoGamma.AddEntry(histos[name].itervalues().next(), sampleTitles[name])

    #fill stacked histogram with the histograms from each sample
    #fill another stacked histogram with each distribution normalized to Z+Jets
    for i, q in enumerate(quantitiesToPlot):
        stack = rt.THStack("stack", titles[i]) #for drawing histograms
        stackNoGamma = rt.THStack("stackNoGamma", titles[i]) #for drawing histograms
        ratiostack = rt.THStack("ratiostack", "")
        ratiostackNoGamma = rt.THStack("ratiostackNoGamma", "")
        ratioHistos = {}
        for name in sampleNames:
            if name != "ZJets": stack.Add(histos[name][q])
            ratioHistos[name] = histos[name][q].Clone()
            ratioHistos[name].Divide(histos["DYJets"][q])
            ratioHistos[name].SetLineWidth(2)
            if name != "ZJets": ratiostack.Add(ratioHistos[name])
            if name != "GJets" and name != "ZJets": 
                stackNoGamma.Add(histos[name][q])
                ratiostackNoGamma.Add(ratioHistos[name])

        #draw and print all histograms
        c.Clear()
        c.cd()
        pad1 = rt.TPad("pad1","pad1",0,0.4,1,1)
        pad1.SetBottomMargin(0)
        pad1.SetLogy()
        pad1.Draw()
        pad1.cd()
        stackNoGamma.Draw("nostack,elp")
        pad1.Modified()
        rt.gPad.Update()
        stackNoGamma.GetXaxis().SetTitle(str(titles[i]))
        stackNoGamma.GetYaxis().SetTitle("Number of events in 4/fb")
        stackNoGamma.GetYaxis().SetLabelSize(0.03)
        stackNoGamma.GetYaxis().SetTitleOffset(0.45)
        stackNoGamma.GetYaxis().SetTitleSize(0.05)
        legNoGamma.Draw()
        c.cd()
        pad2 = rt.TPad("pad2","pad2",0,0.0,1,0.4)
        pad2.SetTopMargin(0)
        pad2.SetTopMargin(0.008)
        pad2.SetBottomMargin(0.25)
        pad2.SetGridy()
        pad2.Draw()
        pad2.cd()
        ratiostackNoGamma.Draw("nostack,elp")
        pad2.Modified()
        rt.gPad.Update()
        ratiostackNoGamma.SetMinimum(0.7)
        ratiostackNoGamma.SetMaximum(1.3)
        ratiostackNoGamma.GetXaxis().SetTitle(str(titles[i]))
        ratiostackNoGamma.GetYaxis().SetTitle("Ratio with Z->#nu#nu")
        ratiostackNoGamma.GetXaxis().SetLabelSize(0.1);
        ratiostackNoGamma.GetYaxis().SetLabelSize(0.08);
        ratiostackNoGamma.GetYaxis().SetTitleOffset(0.35);
        ratiostackNoGamma.GetXaxis().SetTitleOffset(1.00);
        ratiostackNoGamma.GetYaxis().SetTitleSize(0.08);
        ratiostackNoGamma.GetXaxis().SetTitleSize(0.08);

        c.Print("controlSample"+q+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+(selectByGenMuons*("SelectByGenMuons"))+".pdf")
        c.Print("controlSample"+q+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+(selectByGenMuons*("SelectByGenMuons"))+".root")

    #return the histograms
    return histos

##begin main program

#set ROOT to batch mode
rt.gROOT.SetBatch()

#load the efficiency histogram file; load the acceptance and efficiency histograms
efficiencyFile = rt.TFile("Phys14LeptonPhotonEfficiencyNoteIncorrectErrors.root")
effHistos = {}
accHistos = {}
effHistos["DYJets"] = efficiencyFile.Get("WJetsMuonEfficiency")
effHistos["WJets"] = efficiencyFile.Get("WJetsMuonEfficiencyTight")
effHistos["GJets"] = efficiencyFile.Get("GJetsPhotonEfficiency")
accHistos["DYJets"] = efficiencyFile.Get("DYJetsMuonAcceptance")

#load TFiles and initialize trees
datanames = [
             "DYJets100", "DYJets200", "DYJets400", "DYJets600", 
             "ZJets100", "ZJets200", "ZJets400", "ZJets600", 
             "WJets100", "WJets200", "WJets400", "WJets600", 
             "GJets100", "GJets200", "GJets400", "GJets600"
             ]
prefix = "control"
postfix = "_4000pb_weighted.root"
files = {}
for sample in datanames: files[sample] = rt.TFile(prefix+sample+postfix)
trees = {}
treeName = "RazorInclusive"
for sample in datanames: trees[sample] = files[sample].Get(treeName)

#file for saving histograms
outfile = rt.TFile("controlSampleHistograms.root", "recreate")
outfile.cd()

#plots: MET, Rsq, MR
results = {}
reweighingHistos = {}

#without reweighing
results["histosNoReweighingGenMuons"] = compareControlDistributions(trees, ["met", "Rsq", "MR", "HT", "numJets"], ["MET (GeV)", "R^{2}", "M_{R} (GeV)", "HT (GeV)", "Number of jets"], [0., 0.15, 300, 0., 0], [1000, 1, 1500, 3000, 11], [25, 25, 25, 25, 11], selectByGenMuons=True)
for name in results["histosNoReweighingGenMuons"]: reweighingHistos[name] = results["histosNoReweighingGenMuons"][name]["met"] #set up reweighing histogram using the latest MET distribution
results["histosMetCorrGenMuons"] = compareControlDistributions(trees, ["met", "Rsq", "MR", "HT", "numJets"], ["MET (GeV)", "R^{2}", "M_{R} (GeV)", "HT (GeV)", "Number of jets"], [0.0, 0.15, 300, 0.0, 0], [1000, 1, 1500, 3000, 11], [25, 25, 25, 25, 11], reweighByHistogram = reweighingHistos, reweighBy = "met", selectByGenMuons=True)

#results["histosNoReweighing"] = compareControlDistributions(trees, ["met", "Rsq", "MR", "HT", "numJets"], ["MET (GeV)", "R^{2}", "M_{R} (GeV)", "HT (GeV)", "Number of jets"], [0.0, 0.0, 0.0, 0.0, 0], [1000, 1, 1500, 3000, 11], [25, 25, 25, 25, 11])
#reweigh by lepton/photon selection efficiency
#results["histosEffCorr"] = compareControlDistributions(trees, ["met", "Rsq", "MR", "HT", "numJets"], ["MET (GeV)", "R^{2}", "M_{R} (GeV)", "HT (GeV)", "Number of jets"], [0.0, 0.0, 0.0, 0.0, 0], [1000, 1, 1500, 3000, 11], [25, 25, 25, 25, 11], effHistogram = effHistos)
#correct for lepton acceptance
#results["histosEffCorrAccCorr"] = compareControlDistributions(trees, ["met", "Rsq", "MR", "HT", "numJets"], ["MET (GeV)", "R^{2}", "M_{R} (GeV)", "HT (GeV)", "Number of jets"], [0.0, 0.0, 0.0, 0.0, 0], [1000, 1, 1500, 3000, 11], [25, 25, 25, 25, 11], effHistogram = effHistos, accHistogram = accHistos)
#after reweighing by MET
#for name in results["histosEffCorrAccCorr"]: reweighingHistos[name] = results["histosEffCorrAccCorr"][name]["met"] #set up reweighing histogram using the latest MET distribution
#results["histosEffCorrAccCorrMetCorr"] = compareControlDistributions(trees, ["met", "Rsq", "MR", "HT", "numJets"], ["MET (GeV)", "R^{2}", "M_{R} (GeV)", "HT (GeV)", "Number of jets"], [0.0, 0.0, 0.0, 0.0, 0], [1000, 1, 1500, 3000, 11], [25, 25, 25, 25, 11], effHistogram = effHistos, accHistogram = accHistos, reweighBy = "met", reweighByHistogram = reweighingHistos)

#save histograms to ROOT file
for result in results:
    for key in results[result]:
        for q in results[result][key]:
            results[result][key][q].Write()
