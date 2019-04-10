import sys, os, copy
import ROOT as rt

from macro import macro, razorWeights
from macro.razorAnalysis import Analysis, make_parser
from macro.razorMacros import makeControlSampleHistsForAnalysis, makeVetoLeptonCorrectionHist
import BTagClosureTestMacro as bclosure

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    parser = make_parser()
    parser.add_argument("--muons", help="require muons", action='store_true')
    parser.add_argument("--electrons", help="require electrons", action='store_true')
    parser.add_argument("--tight", help="require tight leptons", action='store_true')
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    boostCuts = not args.noBoostCuts

    #load the MT cut efficiency as a function of lepton pt
    mtHists = {}
    mtLepFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/VetoLeptonMTCutEfficiency.root")
    mtHists["VetoLeptonPt"] = mtLepFile.Get("VetoLeptonMTCutEfficiencyVsPt")
    mtHists["VetoLeptonEta"] = mtLepFile.Get("VetoLeptonMTCutEfficiencyVsEta")
    #load the MT cut efficiency as a function of tau pt
    mtTauFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/VetoTauMTCutEfficiency.root")
    mtHists["VetoTauPt"] = mtTauFile.Get("VetoTauMTCutEfficiencyVsPt")
    mtHists["VetoTauEta"] = mtTauFile.Get("VetoTauMTCutEfficiencyVsEta")
    #load the dPhi cut efficiency as a function of lepton pt
    dphiHists = {}
    dphiLepFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/DPhiCutEfficiencyForLostLepton.root")
    dphiHists["VetoLeptonPt"] = dphiLepFile.Get("VetoLeptonDPhiCutEfficiencyVsPt")
    dphiHists["VetoLeptonEta"] = dphiLepFile.Get("VetoLeptonDPhiCutEfficiencyVsEta")
    #load the dPhi cut efficiency as a function of tau pt
    dphiTauFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/DPhiCutEfficiencyForLostTau.root")
    dphiHists["VetoTauPt"] = dphiTauFile.Get("VetoTauDPhiCutEfficiencyVsPt")
    dphiHists["VetoTauEta"] = dphiTauFile.Get("VetoTauDPhiCutEfficiencyVsEta")

    #initialize
    plotOpts = { "comment":False, 'SUS15004CR':True } 
    regionsOrder = []
    regions = {}
    regionsCorrespondence = {} #get the control region corresponding to each signal region
    regionMtHists = {}
    regionDphiHists = {}
    for ltype in ["VetoLepton","VetoTau"]:
        for jtype,jets in {
                "DiJet":(2,3),
                "MultiJet":(4,6),
                "SevenJet":(7,-1)
                }.iteritems():
            #veto lepton/tau control region
            regionsOrder.append(ltype+jtype) 
            regions[ltype+jtype] = Analysis(ltype+"ControlRegion",tag=tag,
                    njetsMin=jets[0],njetsMax=jets[1],boostCuts=boostCuts)
            #corresponding signal region 
            sigR = jtype+"For"+ltype
            regionsOrder.append(sigR) 
            regions[jtype+"For"+ltype] = Analysis(sigR+"ControlRegion",
                    tag=tag,njetsMin=jets[0],njetsMax=jets[1],boostCuts=boostCuts)
            regionsCorrespondence[sigR] = ltype+jtype
            regionMtHists[sigR] = mtHists[ltype+"Pt"]
            regionDphiHists[sigR] = dphiHists[ltype+"Pt"]
            #veto lepton/tau control region
            regionsOrder.append(ltype+jtype+"PtCorr") 
            regions[ltype+jtype+"PtCorr"] = Analysis(ltype+"ControlRegion",tag=tag,
                    njetsMin=jets[0],njetsMax=jets[1],boostCuts=boostCuts)
            #corresponding signal region 
            sigRPtCorr = jtype+"For"+ltype+"PtCorr"
            regionsOrder.append(sigRPtCorr) 
            regions[sigRPtCorr] = Analysis(sigR+"ControlRegion",
                    tag=tag,njetsMin=jets[0],njetsMax=jets[1],boostCuts=boostCuts)
            regionsCorrespondence[sigRPtCorr] = ltype+jtype+"PtCorr"
            regionMtHists[sigRPtCorr] = mtHists[ltype+"Eta"]
            regionDphiHists[sigRPtCorr] = dphiHists[ltype+"Eta"]

    sfHists = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), 
            processNames=regions["VetoLeptonSevenJet"].samples, scaleFactorNames={ "ZInv":"GJetsInv" },
            debugLevel=debugLevel)
    sfHistsSignal = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag),
            processNames=regions["SevenJetForVetoLepton"].samples, scaleFactorNames={ "ZInv":"GJetsInv", "TTJets1L":"TTJets", "TTJets2L":"TTJets" }, debugLevel=debugLevel)
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    for h in [sfHists, sfHistsSignal]:
        h['NJetsTTJets'] = sfNJetsFile.Get("TTJetsScaleFactors")
        h['NJetsWJets'] = sfNJetsFile.Get("WJetsScaleFactors")
        h['NJetsInv'] = sfNJetsFile.Get("GJetsInvScaleFactors")
        bclosure.loadScaleFactors(h, tag=tag)
    sfVars = ("MR","Rsq")
    if not args.noSave:
        #recreate output file to avoid confusion
        outfile = rt.TFile("data/ScaleFactors/RazorMADD2015/RazorVetoLeptonClosureTests_%s.root"%(tag), "RECREATE")
        outfile.Close()

    hists = {}
    for region in regionsOrder:
        analysis = regions[region]
        outdir = 'Plots/'+tag+'/'+region
        if args.tight:
            outdir = outdir.replace("Veto","Tight")
            analysis.samples.remove('ZInv')
            analysis.samples.remove('QCD')
            analysis.cutsData += " && lep1PassTight && ((abs(lep1Type) == 11 && lep1.Pt() > 30) || (abs(lep1Type) == 13 && lep1.Pt() > 25)) && MET > 30"
            analysis.cutsData = analysis.cutsData.replace('&& NJets80 >= 2','')
            analysis.cutsMC += " && lep1PassTight && ((abs(lep1Type) == 11 && lep1.Pt() > 30) || (abs(lep1Type) == 13 && lep1.Pt() > 25)) && MET > 30"
            analysis.cutsMC = analysis.cutsMC.replace('&& NJets80 >= 2','')
        #modify cuts and output dir if running with ele/mu separately
        if args.electrons:
            analysis.cutsData += " && abs(lep1Type) == 11"
            analysis.cutsMC += " && abs(lep1Type) == 11"
            outdir = outdir.replace('Lepton','Electron')
            if 'Tau' in region: continue
        if args.muons:
            analysis.cutsData += " && abs(lep1Type) == 13"
            analysis.cutsMC += " && abs(lep1Type) == 13"
            outdir = outdir.replace('Lepton','Muon')
            if 'Tau' in region: continue
        os.system('mkdir -p '+outdir)
        #set up analysis
        if region.startswith('Veto'):
            sfHistsToUse = sfHists
            auxSFsToUse = razorWeights.getNJetsSFs(analysis,jetName='NJets40')
            treeName = "ControlSampleEvent"
            bjetsName = 'NBJetsMedium'
        else:
            sfHistsToUse = sfHistsSignal
            auxSFsToUse = razorWeights.getNJetsSFs(analysis,jetName='nSelectedJets')
            treeName = "RazorInclusive"
            bjetsName = 'nBTaggedJets'
        auxSFsToUse = razorWeights.addAllBTagSFs(analysis, auxSFsToUse, bjetsName=bjetsName)
        bclosure.adjustForRegionBInclusive(analysis, sfHistsToUse, auxSFsToUse)
        #set up lepton pt correction
        varForCorrection = "lep1.Pt()"
        # The non-closure on the lepton pt distribution is taken as a 
        # correlated shape systematic. For lepton eta, on the other hand,
        # in the absence of any systematic non-closure we take the
        # uncertainty on the data/MC ratio as a systematic in each bin
        useUncertainty = False 
        sigVarForCorrection = "leadingGenLeptonPt"
        if "PtCorr" in region: 
            varForCorrection = "abs(lep1.Eta())"
            useUncertainty = True
            sigVarForCorrection = "abs(leadingGenLeptonEta)"
            if region.startswith('VetoLepton'):
                for proc in analysis.samples:
                    auxSFsToUse[proc][region.replace('PtCorr','')] = ("lep1.Pt()", 
                            "abs(lep1Type) == 11 || abs(lep1Type) == 13")
            elif region.startswith('VetoTau'):
                for proc in analysis.samples:
                    auxSFsToUse[proc][region.replace('PtCorr','')] = ("lep1.Pt()", "abs(lep1Type) == 15")
            elif 'VetoLepton' in region:
                for proc in analysis.samples:
                    auxSFsToUse[proc][region.replace('PtCorr','')] = ("leadingGenLeptonPt", 
                            "abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13")
            elif 'VetoTau' in region:
                for proc in analysis.samples:
                    auxSFsToUse[proc][region.replace('PtCorr','')] = ("leadingGenLeptonPt", "abs(leadingGenLeptonType) == 15")
        #sanity check
        print "\nRegion:",region
        print "Tree name:",treeName
        print "Aux SFs to use:",auxSFsToUse
        print "Variable for correction:",varForCorrection
        print "Signal region variable for correction:",sigVarForCorrection,"\n"
        #perform analysis
        hists[region] = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, 
                sfHists=sfHistsToUse, sfVars=sfVars, printdir=outdir, auxSFs=auxSFsToUse, 
                treeName=treeName, debugLevel=debugLevel, noFill=args.noFill )
        if not args.noSave:
            #export histograms
            macro.exportHists(hists[region], outFileName='controlHistograms'+region+'.root', 
                    outDir=outdir, debugLevel=debugLevel, delete=False)
        #compute correction factors
        if region.startswith('Veto'):
            #make control region scale factors
            sfHists[region] = makeVetoLeptonCorrectionHist(hists[region], lumiData=analysis.lumi, 
                    debugLevel=debugLevel, var=varForCorrection, 
                    regionName=region, doDataOverMC=True, sfHists=sfHistsToUse, 
                    printdir=outdir)
            sfHistsSignal[region] = sfHists[region]
        else:
            #make signal region scale factors
            controlRegionHists = hists[regionsCorrespondence[region]]
            mtHistToUse = regionMtHists[region]
            dphiHistToUse = regionDphiHists[region]
            sfHists[region] = makeVetoLeptonCorrectionHist(controlRegionHists, 
                    lumiData=analysis.lumi, debugLevel=debugLevel, var=varForCorrection, 
                    regionName=region, doDataOverMC=False, useUncertainty=useUncertainty,
                    histsToCorrect=hists[region], signalRegionVar=sigVarForCorrection, 
                    mtEfficiencyHist=mtHistToUse, dPhiEfficiencyHist=dphiHistToUse, 
                    sfHists=sfHistsToUse, printdir=outdir)
            sfHistsSignal[region] = sfHists[region]
            #write out to file
            sfHistClone = sfHists[region].Clone()
            if not args.noSave:
                print "Writing scale factor histogram",sfHistClone.GetName(),"to file"
                outfile = rt.TFile("data/ScaleFactors/RazorMADD2015/RazorVetoLeptonClosureTests_%s.root"%(tag), "UPDATE")
                outfile.Append(sfHistClone)
                outfile.Write()
                outfile.Close()

