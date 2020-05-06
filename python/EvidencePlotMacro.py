import sys,os
import copy
import argparse
import ROOT as rt
import glob

#local imports
from macro.razorAnalysis import xbinsSignal, colsSignal
from macro.razorMacros import makeRazorBinEvidencePlots, makeRazorMCTotalPlusSignalPlot
from RunCombine import exec_me
import macro.macro as macro
import macro.razorAnalysis as razor
import SignalRegionMacro as signal
import SMSTemplates as sms

BACKGROUND_DIR = "/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorMADD2016_08Apr2018"

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('--box', help="choose a box", required=True)
    parser.add_argument('--dir', help="output directory", default="SignalRegionPlots", dest='outDir')
    parser.add_argument('-m','--model', default="T1bbbb", help="signal model name")
    parser.add_argument('--mGluino',default=-1,type=int, help="mass of gluino")
    parser.add_argument('--mStop',default=-1,type=int, help="mass of stop")
    parser.add_argument('--mLSP',default=-1,type=int, help="mass of LSP")

    args = parser.parse_args()
    outDir = args.outDir
    box = args.box
    boxList = box.split('_') #interpret box1_box2_... as a list of boxes to combine
    tag = 'Razor2016_MoriondRereco'

    dirToUse = razor.razorSignalDirs[tag]

    #make output directory
    os.system('mkdir -p '+outDir)

    c = rt.TCanvas("c","c",800,600)
    for curBox in boxList:
        #get configuration for this box
        nbmax = 3
        if curBox == 'DiJet':
            nbmax = 2
        analyses = [razor.Analysis(curBox, tag=tag,
            nbMin=nb) for nb in range(nbmax+1)]
        unrollBins = [a.unrollBins for a in analyses]
        uncerts = copy.copy(sms.signalShapeUncerts)
        signal.adjustForFineGrainedMCPred(analyses[0], {},
                {}, {}, {})
        samples = analyses[0].samples

        #call SMS template maker
        if 'T1' in args.model or 'T5' in args.model:
            modelName = 'SMS-'+args.model+'_'+str(args.mGluino)+'_'+str(args.mLSP)
        else:
            modelName = 'SMS-'+args.model+'_'+str(args.mStop)+'_'+str(args.mLSP)
        signalFilename=dirToUse+'/'+modelName+'.root'
        brString = ""
        xBR=-1
        yBR=-1
        if 'T1x' in args.model:
            xBR = float(args.model[args.model.find('x')+1:args.model.find('y')].replace('p','.'))
            yBR = float(args.model[args.model.find('y')+1:].replace('p','.'))
            brString = '--xBR %.2f --yBR %.2f'%(xBR,yBR)
            signalFilename = dirToUse+'/SMS-T1ttbb_'+str(args.mGluino)+'_'+str(args.mLSP)+'.root'
        smsOpts = sms.SMSOpts(doNPVExtra=True,
                xBR=xBR, yBR=yBR,
                doGenMetVsPFMet=True)
        signalHists = sms.makeSMSTemplates(curBox, signalFilename, 
                uncertainties=uncerts,
                opts=smsOpts, boostCuts=True)
        #update with correct names
        for x,h in signalHists.items():
            h.SetName(h.GetName().replace(curBox+'_'+args.model,modelName))
            signalHists[h.GetName()] = signalHists.pop(x)
        signalHist = signalHists['Signal']

        #make combined unrolled histograms for background
        makeRazorBinEvidencePlots(curBox, samples=samples, 
                inDir=BACKGROUND_DIR, outDir=outDir, unblind=True,
                signalHist=signalHist, unrollBins=unrollBins, zmin=1e-3)
        #draw signal and background in unrolled format
        if args.model == 'T2tt':
            signalString = 'pp #rightarrow #tilde{t}#tilde{t}, #mu = 1.0, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1}'
        elif args.model == 'T1bbbb':
            signalString = 'pp #rightarrow #tilde{g}#tilde{g}, #mu = 1.0, #tilde{g} #rightarrow b#bar{b}#tilde{#chi}^{0}_{1}'
        else:
            print "Model %s is not yet implemented in EvidencePlotMacro; using default signal string"%args.model
            signalString = "Signal"
        makeRazorMCTotalPlusSignalPlot(curBox, samples, 
                inDir=BACKGROUND_DIR, signalHist = signalHist, 
                outDir=outDir, unrollBins=unrollBins, signalString=signalString, 
                modelName=modelName, unblind=True) 
                

