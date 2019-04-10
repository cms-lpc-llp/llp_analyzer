import sys
import os.path
import argparse
import ROOT as rt

#local imports
from framework import Config
from GChiPairs import gchipairs
from SidebandMacro import LUMI
from CheckSignalContamination import checkSignalContamination

def usage():
    print "usage: ScanSignalContamination.py modelName"

if __name__ == '__main__':
    rt.gROOT.SetBatch()

    if len(sys.argv) != 2:
        usage()
        sys.exit()
    model = sys.argv[1]
    pairs = gchipairs(model)

    contam = rt.TH2F("contam", "Signal contamination", 80, 0, 2000, 40, 0, 1000)
    for p in pairs:
        mLSP = p[1]
        if 'T2' in model:
            mGluino = -1
            mStop   = p[0]
            outFileName = 'contamination_%s_%d_%d.root'%(model,mStop,mLSP)
        else:
            mGluino = p[0]
            mStop   = -1
            outFileName = 'contamination_%s_%d_%d.root'%(model,mGluino,mLSP)
        try:
            h = checkSignalContamination("config/run2_20151229_ControlRegion.config", outDir='SMSPlots', lumi=LUMI, box="TTJetsSingleLeptonControlRegion", model=model, mGluino=mGluino, mStop=mStop, mLSP=mLSP, mergeBins=True, treeName="RazorInclusive")
        except:
            print "Problem getting histogram for",p
            continue
        #get maximum signal contamination
        maxContam = 0.0
        for bx in range(1,h.GetSize()):
            maxContam = max(maxContam, h.GetBinContent(bx)) 
        if 'T2' in model:
            contam.Fill(mStop,mLSP,maxContam)
        else:
            contam.Fill(mGluino,mLSP,maxContam)
        print p, "signal contamination", maxContam
        #output
        f = rt.TFile(outFileName,'recreate')
        h.Write()
        f.Close()
        h.Delete()

    c = rt.TCanvas("c", "c", 800, 600)
    if 'T2' in model:
        contam.GetXaxis().SetTitle("Stop mass")
    else:
        contam.GetXaxis().SetTitle("Gluino mass")
    contam.GetYaxis().SetTitle("LSP mass")
    contam.SetStats(0)
    contam.Draw("colz")
    c.Print("signalContaminationScan"+model+".root")
