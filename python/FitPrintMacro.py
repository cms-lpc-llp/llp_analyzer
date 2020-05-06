### Load fit results from toy file and print out the yields bin by bin 

import sys,os
import argparse
import ROOT as rt

#local imports
import macro.macro as macro
from macro.razorMacros import *

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("inFile", help="Toy file used to load fit results")
    parser.add_argument("-c", "--config", help="Config file to use", default="config/run2_sideband.config")
    parser.add_argument("-b", "--box", help="Analysis box", default="MultiJet")
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument("-s", "--signal", help="Signal ntuple (razor format)", default=None)
    parser.add_argument("-m", "--model", help="Signal model", default="T1bbbb")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug

    #get signal histogram, if any
    signalHist = None
    threeDeeSignalHist = None
    if args.signal is not None:
        print "Opening signal file",args.signal
        signalFile = rt.TFile.Open(args.signal)
        assert signalFile
        signalHist = signalFile.Get(args.box+"_"+args.model)
        assert signalHist
        from PlotFit import get3DHistoFrom1D
        cfg = Config.Config(args.config)
        binsX = array('d', cfg.getBinning(args.box)[0])
        binsY = array('d', cfg.getBinning(args.box)[1])
        binsZ = array('d', cfg.getBinning(args.box)[2])
        threeDeeSignalHist = get3DHistoFrom1D(signalHist,binsX,binsY,binsZ,"signalHist3D") 

    #fit yield table with uncertainties
    fitHist = get3DRazorFitHistogram(args.config, args.inFile, args.box, debugLevel=debugLevel)
    makeRazor3DTable(fitHist, args.box, threeDeeSignalHist)

    #fit correlation matrix
    corrMatrix = getFitCorrelationMatrix(args.config, args.box, args.inFile, debugLevel=debugLevel)
    c = rt.TCanvas('corr','corr',800,600)
    macro.draw2DHist(c, corrMatrix, "Bin number", "Bin number", "Correlation coefficient", zmin=-1.0, zmax=1.0, printstr='fitCorrelation'+args.box, dotext=False, logx=False, logy=False, logz=False, saveroot=True)

