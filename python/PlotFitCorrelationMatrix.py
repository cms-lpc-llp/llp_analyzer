import sys
import ROOT as rt
import argparse
import array

from PlotFit import getCorrelationMatrix
from framework import Config
from macro.plotting import draw2DHist

config = "config/run2_20151108_Preapproval_2b3b_data.config"

if __name__ == '__main__':
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("toyFile", help="Input file with toys from fit result")
    parser.add_argument("xMin", help="Bin number of first MR bin", type=int)
    parser.add_argument("xMax", help="Bin number of last MR bin", type=int)
    parser.add_argument("yMin", help="Bin number of first Rsq bin", type=int)
    parser.add_argument("yMax", help="Bin number of last Rsq bin", type=int)
    parser.add_argument("--btags", help="choose number of btags", type=int)
    parser.add_argument('--box', help="choose a box")
    args = parser.parse_args()

    if args.box is None:
        print "Please provide the name of the box using --box"
        sys.exit()
    if args.btags is None:
        print "Please provide the number of btags using --btags"
        sys.exit()

    zMin = args.btags
    zMax = args.btags+1

    print "Creating fit correlation matrix"

    cfg = Config.Config(config)
    x = array.array('d', cfg.getBinning(args.box)[0]) # MR binning
    y = array.array('d', cfg.getBinning(args.box)[1]) # Rsq binning
    z = array.array('d', cfg.getBinning(args.box)[2]) # nBtag binning

    #load fit information, including toys
    toyFile = rt.TFile.Open(args.toyFile)
    assert toyFile
    print "Opened file",args.toyFile,"to get Bayesian toy results"
    toyTree = toyFile.Get("myTree")
    assert toyTree
    print "Got tree myTree"

    corrMatrix = getCorrelationMatrix(toyTree,'zyx',args.xMin,args.xMax,args.yMin,args.yMax,zMin,zMax,x,y,z)

    c = rt.TCanvas("c", "c", 800, 600)
    printstr = "fitCorrelation"+args.box+str(args.btags)+'B_'+('_'.join([str(args.xMin),str(args.xMax),str(args.yMin),str(args.yMax)]))

    draw2DHist(c, corrMatrix, xtitle="Bin", ytitle="Bin", printstr=printstr, logx=False, logy=False, logz=False, zmin=-1.0, zmax=1.0, numDigits=2, drawCMSPreliminary=False, savepdf=False, saveroot=False, savec=False, palette="FF")
