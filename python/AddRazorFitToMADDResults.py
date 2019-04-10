import sys
import argparse
import ROOT as rt

import macro.macro as macro
from macro.razorAnalysis import razorBinning, razorFitFiles
from macro.razorMacros import import2DRazorFitHistograms

if __name__ == '__main__':
    rt.gROOT.SetBatch()
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", default="Razor2016_MoriondRereco",
            help="Analysis tag, e.g. Razor2015")
    parser.add_argument('--box', help="choose a box", required=True)
    parser.add_argument('--force', help="override existing fit histograms",
            action='store_true')
    # Specifying a different input directory of fit file
    parser.add_argument('--no-sys', help='no systematics', 
            action='store_true', dest='noSys')
    parser.add_argument('--in-dir', help='input directory',
            dest='inDir')
    parser.add_argument('--fit-file', help='file with fit toys',
            dest='fitFile')
    args = parser.parse_args()
    
    c = rt.TCanvas("c","c",400,300) # dummy canvas
    bins = razorBinning[args.box]
    btagsMax = 3
    if args.box in ['DiJet','LeptonJet']:
        btagsMax = 2

    # Input directory
    baseDir = "Plots/%s"%args.tag
    if args.inDir is not None:
        baseDir = args.inDir

    # Fit file
    fitFile = razorFitFiles[args.tag][args.box]
    if args.fitFile is not None:
        fitFile = args.fitFile

    for btags in range(btagsMax+1):
        inDir = "%s/%s%dB"%(baseDir, args.box, btags)
        if args.noSys:
            inDir += 'NoSys'
        inFile = "razorHistograms%s%dB.root"%(args.box, btags)
        if args.noSys:
            inFile = inFile.replace(".root","NoSys.root")
        hists = macro.importHists(inDir+'/'+inFile)

        import2DRazorFitHistograms(hists, bins, fitFile, c, 
                "Data", btags, btagsMax, noStat=True, force=args.force)
        macro.exportHists(hists, outFileName=inFile, outDir=inDir)
