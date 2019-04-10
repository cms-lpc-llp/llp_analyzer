### This script combines the toys from each b-tag fit together into a single tree.

import os.path
import argparse
from array import array
import ROOT as rt

from macro.razorAnalysis import razorFitFiles, signalConfig
from framework import Config
from RunToys import getTree

def getConfigBinning(configFile, box):
    """Returns (x, y, z, nBins), where x, y, and z are the 
        binning along each axis as specified by the config
        and nBins is the total number of bins"""
    cfg = Config.Config(args.config)
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    return x, y, z, nBins

def setTreeStructBin(source, target, sourceBin, targetBin):
    toy = getattr(source, 'b{}'.format(sourceBin))
    setattr(target, 'b{}'.format(targetBin), toy)

def setCombinedTreeStruct(s1, trees, x, y, z):
    """Inserts the bin contents of the individual trees
        into the combined tree in the correct order"""
    # For 3b case the results are stored in the 2b tree
    has3b = False
    if len(z) > 4:
        has3b = True
    # This keeps track of each tree's bin numbers
    curTreeBins = [ -1 for _ in range(len(trees)) ] 
    iBinX = -1
    for ix in range(1, len(x)):
        for iy in range(1, len(y)):
            for nbtags, tree in enumerate(trees):
                iBinX += 1
                curTreeBins[nbtags] += 1
                setTreeStructBin(tree, s1, curTreeBins[nbtags], iBinX)
                if has3b and nbtags == 2:
                    iBinX += 1
                    curTreeBins[nbtags] += 1
                    setTreeStructBin(tree, s1, 
                            curTreeBins[nbtags], iBinX)


if __name__ == '__main__':
    rt.gROOT.SetBatch()
    parser = argparse.ArgumentParser()
    parser.add_argument('box', help='box name')
    parser.add_argument('--tag', help='analysis tag',
            default='Razor2016_MoriondRereco')
    parser.add_argument('--config', help='config file path',
            default=signalConfig)
    args = parser.parse_args()

    # Important: this config has to have the same binning
    # as the configs used for the individual b-tag fits.
    x, y, z, nBins = getConfigBinning(args.config, args.box)

    btagsMax = 2 # 3b fit is included with 2b
    trees = []
    for btags in range(btagsMax+1):
        subBox = "{}_{}b".format(args.box, btags)
        fitFile = razorFitFiles[args.tag][subBox]
        print "Opening file",fitFile
        toysFile = rt.TFile.Open(fitFile)
        trees.append(toysFile.Get("myTree"))

    # The combined tree will not include the function parameters,
    # only the bin yields for each toy
    combinedTree = rt.TTree("myTree", "myTree")
    treeStruct = getTree(combinedTree, paramNames=[], nBins=nBins, 
            box=args.box, z=z)

    print "Combining toys for {} fits".format(args.box)
    numEntries = trees[0].GetEntries()
    for entry in range(numEntries):
        if not entry % 1000: print "Entry ",entry
        for tree in trees:
            tree.GetEntry(entry)
        # This function does all the work
        setCombinedTreeStruct(treeStruct, trees, x, y, z)
        combinedTree.Fill()
    
    outFileName = os.path.basename(razorFitFiles[args.tag][args.box]
            ).replace('.root','_combinedBtags.root')
    outFile = rt.TFile.Open(outFileName, "recreate") 
    combinedTree.Write()
    outFile.Close()
