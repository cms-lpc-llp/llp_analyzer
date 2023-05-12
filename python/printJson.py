import json
from optparse import OptionParser
import ROOT as rt
import sys
from array import *
import os
from itertools import *
from operator import *
import pickle

def walk(top, topdown=True):
    """
    os.path.walk like function for TDirectories.
    Return 4-tuple: (dirpath, dirnames, filenames, top)
        dirpath = 'file_name.root:/some/path' # may end in a '/'?
        dirnames = ['list', 'of' 'TDirectory', 'keys']
        filenames = ['list', 'of' 'object', 'keys']
        top = this level's TDirectory
    """
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    assert isinstance(top, rt.TDirectory)
    names = [k.GetName() for k in top.GetListOfKeys()]
    dirpath = top.GetPath()
    dirnames = []
    filenames = []
    ## filter names for directories
    for k in names:
        d = top.Get(k)
        if isinstance(d, rt.TDirectory):
            dirnames.append(k)
        else:
            filenames.append(k)
    ## sort
    dirnames.sort()
    filenames.sort()
    ## yield
    if topdown:
        yield dirpath, dirnames, filenames, top
    for dn in dirnames:
        d = top.Get(dn)
        for x in walk(d, topdown):
            yield x
    if not topdown:
        yield dirpath, dirnames, filenames, top

    
def convertTree2Dict(runLumiDict,tree,lumiBranch,runBranch):

    if not (hasattr(tree,lumiBranch) and hasattr(tree,runBranch)):
        print "tree does not contain run and lumi branches, returning empty json"
        return runLumiDict
    # loop over tree to get run, lumi "flat" dictionary
    tree.Draw('>>elist','','entrylist')        
    elist = rt.gDirectory.Get('elist')    
    if tree.GetEntries()==0:return runLumiDict
    entry = -1;
    while True:
        entry = elist.Next()
        if entry == -1: break
        tree.GetEntry(entry)
        if entry%10000==0:
            print "processing entry %i"%entry
        if '%s'%(getattr(tree,runBranch)) in runLumiDict.keys():
            currentLumi = runLumiDict['%s'%(getattr(tree,runBranch))]
            if int(getattr(tree,lumiBranch)) in currentLumi:
                pass
            else:                
                currentLumi.append(int(getattr(tree,lumiBranch)))
                runLumiDict.update({'%s'%(getattr(tree,runBranch)):currentLumi})
        else:
            runLumiDict['%s'%(getattr(tree,runBranch))] = [int(getattr(tree,lumiBranch))]
        
    return runLumiDict

def fixDict(runLumiDict):
    # fix run, lumi list by grouping consecutive lumis
    for run in runLumiDict.keys():
        lumiGroups = []
        for k, g in groupby(enumerate(runLumiDict[run]), lambda (i,x):i-x):
            consecutiveLumis = map(itemgetter(1), g)
            lumiGroups.append([consecutiveLumis[0],consecutiveLumis[-1]])
        runLumiDict.update({run:lumiGroups})
        
    return runLumiDict
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-o','--output',dest="output",type="string",default="test.json",
                  help="Name of the json file to write to")
    parser.add_option('-l','--lumi-branch',dest="lumiBranch",type="string",default="lumi",
                  help="Name of lumi branch in tree")
    parser.add_option('-r','--run-branch',dest="runBranch",type="string",default="run",
                  help="Name of run branch in tree")
    
    parser.add_option('-i','--input-list',dest="inputList",type="string",default="test.list",
                  help="Name of text file containing list of ROOT files")
    (options,args) = parser.parse_args()


    rootFiles = []
    #for f in args:
    #    if f.lower().endswith('.root'):
    #        rootFile = rt.TFile.Open(f)
    #        rootFiles.append(rootFile)
   

    file1 = open(options.inputList, 'r')
    Lines = file1.readlines()
     
    for line in Lines:
	rootFile = rt.TFile.Open(line.strip())
        rootFiles.append(rootFile)



    trees = []
    # crawl root file to look for trees
    for rootFile in rootFiles:
        for dirpath, dirnames, filenames, tdirectory in walk(rootFile):
            for filename in set(filenames):
                obj = tdirectory.Get(filename)
                if isinstance(obj, rt.TTree) and obj!=None:
                    print "found tree %s in directory %s in file %s"%(obj.GetName(),tdirectory.GetName(),rootFile.GetName())
                    trees.append(obj)

    # loop over trees found
    runLumiDict = {}
    for tree in trees:
        runLumiDict = convertTree2Dict(runLumiDict,tree,options.lumiBranch,options.runBranch)
    runLumiDict = fixDict(runLumiDict)
    output = open(options.output,'w')
    json.dump(runLumiDict,output,sort_keys=True)
    output.close()
    print '\njson dumped to file %s:'%options.output
    os.system('cat %s'%options.output)
    print '\n'
            
    
