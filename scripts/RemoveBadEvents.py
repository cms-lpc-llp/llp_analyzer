### Open ntuple and remove a specific list of events from it

import sys
import os
import argparse
import ROOT as rt

def MakeEventsDict(filename):
    """file should have lines of the form run:lumi:event"""
    d = {}
    with open(filename) as f:
        for line in f:
            #parse line
            splitLine = line.replace('\n','').split(':')
            run = long(splitLine[0])
            lumi = long(splitLine[1])
            event = long(splitLine[2])
            d[(run,lumi,event)] = True
    return d

def RemoveBadEvents(ntupleFile, badEventsDict, suffix):
    #sanity check
    if '.root' not in ntupleFile:
        print "Error: please provide a ROOT file."
        return

    #get output file name
    outputName = ''.join(ntupleFile.split('.')[:-1])+suffix+'.root'
    print "Output file:",outputName

    #open output and input files
    outputFile = rt.TFile(outputName, 'RECREATE')
    inputFile = rt.TFile(ntupleFile, 'READ')
    assert inputFile

    #find tree and filter it
    previous = None
    for key in inputFile.GetListOfKeys():

        #avoid unwanted TTree cycles by only taking the first key with a given name
        if previous is not None and previous.GetName() == key.GetName(): 
            print "Skipping unwanted cycle."
            continue
        previous = key

        #check if TTree
        className = key.GetClassName()
        print "Getting key from file.  Class type:",className
        obj = key.ReadObj()
        if className != "TTree":
            print "Not a TTree...writing directly to output"
            obj.Write()
            continue
        print "Processing tree",obj.GetName()
        print "Events in the ntuple:",obj.GetEntries()

        #create new tree for output
        outputFile.cd()
        outputTree = obj.CloneTree(0)

        #skim the tree
        eventsFailed = 0
        for i,e in enumerate(obj):
            if i % 1000000 == 0: print "Processing event ",i

            run = e.run
            lumi = e.lumi
            event = e.event
            if (run,lumi,event) not in badEventsDict:
                outputTree.Fill()
            else:
                #print "Rejected event",run,lumi,event
                eventsFailed += 1
        
        #finish
        print "Fraction of events removed:",eventsFailed,"/",obj.GetEntries(),'=',eventsFailed*1.0/obj.GetEntries()
        outputTree.Write()

    outputFile.Close()
    print "Closing output file."

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("ntupleFile", help="razor ntuple file")
    parser.add_argument("maskFile", help="file listing bad run:lumi:event")
    parser.add_argument("-s", "--suffix", help="suffix for output file name", default="_Filtered")
    args = parser.parse_args()

    ntuple = args.ntupleFile
    mask = args.maskFile
    suffix = args.suffix

    badEvents = MakeEventsDict(mask)
    RemoveBadEvents(ntuple, badEvents, suffix)
