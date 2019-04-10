import sys
import os
import argparse
import ROOT as rt

from RunCombine import exec_me

def makeFileLists(inDir, smsName, OneDScan=False):
    """
    inDir: directory to search 
    smsName: name of signal model
    OneDScan: parse only gluino mass, not LSP mass, from filename
    
    Returns: dictionary in which keys are (mGluino, mLSP) pairs
        and values are lists of ntuple files for the corresponding mass point
    """
    inFiles = os.listdir(inDir)

    #build dict of files associated with the different signal mass points
    fileLists = {}
    for f in inFiles:

        #skip files not corresponding to selected smsName
        if smsName not in f:
            continue

        #parse filename to get gluino and LSP masses
        if '.root' not in f: 
            print "Skipping non-ROOT file/directory",f
            continue
        splitF = f.replace('.root','').split('_')

        #check sanity
        if len(splitF) < 2:
            print "Unexpected file",f,": skipping"
            continue

        if not OneDScan:
            try:
                int(splitF[-2])
                mGluino = splitF[-2]
            except ValueError:
                print "Cannot parse gluino mass from",f,": skipping"
                continue

            try:
                int(splitF[-1])
                mLSP = splitF[-1]
            except ValueError:
                print "Cannot parse LSP mass from",f,": skipping"
                continue

            pair = (mGluino, mLSP)

            #add to dictionary if not present
            if pair not in fileLists:
                fileLists[pair] = []

            #add this file to appropriate list
            fileLists[pair].append(f)
            print "Adding",f,"to list of files for model point",pair

        else:
            try:
                int(splitF[-1])
                mGluino = splitF[-1]
            except ValueError:
                print "Cannot parse gluino mass from",f,": skipping"
                continue

            if mGluino not in fileLists:
                fileLists[mGluino] = []

            #add this file to appropriate list
            fileLists[mGluino].append(f)
            print "Adding",f,"to list of files for model point",mGluino

    return fileLists

def haddFastsimFiles(fileLists, smsName, inDir, outDir, OneDScan=False, dryRun=True):
    """
    fileLists: dictionary in which keys are (mGluino, mLSP) pairs
        and values are lists of ntuple files for the corresponding mass point
    smsName: name of signal model
    outDir: directory where output files should go

    Runs hadd on the files for each signal mass point and puts the output in
        the specified directory
    """
    for pair in fileLists:
        print "Signal:",smsName,pair
        if not OneDScan:
            if (os.path.isfile(outDir+'/'+smsName+'_'+pair[0]+'_'+pair[1]+'.root')):
                print "output file exists already. skip"
            else:
                exec_me('hadd -f '+outDir+'/'+smsName+'_'+pair[0]+'_'+pair[1]+'.root '
                        +' '.join([inDir+'/'+f+' ' for f in fileLists[pair]]), dryRun)
        else:
            if (os.path.isfile(outDir+'/'+smsName+'_'+pair+'.root')):
                print "output file exists already. skip"
            else:
                exec_me('hadd -f '+outDir+'/'+smsName+'_'+pair+'.root '
                        +' '.join([inDir+'/'+f+' ' for f in fileLists[pair]]), dryRun)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("inDir", help="input path")
    parser.add_argument("-s", "--smsName", required=True, help="SMS name, e.g. T1bbbb")
    parser.add_argument("--dryRun", action="store_true")
    parser.add_argument("--OneDScan", action="store_true")
    args = parser.parse_args()

    inDir=args.inDir
    smsName=args.smsName
    dryRun=args.dryRun
    OneDScan = args.OneDScan

    fileLists = makeFileLists(inDir, smsName, OneDScan)

    #output directory
    outDir = inDir+'/combined'
    exec_me('mkdir '+outDir, dryRun)

    haddFastsimFiles(fileLists, smsName, inDir, outDir, OneDScan, dryRun=dryRun)
