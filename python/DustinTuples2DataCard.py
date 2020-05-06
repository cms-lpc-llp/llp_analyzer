import os
import argparse
import ROOT as rt
import sys
from array import *

#local imports
import macro.macro as macro

lumiUncertainty = 0.026

def getTheoryCrossSectionAndError(mGluino=-1, mStop=-1):
    thyXsec = -1
    thyXsecErr = -1

    if mGluino!=-1:
        for line in open('data/gluino13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(mGluino))==line.split(',')[0]:
                thyXsec = float(line.split(',')[1]) #pb
                thyXsecErr = 0.01*float(line.split(',')[2])
    if mStop!=-1:
        for line in open('data/stop13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(mStop))==line.split(',')[0]:
                thyXsec = float(line.split(',')[1]) #pb
                thyXsecErr = 0.01*float(line.split(',')[2]) 

    return thyXsec,thyXsecErr

def uncorrelate(hists, sysName, suppressLevel=None):
    """Replaces each histogram whose name contains 'sysName' with many copies that represent uncorrelated bin-by-bin systematics.
    suppressLevel: if provided, new histograms will only be created for bins that differ from nominal by a fractional amount greater than suppressLevel."""
    #get all histograms that match the input string
    toUncorrelate = [name for name in hists if sysName in name]
    print "Treating the following distributions as uncorrelated for",sysName,": "
    for name in toUncorrelate: print name
    
    #get names of individual systematics
    systNames = []
    for name in toUncorrelate:
        systName = name.replace("Up","").replace("Down","")
        if systName not in systNames:
            systNames.append(systName)

    for name in systNames:
        print("Uncorrelating "+name)
        #get histogram with central values
        centerName = name.split("_")[:-1]
        centerName = '_'.join(centerName)
        #get up and down variants
        upName = name+'Up'
        downName = name+'Down'
        uncName = name.split("_")[-1]
        print("Central values taken from "+centerName)
        #for each bin create a new histogram in which that bin is up/down and the rest are centered
        for b in range(1,hists[centerName].GetNbinsX()+1):
            newHistUpName = centerName+"_"+uncName+str(b)+"Up"
            newHistDownName = centerName+"_"+uncName+str(b)+"Down"

            #check level of agreement with the nominal
            if suppressLevel is not None:
                #get percent difference from nominal
                if hists[centerName].GetBinContent(b) > 0:
                    percDifferenceUp = abs(hists[upName].GetBinContent(b)-hists[centerName].GetBinContent(b))/hists[centerName].GetBinContent(b)
                    percDifferenceDown = abs(hists[downName].GetBinContent(b)-hists[centerName].GetBinContent(b))/hists[centerName].GetBinContent(b)
                    percDifference = max(percDifferenceUp, percDifferenceDown)
                    if percDifference <= suppressLevel: 
                        #print "Suppressing nuisance in bin",b,"(agrees at",percDifference,"level)"
                        continue
                elif hists[upName].GetBinContent(b) == hists[centerName].GetBinContent(b) and hists[downName].GetBinContent(b) == hists[centerName].GetBinContent(b): 
                        #print "Suppressing nuisance in bin",b,"because there is no change from the nominal"
                        continue

            #new up histogram
            hists[newHistUpName] = hists[centerName].Clone(newHistUpName)
            hists[newHistUpName].SetDirectory(0)
            hists[newHistUpName].SetBinContent(b, hists[upName].GetBinContent(b)) #new hist has the unperturbed value in every bin except one
            hists[newHistUpName].SetBinError(b, hists[upName].GetBinError(b))

            #new down histogram
            hists[newHistDownName] = hists[centerName].Clone(newHistDownName)
            hists[newHistDownName].SetDirectory(0)
            hists[newHistDownName].SetBinContent(b, hists[downName].GetBinContent(b)) #new hist has the unperturbed value in every bin except one
            hists[newHistDownName].SetBinError(b, hists[downName].GetBinError(b))

        #remove the original histogram
        del hists[upName]
        del hists[downName]

def makeNewHistogramForUncorrelateSFs(name, centerName, systName, number, hists):
    """Set up a new histogram, for use with the uncorrelateSFs function"""
    if "Up" in name: 
        newHistName = centerName+"_"+systName+str(number)+"Up"
    elif "Down" in name:
        newHistName = centerName+"_"+systName+str(number)+"Down"
    else: 
        print("Error: shape histogram name "+name+" needs to contain 'Up' or 'Down'")
        return
    if newHistName not in hists:
        hists[newHistName] = hists[centerName].Clone(newHistName)
        hists[newHistName].SetDirectory(0)
    return newHistName

def setBinContentsForUncorrelateSFs(mrCenter, rsqCenter, refBN, sigBN, sysHist, newHist, referenceHist, oneD=False):
    """If the signal bin is inside the reference bin, perturb the signal bin"""
    #correct MR or Rsq if they lie outside the reference histogram
    if mrCenter > referenceHist.GetXaxis().GetXmax(): 
        mrCenter = referenceHist.GetXaxis().GetXmax() - 1
    if not oneD and rsqCenter > referenceHist.GetYaxis().GetXmax():
        rsqCenter = referenceHist.GetYaxis().GetXmax() - 0.01
    #if the bin matches the current reference histogram bin, update the contents
    if oneD:
        matches = (referenceHist.FindBin(mrCenter) == refBN)
    else:
        matches = (referenceHist.FindBin(mrCenter, rsqCenter) == refBN)
    if matches:
        newHist.SetBinContent(sigBN, sysHist.GetBinContent(sigBN)) #new hist has the unperturbed value in every bin except one
        newHist.SetBinError(sigBN, sysHist.GetBinError(sigBN))
        return True
    else:
        return False

def uncorrelateSFs(hists, sysName, referenceHists, cfg, box, unrollBins=None):
    """Same as uncorrelate(), but treats bins as correlated if they lie inside the same bin in the reference histogram.
    Needs a config and a box name, to get the correct bin configuration for the razor histogram"""
    #get all histograms that match the input string
    toUncorrelate_tmp = [name for name in hists if sysName in name]
    # VERY AD HOC
    # to exclude uncertainties with overlapping names
    toUncorrelate = []
    badStrs = ['MR', 'NJets', 'NBTags']
    for name in toUncorrelate_tmp:
        ok = True
        for s in badStrs:
            if s in name:
                print "Excluding {} because its name contains {}".format(
                        name, s)
                ok = False
        if ok: 
            toUncorrelate.append(name)
    print "Uncorrelate SFs:",sysName
    print("Treating the following distributions as uncorrelated: ")
    for name in toUncorrelate: print name

    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    #make histogram with razor binning
    myTH3 = rt.TH3D("razor3d"+name,"razor3d",len(x)-1,x,len(y)-1,y,len(z)-1,z)

    for name in toUncorrelate:
        print("Using reference histogram to determine bin correlations for "+name)
        #get histogram with central values
        centerName = name.split("_")[:-1]
        centerName = '_'.join(centerName)
        systName = name.split("_")[-1].replace("Up","").replace("Down","")
        print("Central values taken from "+centerName)
        #get reference histogram for scale factor binning
        referenceHist = referenceHists[centerName]
        #for each bin create a new histogram in which that bin is up/down and the rest are centered
        if referenceHist.InheritsFrom("TH2Poly"):
            for bn in range(1,referenceHist.GetNumberOfBins()+1):
                matchedAtLeastOneBin = False
                newHistName = makeNewHistogramForUncorrelateSFs(name, centerName, systName, bn, hists)
                #find all bins of signal histogram that are within this bin
                if unrollBins is None:
                    i = 0
                    for ix in range(1,len(x)):
                        for iy in range(1,len(y)):
                            for iz in range(1,len(z)):
                                #i = 1D histogram bin index
                                i += 1
                                #get MR and Rsq at center of bin in 3d histogram
                                mrCenter = myTH3.GetXaxis().GetBinCenter(ix)
                                rsqCenter = myTH3.GetYaxis().GetBinCenter(iy)
                                if setBinContentsForUncorrelateSFs(mrCenter, rsqCenter, refBN=bn,
                                        sigBN=i, sysHist=hists[name], newHist=hists[newHistName],
                                        referenceHist=referenceHist):
                                    #print " bin",i,"(",ix,iy,iz,") matches!"
                                    matchedAtLeastOneBin = True
                else:
                    #make a TH2Poly from each xy slice of the histogram, unroll each one and attach together
                    i = 0
                    for iz in range(1, len(z)):
                        #one xy slice of the histogram 
                        unrollRows = unrollBins[iz-1][0]
                        unrollCols = unrollBins[iz-1][1]
                        poly = macro.makeTH2PolyFromColumns("poly"+str(iz)+name, 'poly', unrollRows, unrollCols)
                        polyBins = poly.GetBins()
                        for sigBN in range(1, poly.GetNumberOfBins()+1):
                            i += 1
                            thisSigBin = polyBins.At(sigBN-1)
                            #get MR and Rsq at center of bin in 3d histogram
                            mrCenter = (thisSigBin.GetXMax() + thisSigBin.GetXMin())/2.0
                            rsqCenter = (thisSigBin.GetYMax() + thisSigBin.GetYMin())/2.0
                            if setBinContentsForUncorrelateSFs(mrCenter, rsqCenter, refBN=bn,
                                    sigBN=i, sysHist=hists[name], newHist=hists[newHistName],
                                    referenceHist=referenceHist):
                                #print " bin",i,"(",sigBN,iz,") matches!"
                                matchedAtLeastOneBin = True
                        poly.Delete()
                #don't save the histogram if there is no change from the nominal
                if not matchedAtLeastOneBin:
                    #print "No matching signal bins -- discarding histogram"
                    del hists[newHistName]
        else: #TH2F case
            for bx in range(1,referenceHist.GetNbinsX()+1):
                for by in range(1,referenceHist.GetNbinsY()+1):
                    b = referenceHist.GetBin(bx,by)
                    newHistName = makeNewHistogramForUncorrelateSFs(name, centerName, systName, b, hists)
                    matchedAtLeastOneBin = False
                    #find bins in hists[name] that lie inside bin b of referenceHist
                    if unrollBins is None:
                        i = 0
                        for ix in range(1,len(x)):
                            for iy in range(1,len(y)):
                                for iz in range(1,len(z)):
                                    #i = 1D histogram bin index
                                    i+= 1
                                    #get MR and Rsq at center of bin in 3d histogram
                                    mrCenter = myTH3.GetXaxis().GetBinCenter(ix)
                                    rsqCenter = myTH3.GetYaxis().GetBinCenter(iy)
                                    if setBinContentsForUncorrelateSFs(mrCenter, rsqCenter, refBN=b,
                                            sigBN=i, sysHist=hists[name], newHist=hists[newHistName],
                                            referenceHist=referenceHist):
                                        #print " bin",i,"(",ix,iy,iz,") matches!"
                                        matchedAtLeastOneBin = True
                    else:
                        #make a TH2Poly from each xy slice of the histogram, unroll each one and attach together
                        i = 0
                        for iz in range(1, len(z)):
                            #one xy slice of the histogram 
                            unrollRows = unrollBins[iz-1][0]
                            unrollCols = unrollBins[iz-1][1]
                            poly = macro.makeTH2PolyFromColumns("poly"+str(iz)+name, 'poly', unrollRows, unrollCols)
                            polyBins = poly.GetBins()
                            for sigBN in range(1, poly.GetNumberOfBins()+1):
                                i += 1
                                thisSigBin = polyBins.At(sigBN-1)
                                #get MR and Rsq at center of bin in 3d histogram
                                mrCenter = (thisSigBin.GetXMax() + thisSigBin.GetXMin())/2.0
                                rsqCenter = (thisSigBin.GetYMax() + thisSigBin.GetYMin())/2.0
                                if setBinContentsForUncorrelateSFs(mrCenter, rsqCenter, refBN=b,
                                        sigBN=i, sysHist=hists[name], newHist=hists[newHistName],
                                        referenceHist=referenceHist):
                                    #print " bin",i,"(",sigBN,iz,") matches!"
                                    matchedAtLeastOneBin = True
                            poly.Delete()
                    #don't save the histogram if there is no change from the nominal
                    if not matchedAtLeastOneBin:
                        #print "No matching signal bins -- discarding histogram"
                        del hists[newHistName]

        #remove the original histogram
        del hists[name]

def uncorrelateSFs1D(hists, sysName, referenceHists, unrollBins, 
        useRsq=True, bInclusive=False, xInclusive=False):
    """
    Same as uncorrelateSFs except that the reference histogram is assumed to
    be a 1D Rsq histogram instead of a 2D MR-Rsq histogram.
    (to decorrelate by MR value instead of Rsq, set useRsq to False)
    Set bInclusive to True to keep different b-tag bins correlated.
    Set xInclusive to True to have one shape uncertainty for each b-tag bin.
    """
    if bInclusive and xInclusive:
        print "bInclusive and xInclusive are both set to True.  Nothing to do..."
        return
    toUncorrelate = [name for name in hists if sysName in name]
    print "Uncorrelate SFs:",sysName
    print("Treating the following distributions as uncorrelated: ")
    for name in toUncorrelate: print name
    for name in toUncorrelate:
        print("Using reference histogram to determine bin correlations for "+name)
        centerName = name.split("_")[:-1]
        centerName = '_'.join(centerName)
        systName = name.split("_")[-1].replace("Up","").replace("Down","")
        print("Central values taken from "+centerName)
        referenceHist = referenceHists[centerName]
        ibin = 0
        for bx in range(1,referenceHist.GetNbinsX()+1):
            sigBNGlobal = 0
            if bInclusive:
                ibin += 1
                newHistName = makeNewHistogramForUncorrelateSFs(name, centerName, systName, ibin, hists)
                matchedAtLeastOneBin = False
            if xInclusive:
                ibin = 0
            for bz in range(len(unrollBins)):
                if not bInclusive:
                    ibin += 1
                    newHistName = makeNewHistogramForUncorrelateSFs(name, centerName, systName, ibin, hists)
                    matchedAtLeastOneBin = False
                unrollRows = unrollBins[bz][0]
                unrollCols = unrollBins[bz][1]
                poly = macro.makeTH2PolyFromColumns("poly"+str(bz)+name, 'poly', unrollRows, unrollCols)
                polyBins = poly.GetBins()
                for sigBN in range(1, poly.GetNumberOfBins()+1):
                    sigBNGlobal += 1
                    thisSigBin = polyBins.At(sigBN-1)
                    if useRsq:
                        binCenter = (thisSigBin.GetYMax() + thisSigBin.GetYMin())/2.0
                    else:
                        binCenter = (thisSigBin.GetXMax() + thisSigBin.GetXMin())/2.0
                    # Note that binCenter is an Rsq value but it is being passed
                    # into the mrCenter field of setBinContentsForUncorrelateSFs.  
                    # This is just for convenience -- that function assumes
                    # that MR values are on the x-axis and Rsq is on y (not true here)
                    if setBinContentsForUncorrelateSFs(binCenter, -1, refBN=bx,
                            sigBN=sigBNGlobal, sysHist=hists[name], newHist=hists[newHistName],
                            referenceHist=referenceHist, oneD=True):
                        matchedAtLeastOneBin = True
                poly.Delete()
                #don't save the histogram if there is no change from the nominal
                if not bInclusive and not xInclusive and not matchedAtLeastOneBin:
                    print "No matches for bin {} {}".format(bx, bz)
                    del hists[newHistName]
            if bInclusive and not matchedAtLeastOneBin:
                    print "No matches for bin {}".format(bx)
                    del hists[newHistName]
        del hists[name]

def writeDataCard_th1(box,txtfileName,hists,bkgs):
    signalName = 'Signal'
    obsRate = hists["data_obs"].Integral()
    nBkgd = len(bkgs)
    rootFileName = txtfileName.replace('.txt','.root')
    rates = [hists[signalName].Integral()]
    rates += [hists[bkg].Integral() for bkg in bkgs]
    processes = [signalName]+bkgs
    lumiErrs = [1+lumiUncertainty] 
    noLumiErrProcs = ['ttjets1l', 'ttjets2l', 'wjets', 'zinv', 'qcd']
    lumiErrs += [1+lumiUncertainty if bkg.lower() 
            not in noLumiErrProcs else 1.0 for bkg in bkgs] 
    mcErrs = {} #dictionary of uncorrelated mc bkgd lnN uncertainties

    #get list of shape uncertainties
    shapeNames = []
    for name in hists:
        if "Down" in name:
            shapeName = (name.split("_")[-1]).replace("Down","") 
            if shapeName not in shapeNames: shapeNames.append(shapeName)
    shapeErrs = { name:["1.0"] if signalName+"_"+name+"Down" in hists 
            else ["-"] for name in shapeNames }
    for name in shapeNames: 
        shapeErrs[name] += ["1.0" if bkg+"_"+name+"Down" in hists 
            else "-" for bkg in bkgs]

    #20% normalization uncertainty on rare backgrounds
    for bkg in bkgs:
        if bkg.lower() in ['ttjets', 'ttjets1l', 'ttjets2l', 
                'wjets', 'dyjets', 'zinv', 'qcd', 
                'singletop', 'other']: continue
        mcErrs[bkg] = [1.00]
        mcErrs[bkg] += [1.00 + 0.20*(bkg==bkg1) for bkg1 in bkgs]
            
    divider = "------------------------------------------------------------\n"
    datacard = "imax 1 number of channels\n" + \
              "jmax %i number of backgrounds\n"%nBkgd + \
               "kmax * number of nuisance parameters\n" + \
               divider + \
               "observation	%.3f\n"%obsRate + \
               divider + \
               "shapes * * %s $PROCESS $PROCESS_$SYSTEMATIC\n"%(
                       rootFileName) + \
               divider
               
    binString = "bin"
    processString = "process"
    processNumberString = "process"
    rateString = "rate"
    lumiString = "lumi\tlnN"
    for i in range(0,len(bkgs)+1):
        binString +="\t%s"%box
        processString += "\t%s"%processes[i]
        processNumberString += "\t%i"%i
        rateString += "\t%.3f" %rates[i]
        lumiString += "\t%.3f"%lumiErrs[i]
    binString+="\n" 
    processString+="\n"
    processNumberString+="\n"
    rateString +="\n"
    lumiString+="\n"
        
    mcErrStrings = {}
    for bkg in bkgs:
        if bkg not in mcErrs: continue
        mcErrStrings[bkg] = "%s_norm\tlnN"%(bkg)
        for i in range(0,len(bkgs)+1):                
            mcErrStrings[bkg] += "\t%.3f"%mcErrs[bkg][i]
        mcErrStrings[bkg]+="\n"
    shapeErrStrings = {name:name+"\tshape" for name in shapeNames}
    for name in shapeNames: 
        for i in range(0, len(bkgs)+1):
            shapeErrStrings[name] += "\t"+shapeErrs[name][i]
        shapeErrStrings[name]+="\n"

    datacard+=(binString+processString+processNumberString+
            rateString+divider+lumiString)
    for bkg in bkgs:
        if bkg in mcErrStrings: datacard+=mcErrStrings[bkg] 
    for name in sorted(shapeNames):
        datacard+=shapeErrStrings[name] #shape uncertainties

    #print out card
    print "\n",datacard,"\n"

    #write card
    txtfile = open(txtfileName,"w")
    txtfile.write(datacard)
    txtfile.close()
