from optparse import OptionParser
import ROOT as rt
import sys
import glob
from math import *
import os

#local imports 
from limits.SMSConfig import sms_models
from GChiPairs import gchipairs
from GetCombine import writeXsecTree
from macro.plotting import draw2DHist

def interpolateXsec(xsecs, mg):
    """
    Performs linear interpolation of the cross section
    between the two points closest to mg.
    Inputs: dictionary {mass:xsection, ...}, target mass (int)
    Returns: the interpolated cross section as well as
    the mass points used for the linear interpolation.
    """
    mg_low = -1
    mg_high = 1000000
    for m in xsecs:
        if m <= mg and m > mg_low:
            mg_low = m
        elif m >= mg and m < mg_high:
            mg_high = m
    interp_xsec = (xsecs[mg_low] + xsecs[mg_high])/2.0
    return interp_xsec, mg_low, mg_high

def plotSignificance(boxInput, model, sigHist, outDir):
    c = rt.TCanvas("c","c",800,600)
    if 'T2' in model:
        xtitle='m_{stop}'
    else:
        xtitle='m_{gluino}'
    commentstr='{} significance ({})'.format(
            model, boxInput.replace('_',' + '))
    printstr='_'.join(['RazorSignificance',model,boxInput])
    draw2DHist(c, sigHist, xtitle=xtitle, ytitle='', ztitle="Significance", printstr=printstr,
            logx=False, logy=False, logz=False, lumistr="36 fb^{-1}", 
            commentstr=commentstr, dotext=True, palette='FF', printdir=outDir, drawCMSPreliminary=False, zmin=-3, zmax=3, textSize=0.8, commentX=0.25)

if __name__ == '__main__':

    
    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Input/Output directory to store output")
    parser.add_option('--signif',dest="doSignificance",default=False,action='store_true',
                  help="for significance instead of limit")
    parser.add_option('--toys',dest="doHybridNew",default=False,action='store_true',
                  help="for toys instead of asymptotic")
    parser.add_option('--combine-name', dest='combineName',
            help='Name prefix used for combine files')
    parser.add_option('--in-dir', dest='inDir', 
            help='Input directory (if different from output directory)')

    (options,args) = parser.parse_args()

    #get options
    boxInput = options.box
    model = options.model
    directory = options.outDir
    if options.inDir is not None:
        inDir = options.inDir
    else:
        inDir = directory
    doHybridNew = options.doHybridNew
    doSignificance = options.doSignificance

    #configuration
    boxes = boxInput.split('_')
    btag = ''

    haddOutputs = []
    thyXsec = {}
    thyXsecErr = {}
    refXsecFile =  "./data/gluino13TeV.txt" 
    if 'T2' in model:
        refXsecFile = './data/stop13TeV.txt'

    #get theory cross sections and errors
    for mg in range(100,3000,25):
        for line in open(refXsecFile,'r'):
            line = line.replace('\n','')
            if str(mg)==line.split(',')[0]:
                thyXsec[mg] = float(line.split(',')[1]) #pb
                thyXsecErr[mg] = 0.01*float(line.split(',')[2])

    #significance histogram
    sigHist = None
    if doSignificance:
        sms = sms_models[model]
        sigHist = rt.TH2F('sigHist','sigHist',50,sms.mgMin,sms.mgMax,
               50,sms.mchiMin,sms.mchiMax)

    #get combine results
    submodels = sms_models[model].submodels
    if submodels is None:
        pairs = gchipairs(model)
    else:
        pairs = []
        for submodel in submodels:
            pairs += gchipairs(submodel)
    for mg, mchi in pairs:
        print "Looking for",mg,mchi
        try:
            refXsec = 1.e3*thyXsec[mg]
        except KeyError:
            refXsec, interpLow, interpHigh = interpolateXsec(thyXsec, mg)
            print ("Warning: couldn't find cross section for mass {} GeV."
                    "Interpolating between {} and {}.".format(
                        mg, interpLow, interpHigh))
        modelName = 'SMS-'+('_'.join([model, str(mg), str(mchi)]))
        
        #open file if present
        if options.combineName is not None:
            combineName = options.combineName+'_'+modelName
        else:
            combineName = 'MADD_'+options.box+'_'+modelName
        if doSignificance and doHybridNew:
            prefix = 'higgsCombineSignif'
            suffix = 'HybridNew'
        elif doHybridNew: 
            prefix = 'higgsCombineToys'
            suffix = 'HybridNew'
        elif doSignificance:
            prefix = 'higgsCombine'
            suffix = 'ProfileLikelihood'
        else:
            prefix = 'higgsCombine'
            suffix = 'Asymptotic'
        filename = inDir+'/'+prefix+combineName+'.'+suffix+'.mH120.root'
        if not glob.glob(filename): 
            print "Didn't find file",filename
            continue
        tFile = rt.TFile.Open(filename)
        print "File:",filename;

        #check file and tree within
        try:
            if tFile.InheritsFrom("TFile") is False:
                print "Isn't a TFile"
                continue
        except:
            print "File not found!"
            continue
        limit = tFile.Get("limit")
        try:
            if limit.InheritsFrom("TTree") is False: 
                tFile.cd()
                tFile.Close()
                print "'limit' doesn't appear to be a TTree"
                continue
        except:
            tFile.cd()
            tFile.Close()
            print "Error getting tree from file"
            continue
        if doSignificance and limit.GetEntries() < 1: 
            tFile.cd()
            tFile.Close()
            print "Not enough entries in tree!"
            continue
        if (not doSignificance) and limit.GetEntries() < 6: 
            tFile.cd()
            tFile.Close()
            print "Not enough entries in tree!"
            continue

        #get limits out
        limit.Draw('>>elist','','entrylist')
        elist = rt.gDirectory.Get('elist')
        entry = elist.Next()
        limit.GetEntry(entry)
        limits = []
        while True:
            if entry == -1: break
            limit.GetEntry(entry)
            if doSignificance:
                limits.append(max(0.0,limit.limit))
            else:
                limits.append(refXsec*(1.e-3)*limit.limit)
                print limit.limit
            entry = elist.Next()
        tFile.cd()
        tFile.Close()
            
        limits.reverse()
        #print limits
        
        if doSignificance:
            sigHist.SetBinContent(sigHist.FindBin(mg,mchi), limits[0])
        else:
            try:
                haddOutput = writeXsecTree(boxInput, model, directory, mg, mchi, [limits[0]],[limits[1]],[limits[2]],[limits[3]],[limits[4]],[limits[5]])
            except ReferenceError: # happens if a file is corrupt
                continue
            haddOutputs.append(haddOutput)

    if doSignificance:
        plotSignificance(boxInput, model, sigHist, directory)
    elif doHybridNew:
        os.system("hadd -f -k %s/xsecUL_HybridNew_%s.root %s"%(directory,boxInput," ".join(haddOutputs)))
    else:
        os.system("hadd -f -k %s/xsecUL_Asymptotic_%s.root %s"%(directory,boxInput," ".join(haddOutputs)))
