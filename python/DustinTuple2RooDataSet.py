from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
import sys
from array import *
import csv

from macro.razorAnalysis import Analysis

def initializeWorkspace(w,cfg,box):
    variables = cfg.getVariablesRange(box,"variables",w)
    
    w.factory('W[1.]')
    w.set('variables').add(w.var('W'))
    return w

def convertTree2Dataset(tree, cfg, box, workspace, useWeight,
        globalScaleFactor, treeName='RMRTree',isData=False, 
        tag="Razor2016_MoriondRereco"):
    """This defines the format of the RooDataSet"""
    
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    args = workspace.set("variables")
    data = rt.RooDataSet(treeName,'Selected R and MR',args)
    
    btagCutoff = 3
    if box in ["MuEle", "MuMu", "EleEle"]:
        btagCutoff = 1
    elif box in ["DiJet","LeptonJet","EleJet","MuJet", "DiJet_2b", "LeptonJet_2b"]:
        btagCutoff = 2

    # Use the Analysis object to store cuts
    # to keep synchronized with the MADD analysis
    analysis = Analysis(box, tag=tag)
    if isData:
        cuts = analysis.cutsData
    else:
        cuts = analysis.cutsMC
    
    print "Cuts:",cuts
    tree.Draw('>>elist',cuts,'entrylist')
    elist = rt.gDirectory.Get('elist')
    
    entry = -1;
    while True:
        entry = elist.Next()
        if entry == -1: break
        tree.GetEntry(entry)

        #set the RooArgSet and save
        a = rt.RooArgSet(args)
        a.setRealValue('MR',tree.MR)
        a.setRealValue('Rsq',tree.Rsq)
        a.setRealValue('nBtag',min(tree.nBTaggedJets,btagCutoff))
        if useWeight:
            a.setRealValue('W',tree.weight*globalScaleFactor)
        else:
            a.setRealValue('W',1.0)
        data.add(a)
       
    numEntries = data.numEntries()
    
    wdata = rt.RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"MR>=0.","W")
    numEntriesByBtag = []
    sumEntriesByBtag = []
    for i in z[:-1]:
        numEntriesByBtag.append(wdata.reduce('nBtag==%i'%i).numEntries())
        sumEntriesByBtag.append(wdata.reduce('nBtag==%i'%i).sumEntries())
        
    print "Number of Entries [ %s ] ="%(box), numEntriesByBtag
    print "Sum of Weights    [ %s ] ="%(box), sumEntriesByBtag

    return wdata

    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store datasets")
    parser.add_option('-w','--weight',dest="useWeight",default=False,action='store_true',
                  help="use weight")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('--lumi-in',dest="lumi_in", default=1.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('--data',dest="isData", default=False,action='store_true',
                  help="flag to use trigger decision and MET flags")


    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box =  options.box
    lumi = options.lumi
    lumi_in = options.lumi_in
    useWeight = options.useWeight
    
    print 'Input files are %s' % ', '.join(args)
    
    w = rt.RooWorkspace("w"+box)

    variables = initializeWorkspace(w,cfg,box)    
    
    ds = []
        
    btagMin =  w.var('nBtag').getMin()
    btagMax =  w.var('nBtag').getMax()

    z = array('d', cfg.getBinning(box)[2]) # nBtag binning

    for i, f in enumerate(args):
        if f.lower().endswith('.root'):
            rootFile = rt.TFile.Open(f)
            tree = rootFile.Get('RazorInclusive')
            if f.lower().find('sms')==-1:
                ds.append(convertTree2Dataset(tree, cfg, box, w, useWeight, lumi/lumi_in,  'RMRTree_%i'%i, options.isData))
                
            else:
                modelString = f.split('/')[-1].split('.root')[0].split('_')[0]
                model = modelString.split('-')[-1]
                massPoint = '_'.join(f.split('/')[-1].split('.root')[0].split('_')[1:3])
                               
                thyXsec = -1
                thyXsecErr = -1
                mGluino = -1
                mStop = -1
                if "T1" in model or "T5" in model:
                    mGluino = massPoint.split("_")[0]
                if "T2" in model:
                    mStop = massPoint.split("_")[0]
    
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

                if isinstance( rootFile.Get('NEvents'), rt.TH1 ):
                    nEvents = rootFile.Get('NEvents').Integral()
                    globalScaleFactor = thyXsec*lumi/lumi_in/nEvents # FastSim samples
                else:
                    globalScaleFactor = lumi/lumi_in # FullSim samples
                
                ds.append(convertTree2Dataset(tree, cfg, box, w, useWeight, globalScaleFactor, 'signal'))
                
    wdata = ds[0].Clone('RMRTree')
    for ids in range(1,len(ds)):
        wdata.append(ds[ids])
            
    rootTools.Utils.importToWS(w,wdata)
    
    inFiles = [f for f in args if f.lower().endswith('.root')]
    
    args = w.set("variables")
    
    #we cut away events outside our MR window
    mRmin = args['MR'].getMin()
    mRmax = args['MR'].getMax()

    #we cut away events outside our Rsq window
    rsqMin = args['Rsq'].getMin()
    rsqMax = args['Rsq'].getMax()

    btagMin =  args['nBtag'].getMin()
    btagMax =  args['nBtag'].getMax()
    
    if len(inFiles)==1:
        if btagMax>btagMin+1:
            outFile = inFiles[0].split('/')[-1].replace('.root','_lumi-%.3f_%i-%ibtag_%s.root'%(lumi/1000.,btagMin,btagMax-1,box))
            outFile = outFile.replace('_1pb','')
            if not useWeight:
                outFile = outFile.replace("weighted","unweighted")
        else:
            outFile = inFiles[0].split('/')[-1].replace('.root','_lumi-%.3f_%ibtag_%s.root'%(lumi/1000.,btagMin,box))
            outFile = outFile.replace('_1pb','')
            if not useWeight:
                outFile = outFile.replace("weighted","unweighted")
    else:
        if btagMax>btagMin+1:
            if useWeight:
                outFile = 'RazorInclusive_SMCocktail_weighted_lumi-%.3f_%i-%ibtag_%s.root'%(lumi/1000.,btagMin,btagMax-1,box)
            else:
                outFile = 'RazorInclusive_SMCocktail_unweighted_lumi-%.3f_%ibtag_%s.root'%(lumi/1000.,btagMin,box)
        else:
            if useWeight:
                outFile = 'RazorInclusive_SMCocktail_weighted_lumi-%.3f_%ibtag_%s.root'%(lumi/1000.,btagMin,box)
            else:
                outFile = 'RazorInclusive_SMCocktail_unweighted_lumi-%.3f_%ibtag_%s.root'%(lumi/1000.,btagMin,box)
        

    numEntriesByBtag = []
    sumEntriesByBtag = []
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    for i in z[:-1]:
        numEntriesByBtag.append(wdata.reduce('nBtag==%i'%i).numEntries())
        sumEntriesByBtag.append(wdata.reduce('nBtag==%i'%i).sumEntries())
        
    print "Output File: %s"%(options.outDir+"/"+outFile)
    print "Number of Entries Total [ %s ] ="%(box), numEntriesByBtag
    print "Sum of Weights Total    [ %s ] ="%(box), sumEntriesByBtag
    
    outFile = rt.TFile.Open(options.outDir+"/"+outFile,'recreate')
    outFile.cd()
    w.Write()
    outFile.Close()
    
