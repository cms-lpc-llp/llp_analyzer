from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import glob
from RunCombine import exec_me

if __name__ == '__main__':


    boxes = ['TTBarSingleLepton','WSingleLepton']
    #fits = ['Sideband']
    #fits = ['Full']
    fits = ['Sideband','Full']
    configs = ['config/controlsample.config']
    
    
    lumi = 23.8

    dryRun=False
    
    for box in boxes:
        for cfg in configs:            

            exec_me('python python/CRTuple2RooDataSet.py -c %s -b %s -d ControlSampleFits/ ControlSampleFits/SingleMuonAndElectron_Run2015B-GOLDEN.root -l %f'%(cfg,box,lumi),dryRun)
            #exec_me('python python/CRTuple2RooDataSet.py -c %s -b %s -d ControlSampleFits/ ControlSampleFits/RunTwoRazorControlRegions_OneLeptonFull_Run2015B_GoodLumiDCS_NoDuplicates.root -l %f'%(cfg,box,lumi),dryRun)
            
            for fit in fits:                
                if box in ['TTBarSingleLepton']:
                    btag = '1btag'
                else:
                    btag = '0btag'
                
                lumiString = '%.4f'%(lumi/1000)
                #lumiString = lumiString.replace('.','p')

                fitString = ''
                if fit=='Sideband':
                    fitString = '--fit-region LowMR,LowRsq'
                    
                outDir = "ControlSampleFits/%s/"%(box)
                exec_me('mkdir -p %s' %(outDir),dryRun)
                outDir = "ControlSampleFits/%s/%s/"%(box,fit)
                exec_me('mkdir -p %s' %(outDir),dryRun)
                dsName = 'ControlSampleFits/SingleMuonAndElectron_Run2015B-GOLDEN_lumi-%s_%s_%s.root'%(lumiString,btag,box)
                #dsName = 'ControlSampleFits/RunTwoRazorControlRegions_OneLeptonFull_Run2015B_GoodLumiDCS_NoDuplicates_lumi-%s_%s_%s.root'%(lumiString,btag,box)
                exec_me('python python/BinnedFit.py -c %s -b %s -l %f -d %s %s %s --data' %(cfg,box,lumi,outDir,fitString,dsName),dryRun)
                #exec_me('python python/RunToys.py -c %s -b %s -i %s/BinnedFitResults_%s.root -d %s -t 3000'%(cfg,box,outDir,box,outDir),dryRun)
                #exec_me('python python/PlotFit.py -c %s -b %s -l %f -i %s/BinnedFitResults_%s.root -d %s -t %s/toys_Bayes_%s.root --data' %(cfg,box,lumi,outDir,box,outDir,outDir,box),dryRun)
                exec_me('python python/PlotFit.py -c %s -b %s -l %f -i %s/BinnedFitResults_%s.root -d %s %s --data' %(cfg,box,lumi,outDir,box,outDir,fitString),dryRun)
                      
            
            
