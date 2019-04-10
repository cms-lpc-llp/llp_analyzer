from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import glob
from RunCombine import exec_me
import datetime

if __name__ == '__main__':
    
    boxes = ['MultiJet']
    fits = ['Sideband']
    weights = ['unweighted']
    configs = ['config/run2_20151108_Preapproval_2b3b.config']

    dryRun = False
    
    lumi = 17000
    btag = '0-3btag'
    
    dateString = str(datetime.date.today()).replace("-","_")
    dateString = '2015_11_17'
    outDir = "fits_%s"%dateString
    exec_me('mkdir -p %s'%outDir,dryRun)
    
    for box in boxes:        
        for cfg in configs:                
            backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.3f_%s_%s.root'%(lumi/1000.,btag,box)
            #exec_me('python python/DustinTuple2RooDataSet.py -c %s -b %s -d Datasets/ Backgrounds/*.root -l %f -w'%(cfg,box,lumi),dryRun)
            
            if 'unweighted' in weights:
                backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_unweighted_lumi-%.3f_%s_%s.root'%(lumi/1000.,btag,box)
                #exec_me('python python/RooDataSet2UnweightedDataset.py -c %s -b %s -d Datasets/ Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.3f_%s_%s.root'%(cfg,box,lumi/1000.,btag,box),dryRun)
                
            for weight in weights:                
                outDir = "fits_%s/%s_%iifb/"%(dateString,weight,lumi/1000.)
                exec_me('mkdir -p %s'%outDir,dryRun)
                outDir = "fits_%s/%s_%iifb/%s"%(dateString,weight,lumi/1000.,box)
                exec_me('mkdir -p %s'%outDir,dryRun)
                backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_%s_lumi-%.3f_%s_%s.root'%(weight,lumi/1000.,btag,box)
                for fit in fits:                     
                    fitString = ''
                    if fit=='Sideband':
                        fitString = '--fit-region LowMR,LowRsq'

                    
                    outDir = "fits_%s/%s_%iifb/%s/%s/"%(dateString,weight,lumi/1000.,box,fit)
                    exec_me('mkdir -p %s'%outDir,dryRun)
                    #exec_me('python python/BinnedFit.py -c %s -d %s -b %s -l %f %s %s' %(cfg,outDir,box,lumi,backgroundDsName,fitString),dryRun)
                    #exec_me('python python/RunToys.py -c %s -d %s -b %s -l %f -i %s/BinnedFitResults_%s.root -t 10000' %(cfg,outDir,box,lumi,outDir,box),dryRun)
                    #exec_me('python python/PlotFit.py -c %s -d %s -b %s -l %f -i %s/BinnedFitResults_%s.root -t %s/toys_Bayes_%s.root %s' %(cfg,outDir,box,lumi,outDir,box,outDir,box,fitString),dryRun)
                    #exec_me('python python/RunToys.py -c %s -d %s -b %s -l %f -i %s/BinnedFitResults_%s.root -t 10000 --no-sys' %(cfg,outDir,box,lumi,outDir,box),dryRun)
                    exec_me('python python/RunToys.py -c %s -d %s -b %s -l %f -i %s/BinnedFitResults_%s.root -t 1000 --freq -s 0 %s' %(cfg,outDir,box,lumi,outDir,box,fitString),dryRun)
                    #exec_me('python python/RunToys.py -c %s -d %s -b %s -l %f -i %s/BinnedFitResults_%s.root -t 1000 --freq -s 1 %s' %(cfg,outDir,box,lumi,outDir,box,fitString),dryRun)
                    #exec_me('python python/RunToys.py -c %s -d %s -b %s -l %f -i %s/BinnedFitResults_%s.root -t 1000 --freq -s 2 %s' %(cfg,outDir,box,lumi,outDir,box,fitString),dryRun)
                    #exec_me('python python/PlotGOF.py -c %s -d %s -b %s -l %f -t %s/toys_Freq_s0_%s.root,%s/toys_Freq_s1_%s.root,%s/toys_Freq_s2_%s.root %s' %(cfg,outDir,box,lumi,outDir,box,outDir,box,outDir,box,fitString),dryRun)
                    exec_me('python python/PlotGOF.py -c %s -d %s -b %s -l %f -t %s/toys_Freq_s0_%s.root %s' %(cfg,outDir,box,lumi,outDir,box,fitString),dryRun)
                    
                    # for getting fit files without toys:
                    #exec_me('python python/PlotFit.py -c %s -d %s -b %s -l %f -i %s/BinnedFitResults_%s.root %s' %(cfg,outDir,box,lumi,outDir,box,fitString),dryRun)

                if len(fits)==2:
                    outDir1 = "fits_%s/%s_%iifb/%s/"%(dateString,box,lumi/1000.,fits[0])
                    outDir2 = "fits_%s/%s_%iifb/%s/"%(dateString,box,lumi/1000.,fits[1])
                    outDir = "fits_%s/%s_%iifb/"%(dateString,box,lumi/1000.)

                    # for comparing fits
                    #exec_me('python python/CompareFits.py -b %s -c %s -l %f -1 %s/BinnedFitResults_%s.root -2 %s/BinnedFitResults_%s.root -d %s'%(box,cfg,lumi,outDir1,box,outDir2,box,outDir),dryRun)
            
            

