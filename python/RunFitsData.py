import os
import ROOT as rt    
from RunCombine import exec_me
import datetime

if __name__ == '__main__':
    
    boxes = ['MultiJet']
    fits = ['Sideband']
    #lumi = 2185.
    lumi = 201000.
    #config = 'config/run2_20151108_Preapproval_2b3b_data.config'
    config = 'config/run2_20160610_data.config'
    
    #dataset = {'MultiJet':'RazorInclusive_HTMHT_Run2015D_GoodLumiGolden_RazorSkim_Filtered_lumi-2.185_0-3btag_MultiJet.root',
    #           'MuMultiJet':'RazorInclusive_SingleMuon_Run2015D_GoodLumiGolden_RazorSkim_Filtered_lumi-2.185_0-3btag_MuMultiJet.root',
    #           'EleMultiJet':'RazorInclusive_SingleElectron_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter_lumi-2.185_0-3btag_EleMultiJet.root',
    #           'DiJet':'RazorInclusive_HTMHT_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter_lumi-2.185_0-3btag_DiJet.root'
    #           }
        
    dataset = {'MultiJet':'FullRazorInclusive_HTMHT_2016B_PRv2_GoodLumiGoldenJun16.root',
               'DiJet':'FullRazorInclusive_HTMHT_2016B_PRv2_GoodLumiGoldenJun16.root',
               'LeptonMultiJet':'FullRazorInclusive_SingleLepton_2016B_PRv2_GoodLumiGoldenJun16_NoDuplicates.root',
               'LeptonJet':'FullRazorInclusive_SingleLepton_2016B_PRv2_GoodLumiGoldenJun16_NoDuplicates.root',
               }

    
    dateString = str(datetime.date.today()).replace("-","_")
    #dateString = '2015_11_17'
    dateString = '2016_08_23'
    
    dryRun = False
    
    for box in boxes:
        for fit in fits:
            btag = '0-3btag'
            if box in ['DiJet','LeptonJet']:
                btag = '0-2btag'
            
            fitString = '--fit-region Full'
            plotString = '--plot-region Full'
            if fit=='Sideband':
                fitString = '--fit-region LowMR,LowRsq'
                plotString = '--plot-region LowMR,LowRsq'
            #exec_me('python python/DustinTuple2RooDataSet.py -b %s -c %s -d Datasets/ Run2016BCD/%s --data -l %i'% (box, config, dataset[box], lumi ), dryRun )
            outDir = 'fits_%s/Run2016BCD/%s/%s/' % (dateString, box, fit)
            roodataset = dataset[box].replace('.root','_lumi-%.3f_%s_%s.root'%(lumi/1000.,btag,box))
            #exec_me('python python/BinnedFit.py -b %s -c %s -d %s Datasets/%s --data -l %i %s %s'% (box, config, outDir, roodataset, lumi, fitString, plotString), dryRun )
        
            #exec_me('python python/PlotFit.py -b %s -c %s -d %s -i %s/BinnedFitResults_%s.root --data -l %i'% (box, config, outDir, outDir, box, lumi), dryRun )
            varyNString = ''
            #if fit=='Sideband':
            #    varyNString = '--vary-n'
                
            #exec_me('python python/RunToys.py -c %s -d %s -b %s -l %i -i %s/BinnedFitResults_%s.root %s --no-stat -t 10000 %s' %(config,outDir,box,lumi,outDir,box,varyNString,fitString), dryRun )
            #exec_me('python python/RunToys.py -c %s -d %s -b %s -l %i -i %s/BinnedFitResults_%s.root %s -t 10000 %s' %(config,outDir,box,lumi,outDir,box,varyNString,fitString), dryRun )
            exec_me('python python/PlotFit.py -b %s -c %s -d %s -i %s/BinnedFitResults_%s.root --data -l %i -t %s/toys_Bayes_%s.root -s %s/toys_Bayes_noStat_%s.root %s --no-stat'% (box, config, outDir, outDir, box, lumi, outDir, box, outDir,box,fitString), dryRun )
            #exec_me('python python/PlotFit.py -b %s -c %s -d %s -i %s/BinnedFitResults_%s.root --data -l %i -t %s/toys_Bayes_%s.root -s %s/toys_Bayes_noStat_%s.root %s --no-stat --print-errors %s'% (box, config, outDir, outDir, box, lumi, outDir, box, outDir,box,fitString,plotString), dryRun )
    
            #exec_me('python python/RunToys.py -c %s -d %s -b %s -l %i -i %s/BinnedFitResults_%s.root --freq -t 200 -s 0 %s' %(config,outDir,box,lumi,outDir,box,fitString), dryRun )            
            #exec_me('python python/PlotGOF.py -c %s -d %s -b %s -l %f -t %s/toys_Freq_s0_%s.root %s --data' %(config,outDir,box,lumi,outDir,box,fitString),dryRun)
