from optparse import OptionParser
import os
import ROOT as rt
from array import *
from framework import Config
import sys
import glob

def exec_me(command,dryRun=True):
    print command
    if not dryRun: os.system(command)

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('--mGluino',dest="mGluino", default=-1,type="float",
                  help="mass of gluino")
    parser.add_option('--mStop',dest="mStop", default=-1,type="float",
                  help="mass of stop")
    parser.add_option('--mLSP',dest="mLSP", default=-1,type="float",
                  help="mass of LSP")
    parser.add_option('-l','--lumi-array',dest="lumi_array", default="0.2,0.5,1,3,4,7,10",type="string",
                  help="lumi array in fb^-1, e.g.: 0.2,0.5,1,3,4,7,10")
    parser.add_option('--signif',dest="signif",default=False,action='store_true',
                  help="calculate significance instead of limit")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('--fit',dest="fit",default=False,action='store_true',
                  help="Turn on pre-fit")
    parser.add_option('--no-fit',dest="noFit",default=False,action='store_true',
                  help="no fit, just use MC")
    parser.add_option('--min-tol',dest="min_tol",default=0.001,type="float",
                  help="minimizer tolerance (default = 0.001)")
    parser.add_option('--min-strat',dest="min_strat",default=2,type="int",
                  help="minimizer strategy (default = 2)")
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                  help="Just print out commands to run")
    parser.add_option('-u','--unweighted',dest="unweighted",default=False,action='store_true',
                  help="use unweighted dataset")
    parser.add_option('--penalty',dest="penalty",default=False,action='store_true',
                  help="penalty terms on background parameters")
    parser.add_option('--data',dest="isData", default=False,action='store_true',
                  help="changes plots for data")
    parser.add_option('--just-combo',dest="justCombo", default=True,action='store_false',
                  help="just do the combination, on by default")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default='FitResults/BinnedFitResults.root',type="string",
                  help="input fit file")
    parser.add_option('--no-signal-sys',dest="noSignalSys",default=False,action='store_true',
                  help="no signal shape systematic uncertainties")
    parser.add_option('--rMax',dest="rMax",default=-1,type="float",
                  help="maximum r value (for better precision)")
    parser.add_option('--histo-file',dest="histoFile", default=None,type="string",
                  help="input histogram file for MADD/fit systematic")
    #pdf uncertainty options.  current prescription is just to take 10% uncorrelated error on each bin
    #parser.add_option('--num-pdf-weights',dest="numPdfWeights",default=0,type='int',
    #              help='number of pdf nuisance parameters to use')
    #parser.add_option('--compute-pdf-envelope',dest="computePdfEnvelope",default=False,action='store_true',
    #              help="Use the SUS pdf reweighting prescription, summing weights in quadrature")


    (options,args) = parser.parse_args()

    boxes = options.box.split('_')

    signif = options.signif
    
    lumiArray = array('d',[float(lumi) for lumi in options.lumi_array.split(',')])
    cfg = Config.Config(options.config)
    
    model = options.model
    if options.mGluino>-1:
        #massPoint = 'mGl-%i_mLSP-%i'%(options.mGluino,options.mLSP)
        massPoint = '%i_%i'%(options.mGluino,options.mLSP)
    elif options.mStop>-1:
        #massPoint = 'mStop-%i_mLSP-%i'%(options.mStop,options.mLSP)
        massPoint = '%i_%i'%(options.mStop,options.mLSP)

    signalSys = ''
    if options.noSignalSys:
        signalSys = '--no-signal-sys'
    
    fit = ''
    if options.fit:
        fit = '--fit'
    elif options.noFit:
        fit = '--no-fit'
    
    penaltyString = ''
    if options.penalty: penaltyString = '--penalty'
        
    histoString = ''
    if options.histoFile:
        histoString = '--histo-file %s'%(options.histoFile)
        
    dataset = {#'MultiJet':'RazorInclusive_HTMHT_Run2015D_2093pb_GoodLumiGolden_RazorSkim_Filtered',
               'MultiJet':'RazorInclusive_HTMHT_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter',
               'MuMultiJet':'RazorInclusive_SingleMuon_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter',
               'EleMultiJet':'RazorInclusive_SingleElectron_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter',
               'DiJet':'RazorInclusive_HTMHT_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter',
               'LooseLeptonMultiJet':'RazorInclusive_HTMHT_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter'
               }

    eosLocationSMS = {#'T1bbbb': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p21_ForFullStatus20151030/jobs/combined/',
                      #'T1tttt': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/MC/combined/',
                      #'T1bbbb': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForApproval20151208/jobs/combined/',
                      'T1bbbb': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined/',
                      'T1tttt': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined/',
                      'T5ttttDM175': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined/',
                      'T5ttttDM175T2tt': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined/',
                      'T1ttbb': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined/',
                      'T1qqqq': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined/',
                      'T2tt': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined/',
                      'T5qqqqVV': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined/',
                    }

    #eosLocationData = 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/'
    eosLocationData = 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForMoriond20160119/RazorSkim/'

    eosLocationBkg =  'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/'
    
    bkgList = ['RazorInclusive_DYJetsToLL_M-50_HTBinned_1pb_weighted.root',
               'RazorInclusive_DYJetsToLL_M-5to50_HTBinned_1pb_weighted.root',
               'RazorInclusive_SingleTop_1pb_weighted.root',
               'RazorInclusive_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root',
               'RazorInclusive_TTV_1pb_weighted.root'
               'RazorInclusive_VV_1pb_weighted.root',
               'RazorInclusive_WJetsToLNu_HTBinned_1pb_weighted.root',
               'RazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted.root']

    bkgString = ' '.join(['%s/%s'%(eosLocationBkg,bkg) for bkg in bkgList])
    
    exec_me('mkdir -p Datasets',options.dryRun)        
    exec_me('mkdir -p %s'%options.outDir,options.dryRun)
    for lumi in lumiArray:
        for box in boxes:            
            z = array('d', cfg.getBinning(box)[2]) # nBtag binning
            btagMin = z[0]
            btagMax = z[-1]
            if btagMax-1>btagMin:          
                btag = '%i-%ibtag'%(btagMin,btagMax-1)
            else:
                btag = '%ibtag'%(btagMin)    


            #signalDsName = 'Datasets/RazorInclusive_SMS-%s_2J_%s_weighted_lumi-%.3f_%s_%s.root'%(model,massPoint,lumi,btag,box)
            signalDsName = 'Datasets/SMS-%s_%s_lumi-%.3f_%s_%s.root'%(model,massPoint,lumi,btag,box)
            #exec_me('python python/DustinTuple2RooDataSet.py -c %s -b %s -d Datasets/ -w Signals/SMS-%s_%s.root -l %f'%(options.config,box,model,massPoint, 1000*lumi),options.dryRun)
            
            #comment out old pdf uncertainty stuff for now -- keep it here in case SUS changes its mind again.
            #computePdfEnvelopeString = ''
            #if options.computePdfEnvelope:
            #    computePdfEnvelopeString = '--compute-pdf-envelope'
            #exec_me('python python/SMSTemplates.py -c %s -b %s -d Datasets/ %s/SMS-%s_%s.root --num-pdf-weights %d %s -l %f %s'%(options.config,box,eosLocationSMS[model],model,massPoint, options.numPdfWeights, computePdfEnvelopeString, 1000*lumi,signalSys),options.dryRun)
            brString = ''
            if 'T1x' in model:
                xBR = float(model[model.find('x')+1:model.find('y')].replace('p','.'))
                yBR = float(model[model.find('y')+1:].replace('p','.'))
                brString = '--xBR %.2f --yBR %.2f'%(xBR,yBR)
                exec_me('python python/SMSTemplates.py -c %s -b %s -d Datasets/ %s/SMS-%s_%s.root -l %f %s %s'%(options.config,box,eosLocationSMS['T1ttbb'],'T1ttbb',massPoint,1000*lumi,signalSys,brString),options.dryRun)
            elif (model.split('p'))>1:
                signals = model.split('p')
                massPoint2 = ''
                if signals[1]=='T2tt':
                    if int(options.mLSP)<=375:
                        massPoint2 = '%i_%i'%(options.mLSP+175,options.mLSP)
                    else:
                        massPoint2 = '%i_%i'%(550,375)
                exec_me('python python/SMSTemplates.py -c %s -b %s -d Datasets/ %s/SMS-%s_%s.root -l %f %s %s'%(options.config,box,eosLocationSMS[signals[0]],signals[0],massPoint,1000*lumi,signalSys,brString),options.dryRun)
                exec_me('python python/SMSTemplates.py -c %s -b %s -d Datasets/ %s/SMS-%s_%s.root -l %f %s %s --mStop %i --mLSP %i'%(options.config,box,eosLocationSMS[signals[1]],signals[1],massPoint2,1000*lumi,signalSys,brString,options.mLSP+175,options.mLSP),options.dryRun)
                # after overriding massPoint name this should be the dataset's new massPoint:
                massPoint2 = '%i_%i'%(options.mLSP+175,options.mLSP)                
                signal1DsName = 'Datasets/SMS-%s_%s_lumi-%.3f_%s_%s.root'%(signals[0],massPoint,lumi,btag,box)
                signal2DsName = 'Datasets/SMS-%s_%s_lumi-%.3f_%s_%s.root'%(signals[1],massPoint2,lumi,btag,box)
                exec_me('hadd -f %s %s %s'%(signalDsName,signal1DsName,signal2DsName),options.dryRun)
            else:                
                exec_me('python python/SMSTemplates.py -c %s -b %s -d Datasets/ %s/SMS-%s_%s.root -l %f %s'%(options.config,box,eosLocationSMS[model],model,massPoint,1000*lumi,signalSys),options.dryRun)
                
            if options.isData:
                backgroundDsName = 'Datasets/%s_lumi-%.3f_%s_%s.root'%(dataset[box],lumi,btag,box)
                exec_me('python python/DustinTuple2RooDataSet.py -b %s -c %s -d Datasets/ %s/%s.root --data -l %f'% (box, options.config, eosLocationData, dataset[box], 1000*lumi), options.dryRun )
            else:                
                backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.3f_%s_%s.root'%(lumi,btag,box)
                exec_me('python python/DustinTuple2RooDataSet.py -c %s -b %s -d Datasets/ -w %s -l %f'%(options.config, box, bkgString, 1000*lumi),options.dryRun)
            
                if options.unweighted:                
                    backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_unweighted_lumi-%.3f_%s_%s.root'%(lumi,btag,box)
                    exec_me('python python/RooDataSet2UnweightedDataSet.py -c %s -b %s -d Datasets/ Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.3f_%s_%s.root'%(options.config,box,lumi,btag,box),options.dryRun)

                
            #comment out old pdf uncertainty stuff for now -- keep it here in case SUS changes its mind again.
            #exec_me('python python/WriteDataCard.py --num-pdf-weights %d %s -i %s -l %f -c %s -b %s -d %s %s %s %s %s %s'%(options.numPdfWeights, computePdfEnvelopeString, options.inputFitFile,1000*lumi,options.config,box,options.outDir,fit,signalDsName,backgroundDsName,penaltyString,signalSys),options.dryRun)
            exec_me('python python/WriteDataCard.py -i %s -l %f -c %s -b %s -d %s %s %s %s %s %s %s'%(options.inputFitFile,1000*lumi,options.config,box,options.outDir,fit,signalDsName,backgroundDsName,penaltyString,signalSys,histoString),options.dryRun)
            
            if (len(boxes)>1 and not options.justCombo) or len(boxes)==1:
                if signif:
                    exec_me('combine -M ProfileLikelihood --signif --expectSignal=1 -t -1 --toysFreq %s/razor_combine_%s_%s_lumi-%.3f_%s.txt -n %s_%s_lumi-%.3f_%s'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,box),options.dryRun)
                    exec_me('mv higgsCombine%s_%s_lumi-%.3f_%s.ProfileLikelihood.mH120.root %s/'%(model,massPoint,lumi,box,options.outDir),options.dryRun)
                else:
                    if options.rMax>-1:
                        exec_me('combine -M Asymptotic %s/razor_combine_%s_%s_lumi-%.3f_%s_%s.txt -n %s_%s_lumi-%.3f_%s_%s --minimizerTolerance %f --minimizerStrategy %i --setPhysicsModelParameterRanges r=0,%f --saveWorkspace'%(options.outDir,model,massPoint,lumi,btag,box,model,massPoint,lumi,btag,box,options.min_tol,options.min_strat,options.rMax),options.dryRun)
                    else:
                        exec_me('combine -M Asymptotic %s/razor_combine_%s_%s_lumi-%.3f_%s_%s.txt -n %s_%s_lumi-%.3f_%s_%s --minimizerTolerance %f --minimizerStrategy %i --saveWorkspace'%(options.outDir,model,massPoint,lumi,btag,box,model,massPoint,lumi,btag,box,options.min_tol,options.min_strat),options.dryRun)
                    exec_me('mv higgsCombine%s_%s_lumi-%.3f_%s_%s.Asymptotic.mH120.root %s/'%(model,massPoint,lumi,btag,box,options.outDir),options.dryRun)                
        if len(boxes)>1:
            for box in boxes: exec_me('cp %s/razor_combine_%s_%s_lumi-%.3f_%s_%s.txt .'%(options.outDir,model,massPoint,lumi,btag,box),options.dryRun)
            cmds = ['%s=razor_combine_%s_%s_lumi-%.3f_%s_%s.txt'%(box,model,massPoint,lumi,btag,box) for box in boxes]
            exec_me('combineCards.py %s > %s/razor_combine_%s_%s_lumi-%.3f_%s_%s.txt'%(' '.join(cmds),options.outDir,model,massPoint,lumi,btag,options.box),options.dryRun)
            exec_me('cat %s/razor_combine_%s_%s_lumi-%.3f_%s_%s.txt'%(options.outDir,model,massPoint,lumi,btag,options.box),options.dryRun)
            if signif:
                exec_me('combine -M ProfileLikelihood --signif --expectSignal=1 -t -1 --toysFreq %s/razor_combine_%s_%s_lumi-%.3f_%s.txt -n %s_%s_lumi-%.3f_%s_%s'%(options.outDir,model,massPoint,lumi,options.box,model,massPoint,lumi,btag,options.box),options.dryRun)
                exec_me('mv higgsCombine%s_%s_lumi-%.3f_%s_%s.ProfileLikelihood.mH120.root %s/'%(model,massPoint,lumi,btag,options.box,options.outDir),options.dryRun)
            else:
                if options.rMax>-1:
                    exec_me('combine -M Asymptotic %s/razor_combine_%s_%s_lumi-%.3f_%s_%s.txt -n %s_%s_lumi-%.3f_%s_%s --minimizerTolerance %f --minimizerStrategy %i --setPhysicsModelParameterRanges r=0,%f --saveWorkspace'%(options.outDir,model,massPoint,lumi,btag,options.box,model,massPoint,lumi,btag,options.box,options.min_tol,options.min_strat,options.rMax),options.dryRun)
                else:
                    exec_me('combine -M Asymptotic %s/razor_combine_%s_%s_lumi-%.3f_%s_%s.txt -n %s_%s_lumi-%.3f_%s_%s --minimizerTolerance %f --minimizerStrategy %i --saveWorkspace'%(options.outDir,model,massPoint,lumi,btag,options.box,model,massPoint,lumi,btag,options.box,options.min_tol,options.min_strat),options.dryRun)
                exec_me('mv higgsCombine%s_%s_lumi-%.3f_%s_%s.Asymptotic.mH120.root %s/'%(model,massPoint,lumi,btag,options.box,options.outDir),options.dryRun)
            for box in boxes: exec_me('rm razor_combine_%s_%s_lumi-%.3f_%s_%s.txt'%(model,massPoint,lumi,btag,box),options.dryRun)
 
