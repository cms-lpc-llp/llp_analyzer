##Weight and trigger utilities for inclusive razor analysis

import math
import numpy as np
import ROOT as rt
import sys

#####################################
### WEIGHT AND TRIGGER INFO
#####################################

def integral_error(fun, mean, cov, *coords):
    n_toys = 1000
    toys = np.zeros(n_toys)
    for itoy in range(n_toys):
        rand = np.random.multivariate_normal(mean, cov)
        fun.SetParameters(rand)
        integral = fun.Integral(*coords)
        toys[itoy] = integral
    fun.SetParameters(mean)
    return toys.std()

def makeQCDBtagSysHists(hists, whists, region, nbtags):
    """
    Make new up/down shape histograms with QCD varied up/down
    according to the b-tag extrapolation uncertainty
    """
    region = region.lower()
    btag_bin = min(4, nbtags+1)
    sf_central = whists['qcdbtags'+region].GetBinContent(btag_bin)
    sf_unc = whists['qcdbtags'+region+'sys'].GetBinContent(btag_bin)
    hists['Sys']['QCD']['qcdbtagsysUp'] = {}
    hists['Sys']['QCD']['qcdbtagsysDown'] = {}
    for var in hists['QCD']:
        center_hist = hists['QCD'][var]
        up_hist = center_hist.Clone(
                center_hist.GetName()+'qcdbtagsysUp')
        down_hist = center_hist.Clone(
                center_hist.GetName()+'qcdbtagsysDown')
        sf_up = (sf_central + sf_unc) / sf_central
        sf_down = max(0, sf_central - sf_unc) / sf_central
        print "Scaling qcdbtagsysUp hist by",sf_up
        print "Scaling qcdbtagsysDown hist by",sf_down
        up_hist.Scale(sf_up)
        down_hist.Scale(sf_down)
        hists['Sys']['QCD']['qcdbtagsysUp'][var] = up_hist
        hists['Sys']['QCD']['qcdbtagsysDown'][var] = down_hist

def getQCDExtrapolationFactor(mrlow, mrhigh, rsqlow, rsqhigh, nbtags,
        fun, mean, cov, btagHist, errorOpt, debugLevel, alt_fun=None):
    """Get QCD extrapolation factor as a function of MR and Rsq"""
    btagBin = min(4, nbtags+1)
    norm = btagHist.GetBinContent(btagBin)
    area = (mrhigh-mrlow)*(rsqhigh-rsqlow)
    integral = fun.Integral(mrlow, mrhigh, rsqlow, rsqhigh)
    sf = norm*integral/area

    # An alternate extrapolation function can be used
    # to quantify uncertainty on the choice of functional form
    if alt_fun is not None:
        alt_integral = alt_fun.Integral(mrlow, mrhigh, 
                rsqlow, rsqhigh)
        alt_sf = norm * alt_integral / area
        alt_sf_frac = alt_sf / sf

    if errorOpt is not None:
        sfErr = norm * integral_error(fun, mean, cov, 
                mrlow, mrhigh, rsqlow, rsqhigh) / area
        normErr = btagHist.GetBinError(btagBin) * integral / area
        if errorOpt == 'qcdnormUp':
            return sf + sfErr
        elif errorOpt == 'qcdnormDown':
            return max(0., sf - sfErr)
        elif errorOpt == 'qcdbtagUp':
            return sf + normErr
        elif errorOpt == 'qcdbtagDown':
            return max(0., sf - normErr)
        elif errorOpt == 'qcdaltfunctionUp':
            return sf * alt_sf_frac
        elif errorOpt == 'qcdaltfunctionDown':
            return sf / alt_sf_frac

    return sf

def read_mean_and_cov(in_name):
    """
    Retrieves mean vector and covariance matrix
    from first two lines of npz file and returns them
    as numpy arrays.
    """
    with open(in_name, 'r') as in_f:
        arrs = np.load(in_f)
        mean = arrs['arr_0']
        cov = arrs['arr_1']
    return mean, cov

def makeQCDExtrapolation(qcd2DHist, wHists, region, btags,
        debugLevel=0, errOpt=None):
    """Performs extrapolation from high to low delta phi region
        for QCD prediction in the signal region.
        Returns corrected histogram with same shape as qcd2DHist"""
    region = region.lower()
    newHist = qcd2DHist.Clone()
    newHist.SetDirectory(0)
    # Sorry, I am hard-coding this path until I have time 
    # to implement a more flexible system for this
    npz_file = 'macros/BackgroundStudies/QCD/qcd_best_fit_2d_{}.npz'.format(region)
    mean, cov = read_mean_and_cov(npz_file)
    fun = wHists['qcdfunction{}'.format(region)]
    alt_fun = None
    if errOpt is not None and 'qcdaltfunction' in errOpt:
        alt_fun = wHists['qcdaltfunction{}'.format(region)]
    btagHist = wHists['qcdbtags{}'.format(region)]
    for bx in range(1, qcd2DHist.GetNbinsX()+1):
        for by in range(1, qcd2DHist.GetNbinsY()+1):
            mrlow = qcd2DHist.GetXaxis().GetBinLowEdge(bx) 
            mrhigh = qcd2DHist.GetXaxis().GetBinUpEdge(bx) 
            rsqlow = qcd2DHist.GetYaxis().GetBinLowEdge(by) 
            rsqhigh = qcd2DHist.GetYaxis().GetBinUpEdge(by) 
            extrapFactor = getQCDExtrapolationFactor(mrlow, mrhigh, rsqlow, rsqhigh, btags,
                    fun, mean, cov, btagHist, errOpt, debugLevel,
                    alt_fun=alt_fun)
            if debugLevel > 0:
                print "QCD extrapolation factor in bin {} {} {} {} is {} (error opt {})".format(
                        mrlow, mrhigh, rsqlow, rsqhigh, extrapFactor, errOpt)
            newHist.SetBinContent(bx, by, qcd2DHist.GetBinContent(
                bx, by) * extrapFactor)
            newHist.SetBinError(bx, by, qcd2DHist.GetBinError(
                bx, by) * extrapFactor)
    return newHist

def getNISRScaleFactor(nisr, d=1.071):
    """See comments at getNISRCorrection()"""
    if nisr == 0:
        return d
    if nisr == 1:
        return d * 0.920
    if nisr == 2:
        return d * 0.821
    if nisr == 3:
        return d * 0.715
    if nisr == 4:
        return d * 0.662
    if nisr == 5:
        return d * 0.561
    return d * 0.511

def getNISRCorrection(event, d=1.071):
    """Computes a weight according to the number of ISR jets in the (ttbar) event.
        Details at https://indico.cern.ch/event/592621/contributions/2398559/
        attachments/1383909/2105089/16-12-05_ana_manuelf_isr.pdf"""
    n = event.NISRJets
    return getNISRScaleFactor(n, d)

def passTrigger(event, triggerNumList):
    """Checks if the event passed any trigger in the list"""
    passes = False
    for i in triggerNumList:
        if event.HLTDecision[i]:
            passes = True
            break
    return passes

def pileupWeight(event, wHists, puBranch="NPV", debugLevel=0):
    if "pileup" in wHists:
        if not hasattr(event, puBranch):
            print "Error in pileupWeight: tree does not have a branch for ",puBranch
            sys.exit()
        pileupWeight = wHists["pileup"].GetBinContent(wHists["pileup"].GetXaxis().FindFixBin(getattr(event, puBranch)))
        if debugLevel > 1: print puBranch," weight:",pileupWeight,"(",puBranch,"=",getattr(event, puBranch),")"
        return pileupWeight
    else:
        print "Error in pileupWeight: pileup reweighting histogram not found!"
        sys.exit()

def reapplyPileupWeight(event, wHists, weightBranch="pileupWeight", debugLevel=0):
    """Get (pileupweight)/(oldpileupweight)"""
    if "pileup" in wHists:
        #use NPU_0 or NPU, whichever exists in the tree
        if not hasattr(event, "NPU_0"):
            if not hasattr(event, "NPU"):
                sys.exit( "Error in pileupWeight: tree does not have a branch for NPU_0 or NPU" )
            else:
                puBranch = "NPU"
        else:
            puBranch = "NPU_0"
        if not hasattr(event, weightBranch):
            sys.exit( "Error in pileupWeight: tree does not have a branch for "+weightBranch )
        pileupWeight = wHists["pileup"].GetBinContent(wHists["pileup"].GetXaxis().FindFixBin(
            getattr(event, puBranch))) / getattr(event,weightBranch)
        if debugLevel > 1: print puBranch," weight:",pileupWeight,"(",puBranch,"=",getattr(event, puBranch),")"
        return pileupWeight
    else:
        print "Error in pileupWeight: pileup reweighting histogram not found!"
        sys.exit()

def reapplyLepEffWeight(event, wHists, eleBranch="eleEffWeight", muBranch="muonEffWeight", debugLevel=0):
    """Get (new lepton weight)/(old lepton weight)"""
    if abs(event.lep1Type) == 13:
        if "muoneff" in wHists:
            if not hasattr(event, muBranch):
                sys.exit( "Error in reapplyLepEffWeight: tree does not have a branch for "+muBranch )
            muWeight = wHists["muoneff"].GetBinContent(
                    wHists["muoneff"].GetXaxis().FindFixBin( max( min(event.lep1.Pt(), 199.9), 20.01 ) ),
                    wHists["muoneff"].GetYaxis().FindFixBin( abs( event.lep1.Eta() ) ) ) / getattr( event, muBranch )
            if debugLevel > 1: 
                print "Muon weight:",muWeight,"(",muBranch,"=",getattr(event, muBranch),")"
            return muWeight
        else:
            print "Error: muon efficiency reweighting histogram not found!"
            sys.exit()
    elif abs(event.lep1Type) == 11:
        if "eleeff" in wHists:
            if not hasattr(event, eleBranch):
                sys.exit( "Error in reapplyLepEffWeight: tree does not have a branch for "+eleBranch )
            eleWeight = wHists["eleeff"].GetBinContent(
                    wHists["eleeff"].GetXaxis().FindFixBin( max( min(event.lep1.Pt(), 199.9), 20.01 ) ),
                    wHists["eleeff"].GetYaxis().FindFixBin( abs( event.lep1.Eta() ) ) ) / getattr( event, eleBranch )
            if debugLevel > 1: 
                print "Electron weight:",eleWeight,"(",eleBranch,"=",getattr(event, eleBranch),")"
            return eleWeight
        else:
            print "Error: electron efficiency reweighting histogram not found!"
            sys.exit()

def reapplyLepTrigWeight(event, wHists, trigBranch="trigWeight1L", debugLevel=0):
    """Get (new trig weight)/(old trig weight)"""
    if not hasattr(event, trigBranch):
        sys.exit( "Error in reapplyLepTrigWeight: tree does not have a branch for "+trigBranch )
    oldWeight = getattr( event, trigBranch )
    if abs(event.lep1Type) == 13:
        if "muontrig" in wHists:
            newWeight = wHists["muontrig"].GetBinContent(
                    wHists["muontrig"].GetXaxis().FindFixBin( abs( event.lep1.Eta() ) ),
                    wHists["muontrig"].GetYaxis().FindFixBin( max( min(event.lep1.Pt(), 199.9), 20.01 ) ) ) 
        else:
            print "Error: muon trigger reweighting histogram not found!"
            sys.exit()
    elif abs(event.lep1Type) == 11:
        if "eletrig" in wHists:
            newWeight = wHists["eletrig"].GetBinContent(
                    wHists["eletrig"].GetXaxis().FindFixBin( abs( event.lep1.Eta() ) ),
                    wHists["eletrig"].GetYaxis().FindFixBin( max( min(event.lep1.Pt(), 199.9), 20.01 ) ) )
        else:
            print "Error: electron trigger reweighting histogram not found!"
            sys.exit()
    if debugLevel > 1: 
        print "Trigger weight:",newWeight,"( old weight =",oldWeight,")"
    return newWeight / oldWeight

def getNBJetsWeight(event, debugLevel=0):
    """Reweight according to number of gen-level b-jets"""
    sfTag = 0.95
    sfNoTag = (1/.68 - sfTag) / (1/.68 - 1) #using 0.68 as the average b-tag efficiency
    bWeight = 1.0
    nTags = event.NBJetsMedium
    nGenB = event.NGenBJets
    for nb in range(nGenB):
        if nTags > nb: #we successfully tagged this many jets
            bWeight *= sfTag
        else: #we didn't tag this many jets
            bWeight *= sfNoTag
    if debugLevel > 1:
        print "Found",nGenB,"gen b-jets in event and tagged",nTags,"; assigning a weight of",bWeight
    return bWeight

def leptonWeight(event, wHists, doLep2=False, debugLevel=0):
    """Weights using leading lepton pt and eta.  If doLep2 is True, weights using subleading lepton too."""
    #check for histograms
    if not ("ele" in wHists):
        print "Error in leptonWeight: electron scale factor histogram not found!"
        sys.exit()
    if not ("muon" in wHists):
        print "Error in leptonWeight: muon scale factor histogram not found!"
        sys.exit()

    leptonWeight = 1.0
    #leading muon 
    if abs(event.lep1Type) == 13 and event.lep1.Pt() > 0:
        leptonWeight *= wHists["muon"].GetBinContent(
                wHists["muon"].GetXaxis().FindFixBin(max(min(event.lep1.Pt(), 199.9),15.01)),
                #wHists["muon"].GetXaxis().FindFixBin(min(event.lep1.Pt(), wHists["muon"].GetXaxis().GetXmax()-1.0)),
                wHists["muon"].GetYaxis().FindFixBin(abs(event.lep1.Eta())))
    #leading electron
    elif abs(event.lep1Type) == 11 and event.lep1.Pt() > 0:
        leptonWeight *= wHists["ele"].GetBinContent(
                wHists["ele"].GetXaxis().FindFixBin(max(min(event.lep1.Pt(), 199.9),15.01)),
                #wHists["ele"].GetXaxis().FindFixBin(min(event.lep1.Pt(), wHists["ele"].GetXaxis().GetXmax()-1.0)),
                wHists["ele"].GetYaxis().FindFixBin(abs(event.lep1.Eta())))
    if doLep2:
        #subleading muon 
        if abs(event.lep2Type) == 13 and event.lep2.Pt() > 0:
            leptonWeight *= wHists["muon"].GetBinContent(
                    wHists["muon"].GetXaxis().FindFixBin(max(min(event.lep2.Pt(), 199.9),15.01)),
                    #wHists["muon"].GetXaxis().FindFixBin(min(event.lep2.Pt(), wHists["muon"].GetXaxis().GetXmax()-1.0)),
                    wHists["muon"].GetYaxis().FindFixBin(abs(event.lep2.Eta())))
        #subleading electron
        elif abs(event.lep2Type) == 11 and event.lep2.Pt() > 0:
            leptonWeight *= wHists["ele"].GetBinContent(
                    wHists["ele"].GetXaxis().FindFixBin(max(min(event.lep2.Pt(), 199.9),15.01)),
                    #wHists["ele"].GetXaxis().FindFixBin(min(event.lep2.Pt(), wHists["ele"].GetXaxis().GetXmax()-1.0)),
                    wHists["ele"].GetYaxis().FindFixBin(abs(event.lep2.Eta())))
    if debugLevel > 1: 
        print "lepton weight:",leptonWeight
    return leptonWeight

def leptonTriggerWeight(event, wHists, doLep2=False, debugLevel=0):
    if not ("eletrig" in wHists):
        print "Error in leptonTriggerWeight: electron trigger scale factor histogram not found!"
        sys.exit()
    if not ("muontrig" in wHists):
        print "Error in leptonTriggerWeight: muon trigger scale factor histogram not found!"
        sys.exit()
    trigWeight = 1.0
    #leading muon 
    if abs(event.lep1Type) == 13 and event.lep1.Pt() > 0:
        trigWeight *= wHists["muontrig"].GetBinContent(
                wHists["muontrig"].GetXaxis().FindFixBin(max(min(event.lep1.Pt(), 199.9),15.01)),
                #wHists["muontrig"].GetXaxis().FindFixBin(min(event.lep1.Pt(), wHists["muontrig"].GetXaxis().GetXmax()-1.0)),
                wHists["muontrig"].GetYaxis().FindFixBin(abs(event.lep1.Eta())))
    #leading electron
    elif abs(event.lep1Type) == 11 and event.lep1.Pt() > 0:
        trigWeight *= wHists["eletrig"].GetBinContent(
                wHists["eletrig"].GetXaxis().FindFixBin(max(min(event.lep1.Pt(), 199.9),15.01)),
                #wHists["eletrig"].GetXaxis().FindFixBin(min(event.lep1.Pt(), wHists["eletrig"].GetXaxis().GetXmax()-1.0)),
                wHists["eletrig"].GetYaxis().FindFixBin(abs(event.lep1.Eta())))
    if doLep2: 
        #subleading muon 
        if abs(event.lep2Type) == 13 and event.lep2.Pt() > 0:
            trigWeight *= wHists["muontrig"].GetBinContent(
                    wHists["muontrig"].GetXaxis().FindFixBin(max(min(event.lep2.Pt(), 199.9),15.01)),
                    #wHists["muontrig"].GetXaxis().FindFixBin(min(event.lep2.Pt(), wHists["muontrig"].GetXaxis().GetXmax()-1.0)),
                    wHists["muontrig"].GetYaxis().FindFixBin(abs(event.lep2.Eta())))
        #subleading electron
        elif abs(event.lep2Type) == 11 and event.lep2.Pt() > 0:
            trigWeight *= wHists["eletrig"].GetBinContent(
                    wHists["eletrig"].GetXaxis().FindFixBin(max(min(event.lep2.Pt(), 199.9),15.01)),
                    #wHists["eletrig"].GetXaxis().FindFixBin(min(event.lep2.Pt(), wHists["eletrig"].GetXaxis().GetXmax()-1.0)),
                    wHists["eletrig"].GetYaxis().FindFixBin(abs(event.lep2.Eta())))
    if debugLevel > 1: 
        if not doLep2:
            print "1-lepton trigger weight:",trigWeight
        else: 
            print "2-lepton trigger weight:",trigWeight
    return trigWeight

def getMTForTTBarDilepton(event, lepID):
    """Computes the MT from the specified lepton (1 or 2, according to lepID) and the MET"""
    met = event.MET
    phi = event.METPhi
    if lepID == 1:
        lep = event.lep1
    elif lepID == 2:
        lep = event.lep2
    else:
        raise ValueError("invalid lepton ID")
    return ( lep.M2() + 2*met*lep.Pt()*( 1 - math.cos(math.acos(math.cos(phi-lep.Phi()))) ) )**(0.5)

def getTTBarDileptonWeight(event):
    """Returns 0.5 times the number of leptons passing the MT cut"""
    weight = 0
    #lepton 1
    lep1MT = getMTForTTBarDilepton(event, 1)
    if lep1MT > 120:
        weight += 0.5
    #lepton 2
    lep2MT = getMTForTTBarDilepton(event, 2)
    if lep2MT > 120:
        weight += 0.5
    return weight

weightTable = {
        "tightmuoneffUp":"sf_muonEffUp",
        "tightmuoneffDown":"sf_muonEffDown",
        "tighteleeffUp":"sf_eleEffUp",
        "tighteleeffDown":"sf_eleEffDown",
        "vetomuoneffUp":"sf_vetoMuonEffUp",
        "vetomuoneffDown":"sf_vetoMuonEffDown",
        "vetoeleeffUp":"sf_vetoEleEffUp",
        "vetoeleeffDown":"sf_vetoEleEffDown",
        "tightmuonfastsimUp":"sf_muonEffFastsimSFUp",
        "tightmuonfastsimDown":"sf_muonEffFastsimSFDown",
        "tightelefastsimUp":"sf_eleEffFastsimSFUp",
        "tightelefastsimDown":"sf_eleEffFastsimSFDown",
        "vetomuonfastsimUp":"sf_vetoMuonEffFastsimSFUp",
        "vetomuonfastsimDown":"sf_vetoMuonEffFastsimSFDown",
        "vetoelefastsimUp":"sf_vetoEleEffFastsimSFUp",
        "vetoelefastsimDown":"sf_vetoEleEffFastsimSFDown",
        "muontrigUp":"sf_muonTrigUp",
        "muontrigDown":"sf_muonTrigDown",
        "eletrigUp":"sf_eleTrigUp",
        "eletrigDown":"sf_eleTrigDown",
        "btagUp":"sf_btagUp",
        "btagDown":"sf_btagDown",
        "bmistagUp":"sf_bmistagUp",
        "bmistagDown":"sf_bmistagDown",
        "btagfastsimUp":"sf_btagFastsimSFUp",
        "btagfastsimDown":"sf_btagFastsimSFDown",
        "pileupUp":"pileupWeightUp",
        "pileupDown":"pileupWeightDown",
        "isrUp":"ISRSystWeightUp",
        "isrDown":"ISRSystWeightDown",
        "facscaleUp":"sf_facScaleUp",
        "facscaleDown":"sf_facScaleDown",
        "renscaleUp":"sf_renScaleUp",
        "renscaleDown":"sf_renScaleDown",
        "facrenscaleUp":"sf_facRenScaleUp",
        "facrenscaleDown":"sf_facRenScaleDown",
        "wtagUp:":"wTagScaleFactor_Tau21Up",
        "wtagDown":"wTagScaleFactor_Tau21Down",
        }

scaleWeights = ['facscaleUp','facscaleDown','renscaleUp','renscaleDown',
        'facrenscaleUp','facrenscaleDown']

def weight_mc(event, wHists, scale=1.0, weightOpts=[], errorOpt=None, debugLevel=0):
    """Apply pileup weights and other known MC correction factors"""
    lweightOpts = map(str.lower, weightOpts)

    #nominal weight
    eventWeight = event.weight*scale
    if debugLevel > 1: 
        print "Weight options:",weightOpts
        print "Error option:",errorOpt
        print "Weight from ntuple:",event.weight
        print "Scale by:",scale

    if len(lweightOpts) > 0:
        #top/W tag weighting
        if 'boost' in lweightOpts:
            eventWeight *= event.wTagScaleFactor * event.topTagScaleFactor

        #reweighting in number of b-jets
        if 'nbjets' in lweightOpts:
            eventWeight *= getNBJetsWeight(event, debugLevel=debugLevel)

        #top pt reweighting
        if 'toppt' in lweightOpts and hasattr(event, 'topPtWeight') and event.topPtWeight > 0:
            if debugLevel > 1:
                print "Top pt weight:",event.topPtWeight
            eventWeight *=  event.topPtWeight
        if 'nisr' in lweightOpts:
            nisrWeight = getNISRCorrection(event)
            if debugLevel > 1:
                print "NISR weight:",nisrWeight,event.NISRJets
            eventWeight *= nisrWeight

        #pileup reweighting
        if str.lower("reapplyNPUWeights") in lweightOpts:
            eventWeight *= reapplyPileupWeight(event, wHists, debugLevel=debugLevel)

        #lepton scale factors
        if str.lower("reapplyLepWeights") in lweightOpts:
            eventWeight *= reapplyLepEffWeight(event, wHists, debugLevel=debugLevel)
        if str.lower("reapplyTrigWeights") in lweightOpts:
            eventWeight *= reapplyLepTrigWeight(event, wHists, debugLevel=debugLevel)
        elif str.lower("removeTrigWeights") in lweightOpts:
            eventWeight /= event.trigWeight1L
            print event.trigWeight1L, event.lep1.Pt(), event.lep1.Eta()
        if str.lower("removePileupWeights") in lweightOpts:
            eventWeight /= event.pileupWeight
        if str.lower("removeBtagWeights") in lweightOpts:
            if hasattr(event, "btagW"):
                btagWeight = event.btagW
            elif hasattr(event, "btagCorrFactor"):
                btagWeight = event.btagCorrFactor
            else:
                raise RuntimeError("Warning: btag weight branch not found in tree!")
            if btagWeight > 0:
                eventWeight /= btagWeight

        #reweighting for ttbar dilepton control region
        if 'ttbardileptonmt' in lweightOpts:
            eventWeight *= getTTBarDileptonWeight(event)

    #up/down corrections for systematics
    normErrFraction=0.2
    if errorOpt is not None:
        if errorOpt in weightTable:
            errorWeight = getattr(event, weightTable[errorOpt])
            if not rt.TMath.IsNaN(errorWeight):
                eventWeight *= errorWeight
                if debugLevel > 1: 
                    print "%s scale factor: %.2f"%(errorOpt, errorWeight)
        elif 'normUp' in errorOpt:
            eventWeight *= (1+normErrFraction)
            if debugLevel > 1: print errorOpt,"scale factor:",1+normErrFraction
        elif 'normDown' in errorOpt:
            eventWeight /= (1+normErrFraction)
            if debugLevel > 1: print errorOpt,"scale factor:",1/(1+normErrFraction)

        # renormalize weights as needed
        if errorOpt in scaleWeights:
            if not ('nevents' in wHists and 'sumscaleweights' in wHists):
                print "Warning: missing weight normalization histograms"
            else:
                weightRescale = (wHists['nevents'].Integral() /
                        wHists['sumscaleweights'].GetBinContent(
                            scaleWeights.index(errorOpt)+1))
                if debugLevel > 1:
                    print "For",errorOpt,"rescaling weight by",weightRescale
                eventWeight *= weightRescale

    if debugLevel > 1: 
        print "event weight:",eventWeight
    return eventWeight

def weight_data(event, wHists, scale=1.0, weightOpts=[], errorOpt=None, debugLevel=0):
    lweightOpts = map(str.lower, weightOpts)

    if len(lweightOpts) == 0:
        eventWeight = scale
        if debugLevel > 1: print("Applying a weight of "+str(eventWeight))
        return eventWeight
    elif 'ttbardileptonmt' in lweightOpts:
        return scale * getTTBarDileptonWeight(event)
    else: 
        print("Warning: data weight options",lweightOpts,
                "may not make sense; see macro.razorWeights.weight_data")
        return eventWeight

def getMTRelUncertainty(MR, bkg, box):
    muonBoxes = ["MuMultiJet", "MuSixJet", "MuFourJet", "MuJet"]
    eleBoxes = ["EleMultiJet", "EleSixJet", "EleFourJet", "EleJet"]
    bkg = bkg.lower()
    unc = 0.0
    if bkg == "ttjets" or bkg == "ttbar":
        if box in muonBoxes:
            if MR < 400:
                unc = 0.0614822 #inclusive value -- no uncertainty is available for MR < 400
            elif MR < 600:
                unc = 0.0843497
            elif MR < 800:
                unc = 0.0977245
            elif MR < 1000:
                unc = 0.124865
            else:
                unc = 0.0950384
        if box in eleBoxes:
            if MR < 400:
                unc = 0.0337985
            elif MR < 600:
                unc = 0.0254958
            elif MR < 800:
                unc = 0.0338464
            elif MR < 1000:
                unc = 0.0376923
            else:
                unc = 0.0367005
    elif bkg == "wjets":
        if box in muonBoxes:
            if MR < 400:
                unc = 0.172134
            elif MR < 600:
                unc = 0.210269
            elif MR < 800:
                unc = 0.159052
            elif MR < 1000:
                unc = 0.242155
            else:
                unc = 0.263298
        if box in eleBoxes:
            if MR < 400:
                unc = 0.147676
            elif MR < 600:
                unc = 0.0656713
            elif MR < 800:
                unc = 0.0871765
            elif MR < 1000:
                unc = 0.0979534
            else:
                unc = 0.164399
    return unc

def applyMTUncertainty1D(hist, process, debugLevel=0):
    """hist is assumed to be a histogram of MR"""
    if process == "": return 

    bkg = process.split('_')[0].lower()
    box = process.split('_')[1]
    if bkg != "ttjets" and bkg != "ttbar" and bkg != "wjets": return

    #propagate MT uncertainty to each bin
    for bx in range(1,hist.GetNbinsX()+1):
        MRAtBinCenter = hist.GetXaxis().GetBinCenter(bx)
        mtRelUnc = getMTRelUncertainty(MRAtBinCenter, bkg, box)
        mtUnc = mtRelUnc*hist.GetBinContent(bx) #convert relative error to absolute
        if debugLevel > 0:
            print "MT uncertainty on bin",bx,"(",bkg,box,") is",mtUnc,("(%1.3f of %1.3f)" % (mtRelUnc,hist.GetBinContent(bx)))
        currUnc = hist.GetBinError(bx)
        hist.SetBinError(bx, (mtUnc*mtUnc + currUnc*currUnc)**(0.5))
        if debugLevel > 0:
            print "Uncertainty on this bin increases from",currUnc,"to",hist.GetBinError(bx)

def applyMTUncertainty2D(hist, process, debugLevel=0):
    """hist is assumed to be a histogram with MR on the x-axis"""
    if process == "": return 

    bkg = process.split('_')[0].lower()
    box = process.split('_')[1]
    if bkg != "ttjets" and bkg != "ttbar" and bkg != "wjets": return

    #propagate MT uncertainty to each bin
    for bx in range(1,hist.GetNbinsX()+1):
        MRAtBinCenter = hist.GetXaxis().GetBinCenter(bx)
        mtRelUnc = getMTRelUncertainty(MRAtBinCenter, bkg, box)
        for by in range(1,hist.GetNbinsY()+1):
            mtUnc = mtRelUnc*hist.GetBinContent(bx,by) #convert relative error to absolute
            if debugLevel > 0:
                print "MT uncertainty on bin",bx,",",by,"(",bkg,box,") is",mtUnc,("(%1.3f of %1.3f)" % (mtRelUnc,hist.GetBinContent(bx,by)))

            currUnc = hist.GetBinError(bx,by)
            hist.SetBinError(bx,by, (mtUnc*mtUnc + currUnc*currUnc)**(0.5))
            if debugLevel > 0:
                print "Uncertainty on this bin increases from",currUnc,"to",hist.GetBinError(bx,by)

def loadWeightHists(filenames={}, histnames={}, debugLevel=0):
    """Returns a dict with necessary reweighting histograms"""
    if "pileup" not in filenames:
        filenames["pileup"] = "data/ScaleFactors/Placeholders/DummyRun2PileupWeights.root"
        histnames["pileup"] = "PUWeight_Run2"
    if "muon" not in filenames:
        filenames["muon"] = "data/ScaleFactors/Placeholders/DummyRun2MuonWeights.root" #dummy file
        histnames["muon"] = "MuonWeight_Run2_Tight" 
    if "ele" not in filenames:
        filenames["ele"] = "data/ScaleFactors/Placeholders/DummyRun2EleWeights.root" #dummy file
        histnames["ele"] = "EleWeight_Run2_Tight" #razor tag&probe
    wHists = {}
    wFiles = {}
    for name in filenames:
        if debugLevel > 0: print "Opening",name,"weight file",filenames[name]
        wFiles[name] = rt.TFile.Open(filenames[name])
        if debugLevel > 0: print "Getting histogram",histnames[name]
        wHists[name] = wFiles[name].Get(histnames[name])
        wHists[name].SetDirectory(0)
        assert wHists[name]
    for name in wFiles:
        wFiles[name].Close()
    return wHists

def getSFHistNameForErrorOpt(errorOpt, name):
    """Returns the key in the scale factor histogram dictionary for the given error option"""
    if errorOpt in ['sfsysttjetsUp','sfsyswjetsUp','sfsyszinvUp']:
        return name+'Up'
    elif errorOpt in ['sfsysttjetsDown','sfsyswjetsDown','sfsyszinvDown']:
        return name+'Down'
    else:
        return name

vetoLeptonAuxCuts="(abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13) && leadingGenLeptonPt > 5"
vetoTauAuxCuts="abs(leadingGenLeptonType) == 15 && leadingGenLeptonPt > 20"
def getAuxSFsForErrorOpt(auxSFs={}, errorOpt="", auxSFsPerProcess=False):
    """
    Returns scale factor histogram names needed to compute the indicated shape uncertainty.
    Format of the input/output is { "HistogramName":("variableName", "cuts"), ... }
    """

    #for building output
    histNames=[]
    varNames=[]
    cuts=[]

    ###supported error options
    #TTJets dilepton control region systematic
    if 'ttcrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("TTJetsDileptonUp")
        elif 'Down' in errorOpt:
            histNames.append("TTJetsDileptonDown")
        varNames.append("MR")
        cuts.append("1")
    #DYJets dilepton control region systematic
    elif 'zllcrosscheckmr' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("DYJetsInvMRUp")
        elif 'Down' in errorOpt:
            histNames.append("DYJetsInvMRDown")
        varNames.append("MR")
        cuts.append("1")
    elif 'zllcrosscheckrsq' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("DYJetsInvRsqUp")
        elif 'Down' in errorOpt:
            histNames.append("DYJetsInvRsqDown")
        varNames.append("Rsq")
        cuts.append("1")
    #Veto lepton scale factors up/down
    elif 'vetolepptcrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("VetoLeptonPtUp")
        elif 'Down' in errorOpt:
            histNames.append("VetoLeptonPtDown")
        varNames.append("leadingGenLeptonPt")
        cuts.append(vetoLeptonAuxCuts)
    elif 'vetolepetacrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("VetoLeptonEtaUp")
        elif 'Down' in errorOpt:
            histNames.append("VetoLeptonEtaDown")
        varNames.append("abs(leadingGenLeptonEta)")
        cuts.append(vetoLeptonAuxCuts)
    #Veto tau scale factors up/down
    elif 'vetotauptcrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("VetoTauPtUp")
        elif 'Down' in errorOpt:
            histNames.append("VetoTauPtDown")
        varNames.append("leadingGenLeptonPt")
        cuts.append(vetoTauAuxCuts)
    elif 'vetotauetacrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("VetoTauEtaUp")
        elif 'Down' in errorOpt:
            histNames.append("VetoTauEtaDown")
        varNames.append("abs(leadingGenLeptonEta)")
        cuts.append(vetoTauAuxCuts)
    #b-tag bins closure test systematic
    elif 'btagcrosscheckmr' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("MRBUp")
        elif 'Down' in errorOpt:
            histNames.append("MRBDown")
        varNames.append("MR")
        cuts.append("1")
    elif 'btagcrosscheckrsq' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("RsqBUp")
        elif 'Down' in errorOpt:
            histNames.append("RsqBDown")
        varNames.append("Rsq")
        cuts.append("1")
    #b-tag closure test systematic for ZInv
    elif 'btaginvcrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append('ZInvBUp')
        elif 'Down' in errorOpt:
            histNames.append('ZInvBDown')
        varNames.append('Rsq')
        cuts.append('1')

    #return dictionary with needed information
    sfsNeeded = { histNames[i]:(varNames[i],cuts[i]) for i in range(len(histNames)) }
    if auxSFsPerProcess:
        for process in auxSFs:
            auxSFs[process].update(sfsNeeded.copy())
    else:
        auxSFs.update(sfsNeeded)

def splitShapeErrorsByType(shapeErrors):
    """Takes a list of shape uncertainties and splits it into two lists: the first is the list of uncertainties applied as per-event scale factors, and the second is the list of uncertainties that require separate processing."""
    supportedShapeUncertainties = { #True: belongs in list 1.  False: belongs in list 2
        'jes':False,
        'ees':False,
        'mes':False,
        'npvextrap':False,
        'genmetvspfmet':False,
        'btag':True,
        'pileup':True,
        'bmistag':True,
        'facscale':True,
        'renscale':True,
        'facrenscale':True,
        'ttcrosscheck':True,
        'zllcrosscheckmr':True,
        'zllcrosscheckrsq':True,
        'btagcrosscheckmr':True,
        'btagcrosscheckrsq':True,
        'btaginvcrosscheck':True,
        'vetolepptcrosscheck':True,
        'vetotauptcrosscheck':True,
        'vetolepetacrosscheck':True,
        'vetotauetacrosscheck':True,
        'singletopnorm':True,
        'othernorm':True,
        'qcdnorm':True,
        'qcdbtag':True,
        'qcdaltfunction':True,
        'sfsysttjets':True,
        'sfsyswjets':True,
        'sfsyszinv':True,
        'sfstatttjets':True,
        'sfstatwjets':True,
        'sfstatzinv':True,
        'vetomuoneff':True,
        'vetoeleeff':True,
        'tightmuoneff':True,
        'tighteleeff':True,
        'tightmuonfastsim':True,
        'tightelefastsim':True,
        'vetomuonfastsim':True,
        'vetoelefastsim':True,
        'btagfastsim':True,
        'isr':True,
        'muontrig':True,
        'eletrig':True,
        'wtag':True,
        }
    sfUncertainties = []
    otherUncertainties = []
    for shape in shapeErrors:
        #if shape is wrapped in a tuple, unwrap it
        if not isinstance(shape, basestring):
            thisShape = shape[0]
        else:
            thisShape = shape
        if thisShape in supportedShapeUncertainties:
            if supportedShapeUncertainties[thisShape]:
                sfUncertainties.append(shape)
            else:
                otherUncertainties.append(shape)
        else:
            print "Warning in splitShapeErrorsByType: error option",thisShape,"is not supported!"

    return sfUncertainties, otherUncertainties

def getNJetsSFs(analysis,jetName='NJets40'):
    """From an Analysis object, get the needed NJets scale factors"""
    auxSFs = { process:{} for process in analysis.samples }
    if "WJets" in analysis.samples:
        auxSFs["WJets"] = {'NJetsWJets':(jetName,'1')}
    for name in ['TTJets','TTJets1L','TTJets2L']:
        if name in analysis.samples:
            auxSFs[name] = {'NJetsTTJets':(jetName,'1')}
    for name in ['WJetsInv','DYJetsInv','ZInv','GJetsInv']:
        if name in analysis.samples:
            auxSFs[name] = {"NJetsInv":(jetName,"1")}
    return auxSFs

def addBTagSFs(analysis, auxSFs={}, var='MR', gjets=False):
    procs = ['WJets', 'TTJets', 'TTJets1L', 'TTJets2L']
    sfKey = 'MR'
    if gjets:
        procs = ['WJetsInv','DYJetsInv','ZInv','GJetsInv']
        sfKey = 'MRInv'
    for proc in procs:
        if proc in analysis.samples:
            if proc not in auxSFs:
                auxSFs[proc] = {}
            auxSFs[proc][sfKey] = (var, '1')
    return auxSFs

def addBTagDoubleRatioSFs(analysis, auxSFs={}, var='nBTaggedJets'):
    proc = 'ZInv'
    if proc in analysis.samples:
        if proc not in auxSFs:
            auxSFs[proc] = {}
        auxSFs[proc]['NBTagsInv'] = (var, '1')
    return auxSFs

def addAllBTagSFs(analysis, auxSFs={}, var='MR', gjets=False, 
        bjetsName='NBJetsMedium'):
    """
    Gets scale factor directives for each b-tag bin separately.
    Use to process a sample inclusive in the number of btags.
    """
    procs = ['WJets', 'TTJets', 'TTJets1L', 'TTJets2L']
    sfKey = 'MR'
    if gjets:
        procs = ['WJetsInv','DYJetsInv','ZInv','GJetsInv']
        sfKey = 'MRInv'
    for proc in procs:
        if proc in analysis.samples:
            if proc not in auxSFs:
                auxSFs[proc] = {}
            nbMax = 2
            if (analysis.njetsMin >= 4 or 'MultiJet' in analysis.region or 'SevenJet' in analysis.region) and not gjets:
                nbMax = 3
            for nb in range(nbMax+1):
                rel = '=='
                if nb == nbMax:
                    rel = '>='
                auxSFs[proc]['{}{}B'.format(sfKey, nb)] = (var, 
                        '{} {} {}'.format(bjetsName, rel, nb))
    return auxSFs

def loadPhotonPurityHists(sfHists={}, tag="Razor2016", debugLevel=0):
    filenames = { 
            "Razor2016_MoriondRereco":"data/ScaleFactors/RazorMADD2015/PhotonCR_Purity_2016_Rereco_36p8ifb.root",
            "Razor2016G_SUSYUnblind_80X":"data/ScaleFactors/RazorMADD2015/PhotonCR_Purity_2016G_SUSYUnblind.root",
            "Razor2016_ICHEP_80X":"/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/ChargedIsoFits/12p9/PhotonControlRegionPlots_MR300Rsq0p15.root"
            }
    if tag not in filenames:
        sys.exit("tag %s is not supported for photon purity measurement"%tag)
    filename = filenames[tag]
    infile = rt.TFile.Open(filename)
    sfHists["PhotonPurityEB"] = infile.Get("histChargedIso_EB_MRRsq")
    sfHists["PhotonPurityEB"].SetDirectory(0)
    sfHists["PhotonPurityEE"] = infile.Get("histChargedIso_EE_MRRsq")
    sfHists["PhotonPurityEE"].SetDirectory(0)
    # Instead of purity values, switch to 1-purity
    for ibin in range(sfHists["PhotonPurityEB"].GetSize()+1):
        sfHists["PhotonPurityEB"].SetBinContent( ibin, 1 - sfHists["PhotonPurityEB"].GetBinContent(ibin) )
        sfHists["PhotonPurityEE"].SetBinContent( ibin, 1 - sfHists["PhotonPurityEE"].GetBinContent(ibin) )
    return sfHists

def getPhotonPuritySFs(auxSFs={}):
    if "QCD" not in auxSFs:
        auxSFs["QCD"] = {}
    auxSFs["QCD"]['PhotonPurityEB'] = (('MR_NoPho','Rsq_NoPho'),'abs(pho1.Eta()) < 1.479')
    auxSFs["QCD"]['PhotonPurityEE'] = (('MR_NoPho','Rsq_NoPho'),'abs(pho1.Eta()) >= 1.479')
    return auxSFs

def getNPVHist(tag):
    """Loads the NPV histogram from the file corresponding to the given
        analysis tag, and returns it."""
    filenames = {
            "Razor2016_MoriondRereco":
                "data/PileupWeights/NPVTarget_2016_35p9ifb.root"
                }   
    if tag not in filenames:
        sys.exit("tag %s is not supported for NPV reweighting"%(tag))
    inFile = rt.TFile.Open(filenames[tag])
    npvHist = inFile.Get("NPV_2016")
    assert npvHist
    npvHist.SetDirectory(0)
    return npvHist

def makeWeightHistDict(fileDict, debugLevel=0):
    """Gets NEvents and SumScaleWeights hists from each file"""
    hists = {}
    for name,f in fileDict.iteritems():
        if debugLevel > 0: 
            print "Loading weight histograms for process",name
        nevents = f.Get("NEvents")
        if nevents:
            nevents.SetDirectory(0)
        sumScaleWeights = f.Get("SumScaleWeights")
        if sumScaleWeights:
            sumScaleWeights.SetDirectory(0)
        hists[name] = {'nevents':nevents, 
                'sumscaleweights':sumScaleWeights}
    return hists

