import array
import ROOT as rt
rt.gROOT.SetBatch()
rt.gStyle.SetPaintTextFormat(".3f")

def get_cuts(box=None, dphiregion=None, mc=False, sideband=False):
    cuts = [
        "MR > 550 && MR < 4000",
        "(Rsq > 0.2 || (MR > 1600 && Rsq > 0.1))",
        "(box==11||box==12||box==14)",
        "(HLTDecision[163] || HLTDecision[166] || HLTDecision[167])",
        "Flag_HBHENoiseFilter", "Flag_HBHEIsoNoiseFilter", "Flag_goodVertices",
        "Flag_eeBadScFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", 
        "Flag_CSCTightHaloFilter", "Flag_badChargedCandidateFilter", 
        "Flag_badMuonFilter" ]
    if sideband:
        cuts.append("Rsq > 0.3 && MR < 650")
    else:
        cuts.append("Rsq < 0.3")
    box_cuts = { 'dijet'   :'nSelectedJets >= 2 && nSelectedJets <= 3',
                 'multijet':'nSelectedJets >= 4 && nSelectedJets <= 6',
                 'sevenjet':'nSelectedJets >= 7'
                 }
    dphi_cuts = { 'lo':'abs(dPhiRazor) < 2.8',
                  'hi':'abs(dPhiRazor) > 2.8', }
    if box is not None:
        cuts.append(box_cuts[box])
    if dphiregion is not None:
        cuts.append(dphi_cuts[dphiregion])
    prefix = '('
    suffix = ')'
    if mc:
        prefix = '35800*mcScaleFactor*weight*'+prefix
    return prefix+' && '.join(cuts)+suffix

def get_binning(box, region):
    if region == 'lo':
        if box == 'sevenjet':
            return { 'MR':[650, 800, 1600], 
                    'Rsq':[0.20, 0.22, 0.24, 0.26, 0.28, 0.30] }
        return { 'MR': [650, 800, 1000, 1200, 1400, 1600],
                'Rsq':[0.20, 0.22, 0.24, 0.26, 0.28, 0.30] }
    elif region == 'mrsideband':
        return { 'MR': [550, 650],
                'Rsq':[0.20, 0.22, 0.24, 0.26, 0.28, 
                    0.30, 0.41, 0.52, 0.64, 1.5] }
    else:
        if box == 'dijet':
            return { 'MR': [1600, 4000],
                    'Rsq':[0.10, 0.12, 0.14, 0.16, 0.18, 
                        0.20, 0.22, 0.24, 0.26, 0.28, 0.30] }
        else:
            return { 'MR': [1600, 4000],
                    'Rsq':[0.10, 0.12, 0.14, 0.16, 0.18, 0.20] }

def get_hist_name(box, region, dphiregion, mc=False):
    name = 'h{}mr{}dphi{}'.format(box, region, dphiregion)
    if mc:
        return name+'mc'
    else:
        return name

def make_hist(name, binning):
    hist = rt.TH2F(name, 'QCD SFs',
            len(binning['MR'])-1, array.array('d', binning['MR']),
            len(binning['Rsq'])-1, array.array('d', binning['Rsq'])
            )
    hist.Sumw2()
    return hist

def get_tf_err(hist_lo, hist_hi, mc_hist_lo, mc_hist_hi,
        ix, iy):
    numer = hist_lo.GetBinContent(ix, iy) - mc_hist_lo.GetBinContent(ix, iy)
    denom = hist_hi.GetBinContent(ix, iy) - mc_hist_hi.GetBinContent(ix, iy)
    err_numer = (hist_lo.GetBinError(ix, iy)**2 + mc_hist_lo.GetBinError(ix, iy)**2)**(0.5)
    err_denom = (hist_hi.GetBinError(ix, iy)**2 + mc_hist_hi.GetBinError(ix, iy)**2)**(0.5)
    # inflate uncertainty by this amount to cover the assumption 
    # that the scale factors do not depend on the number of b-tags
    inflate_factor = 1.3
    return inflate_factor * ((err_numer/denom)**2 
            + (numer*err_denom/(denom*denom))**2)**(0.5)

def set_tf(hist_out, hist_lo, hist_hi, mc_hist_lo, mc_hist_hi, ix, iy):
    numer = max(0, hist_lo.GetBinContent(ix, iy) - mc_hist_lo.GetBinContent(ix, iy))
    denom = max(0, hist_hi.GetBinContent(ix, iy) - mc_hist_hi.GetBinContent(ix, iy))
    if denom == 0:
        print "Denominator is zero in bin {} {}".format(ix, iy)
        hist_out.SetBinContent(ix, iy, 0)
        hist_out.SetBinError(ix, iy, 0) # this will exclude it from the fit
        return
    tf = numer / denom
    tf_err = get_tf_err(hist_lo, hist_hi, mc_hist_lo, mc_hist_hi, ix, iy)
    hist_out.SetBinContent(ix, iy, tf)
    hist_out.SetBinError(ix, iy, tf_err)


if __name__ == '__main__':
    boxes = ['dijet', 'multijet', 'sevenjet']
    regions = ['lo', 'hi', 'mrsideband']
    dphiregions = ['lo', 'hi']
    draw_string = 'Rsq:MR>>'
    mc_draw_string = draw_string+'+'
    
    fD = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_27Nov2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_Data_NoDuplicates_GoodLumiGolden.root","read")
    tree = fD.Get("RazorInclusive")
    bkg_files = {
            'ttjets':rt.TFile.Open('TTJets.root'),
            'wjets':rt.TFile.Open('WJets.root'),
            'zinv':rt.TFile.Open('ZInv.root'),
            }
    bkg_trees = { proc:f.Get("RazorInclusive") for proc, f in bkg_files.iteritems() } 

    out_name = 'qcd_data.root'
    out_f = rt.TFile.Open(out_name, 'recreate')

    for box in boxes:
        for region in regions:
            binning = get_binning(box, region)
            hist_out = make_hist('{}_{}'.format(box, region), binning)
            hists = []
            mc_hists = []
            for dphiregion in dphiregions:
                cuts = get_cuts(box, dphiregion, 
                        sideband=(region=='mrsideband'))
                hist_name = get_hist_name(box, region, dphiregion)
                print "Doing {}".format(hist_name)
                hists.append(make_hist(hist_name, binning))
                print "Cuts",cuts
                nevents = tree.Draw(draw_string+hist_name, cuts)
                mc_cuts = get_cuts(box, dphiregion, mc=True,
                        sideband=(region=='mrsideband'))
                print "MC cuts",mc_cuts
                mc_hist_name = get_hist_name(box, region, dphiregion, mc=True)
                mc_hists.append(make_hist(mc_hist_name, binning))
                for bkg, t in bkg_trees.iteritems():
                    nevents = t.Draw(mc_draw_string+mc_hist_name, mc_cuts)
            for ix in range(1, hist_out.GetNbinsX()+1):
                for iy in range(1, hist_out.GetNbinsY()+1):
                    set_tf(hist_out, hists[0], hists[1],
                        mc_hists[0], mc_hists[1], ix, iy)
            hist_out.Write()
    out_f.Close()
