import argparse
import array
import ROOT as rt
rt.gROOT.SetBatch()
rt.gStyle.SetPaintTextFormat(".3f")

from getInclusiveTransferFactors import set_tf, get_tf_err, make_hist

def get_cuts(njets, dphiregion=None, mc=False,
        exclusive=False, sideband=False):
    cuts = [
        "MR > 550 && MR < 4000",
        "(Rsq > 0.2 || (MR > 1600 && Rsq > 0.1))",
        "(box==11||box==12||box==14)",
        "Flag_HBHENoiseFilter", "Flag_HBHEIsoNoiseFilter", "Flag_goodVertices",
        "Flag_eeBadScFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", 
        "Flag_CSCTightHaloFilter", "Flag_badChargedCandidateFilter", 
        "Flag_badMuonFilter" ]
    if sideband:
        cuts.append("Rsq > 0.3 && MR < 650")
    else:
        cuts.append("Rsq < 0.3")
    rel = ">="
    if exclusive:
        rel = "=="
    box_cut = "nSelectedJets {} {}".format(rel, njets)
    dphi_cuts = { 'lo':'abs(dPhiRazor) < 2.8',
                  'hi':'abs(dPhiRazor) > 2.8', }
    cuts.append(box_cut)
    if dphiregion is not None:
        cuts.append(dphi_cuts[dphiregion])
    prefix = '('
    suffix = ')'
    if mc:
        prefix = '35800*mcScaleFactor*weight*'+prefix
    return prefix+' && '.join(cuts)+suffix

def get_binning(region=''):
    if region == 'mrsideband':
        return { 'MR':[550, 650],
                'Rsq':[0.30, 0.41, 0.52, 0.64, 1.5] }
    elif region == 'mrhi': 
        return { 'MR':[1600, 4000],
                'Rsq':[0.10, 0.12, 0.14, 0.16, 0.18, 0.20] }
    return { 'MR': [550, 650, 800, 1000, 1200, 1400, 1600],
            'Rsq':[0.20, 0.22, 0.24, 0.26, 0.28, 0.30] }

def get_hist_name(njets, dphiregion, mc=False, region=''):
    name = 'h{}jetsdphi{}'.format(njets, dphiregion)
    if region == 'mrsideband':
        name += 'sideband'
    elif region == 'mrhi':
        name += 'mrhi'
    if mc:
        return name+'mc'
    else:
        return name


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exclusive', action='store_true',
            help='Make njets bins exclusive')
    args = parser.parse_args()

    njets = [2, 3, 4, 5, 6, 7]
    regions = ['', 'mrhi', 'mrsideband']
    dphiregions = ['lo', 'hi']
    draw_string = 'Rsq:MR>>'
    mc_draw_string = draw_string+'+'
    
    fD = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_05Oct2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_Data_NoDuplicates_GoodLumiGolden.root","read")
    tree = fD.Get("RazorInclusive")
    bkg_files = {
            'ttjets':rt.TFile.Open('TTJets.root'),
            'wjets':rt.TFile.Open('WJets.root'),
            'zinv':rt.TFile.Open('ZInv.root'),
            }
    bkg_trees = { proc:f.Get("RazorInclusive") for proc, f in bkg_files.iteritems() } 

    tag = 'inclusive'
    if args.exclusive:
        tag = 'exclusive'
    out_name = 'qcd_njets_data_{}.root'.format(tag)
    out_f = rt.TFile.Open(out_name, 'recreate')

    for region in regions:
        sideband = (region == 'mrsideband')
        binning = get_binning(region)
        for nj in njets:
            hist_out = make_hist('nj{}{}'.format(nj, region), binning)
            hists = []
            mc_hists = []
            for dphiregion in dphiregions:
                cuts = get_cuts(nj, dphiregion, 
                        exclusive=args.exclusive, sideband=sideband)
                hist_name = get_hist_name(nj, dphiregion, region=region)
                print "Doing {}".format(hist_name)
                hists.append(make_hist(hist_name, binning))
                tree.Draw(draw_string+hist_name, cuts)
                mc_cuts = get_cuts(nj, dphiregion, mc=True,
                        exclusive=args.exclusive, sideband=sideband)
                mc_hist_name = get_hist_name(nj, dphiregion, mc=True,
                        region=region)
                mc_hists.append(make_hist(mc_hist_name, binning))
                for bkg, t in bkg_trees.iteritems():
                    t.Draw(mc_draw_string+mc_hist_name, mc_cuts)
            for ix in range(1, hist_out.GetNbinsX()+1):
                for iy in range(1, hist_out.GetNbinsY()+1):
                    set_tf(hist_out, hists[0], hists[1],
                        mc_hists[0], mc_hists[1], ix, iy)
            hist_out.Write()
    out_f.Close()
