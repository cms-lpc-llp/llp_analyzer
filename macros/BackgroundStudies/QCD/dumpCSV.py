### This creates the CSV file read as input by "QCD Transfer Factors.ipynb"

import array
import csv
import ROOT as rt
rt.gROOT.SetBatch()
rt.gStyle.SetPaintTextFormat(".3f")

def get_cuts(box, btags, dphiregion, mc=False):
    cuts = [
        "MR > 400 && MR < 4000",
        "(Rsq > 0.2 || (MR > 1600 && Rsq > 0.1))",
        "(box==11||box==12||box==14)",
        "Flag_HBHENoiseFilter", 
        "Flag_HBHEIsoNoiseFilter", 
        "Flag_goodVertices",
        "Flag_eeBadScFilter", 
        "Flag_EcalDeadCellTriggerPrimitiveFilter", 
        "Flag_CSCTightHaloFilter",
        "Flag_badChargedCandidateFilter", 
        "Flag_badMuonFilter"
        ]
    box_cuts = {
            'dijet'   :'nSelectedJets >= 2 && nSelectedJets <= 3',
            'multijet':'nSelectedJets >= 4 && nSelectedJets <= 6',
            'sevenjet':'nSelectedJets >= 7'
            }
    btag_cuts = {
            0:'nBTaggedJets == 0',
            1:'nBTaggedJets == 1',
            2:'nBTaggedJets >= 2',
            }
    dphi_cuts = {
            'lo':'abs(dPhiRazor) < 2.8',
            'hi':'abs(dPhiRazor) > 2.8',
            }
    cuts.append(box_cuts[box])
    cuts.append(btag_cuts[btags])
    cuts.append(dphi_cuts[dphiregion])
    prefix = '('
    suffix = ')'
    if mc:
        prefix = '35800*mcScaleFactor*weight*'+prefix
    return prefix+' && '.join(cuts)+suffix

def get_binning(box, region):
    if box != 'sevenjet':
        if region == 'lo':
            return {
                    'MR': [400, 600, 800, 1000, 1200, 1400, 1600],
                    'Rsq':[0.20, 0.22, 0.24, 0.26, 0.28, 0.30]
                    }
        else:
            return {
                    'MR': [1600, 4000],
                    'Rsq':[0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
                    }
    else:
        if region == 'lo':
            return {
                    'MR': [400, 600, 800, 1000, 1600],
                    'Rsq':[0.2, 0.225, 0.25, 0.275, 0.3]
                    }
        else:
            return {
                    'MR': [1600, 4000],
                    'Rsq':[0.10, 0.13, 0.20]
                    }

def get_hist_name(box, btags, region, dphiregion, mc=False):
    name = 'h{}{}bmr{}dphi{}'.format(box, btags, region, dphiregion)
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
    return ((err_numer/denom)**2 + (numer*err_denom/(denom*denom))**2)**(0.5)

def get_line(hist_lo, hist_hi, mc_hist_lo, mc_hist_hi,
        ix, iy, box, nb, region):
    line = []
    line.append(hist_lo.GetXaxis().GetBinLowEdge(ix))
    line.append(hist_lo.GetXaxis().GetBinUpEdge(ix))
    line.append(hist_lo.GetYaxis().GetBinLowEdge(iy))
    line.append(hist_lo.GetYaxis().GetBinUpEdge(iy))
    line.append(box)
    line.append(nb)
    line.append(region)
    line.append(hist_lo.GetBinContent(ix, iy))
    line.append(mc_hist_lo.GetBinContent(ix, iy))
    line.append(mc_hist_lo.GetBinError(ix, iy))
    line.append(hist_hi.GetBinContent(ix, iy))
    line.append(mc_hist_hi.GetBinContent(ix, iy))
    line.append(mc_hist_hi.GetBinError(ix, iy))
    line.append((hist_lo.GetBinContent(ix, iy) - mc_hist_lo.GetBinContent(ix, iy))
            / (hist_hi.GetBinContent(ix, iy) - mc_hist_hi.GetBinContent(ix, iy)))
    line.append(get_tf_err(hist_lo, hist_hi, mc_hist_lo, mc_hist_hi, ix, iy))
    return line


if __name__ == '__main__':
    boxes = ['dijet', 'multijet', 'sevenjet']
    btags = [0, 1, 2]
    regions = ['lo', 'hi']
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

    header = ['mr_lo', 'mr_hi', 'rsq_lo', 'rsq_hi', 'box', 'nb', 
            'region', 'counts_lo', 'bkg_lo', 'bkg_err_lo', 
            'counts_hi', 'bkg_hi', 'bkg_err_hi', 'transfer', 'transfer_err']
    out_name = 'qcd_data.csv'

    with open(out_name, 'wb') as out_f:
        writer = csv.writer(out_f)
        writer.writerow(header)
        for box in boxes:
            for nb in btags:
                for lohi in regions:
                    binning = get_binning(box, lohi)
                    hists = []
                    mc_hists = []
                    for dphiregion in dphiregions:
                        cuts = get_cuts(box, nb, dphiregion)
                        hist_name = get_hist_name(box, nb, lohi, dphiregion)
                        print "Doing {}".format(hist_name)
                        hists.append(make_hist(hist_name, binning))
                        tree.Draw(draw_string+hist_name, cuts)
                        mc_cuts = get_cuts(box, nb, dphiregion, mc=True)
                        mc_hist_name = get_hist_name(box, nb, lohi, dphiregion, mc=True)
                        mc_hists.append(make_hist(mc_hist_name, binning))
                        for bkg, t in bkg_trees.iteritems():
                            t.Draw(mc_draw_string+mc_hist_name, mc_cuts)
                    ref_hist = hists[0]
                    for ix in range(1, ref_hist.GetNbinsX()+1):
                        for iy in range(1, ref_hist.GetNbinsY()+1):
                            writer.writerow(get_line(hists[0], hists[1],
                                mc_hists[0], mc_hists[1], ix, iy, box, nb, lohi))
