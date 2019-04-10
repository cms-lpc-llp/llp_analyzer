import array
import ROOT as rt

from getInclusiveTransferFactors import get_cuts, get_binning
from mcStudy import normalize_columns

def make_hist(name, binning):
    hist = rt.TH2F(name, 'QCD',
            len(binning['dphi'])-1, array.array('d', binning['dphi']),
            len(binning['nSelectedJets'])-1, array.array('d', binning['nSelectedJets'])
            )
    hist.Sumw2()
    return hist

def get_cuts_var(var, var_min, var_max, mc=False):
    cuts = get_cuts(mc=mc)
    var_cut_str = " && {} >= {} && {} < {}".format(var, var_min, var, var_max)
    # chop off end paren and add additional cuts
    cuts = cuts[:-1] + var_cut_str + ')'
    return cuts

def print_hist(hist):
    name = hist.GetName()
    c = rt.TCanvas('c'+name, 'c', 400, 300)
    hist.SetStats(0)
    hist.GetXaxis().SetTitle('#Delta #phi_{razor}')
    hist.GetYaxis().SetTitle('N_{jets}')
    hist.SetTitle('')
    hist.Draw("colztexte")
    c.Print('QCD_NJetsVsDPhi_{}.pdf'.format(name))

def write_hist(hist, f):
    f.cd()
    name = hist.GetName().replace('lo', 'btags')
    out_hist = hist.Clone(hist.GetName()+'_out')
    out_hist = out_hist.ProjectionY(name)
    out_hist.Scale(1./out_hist.Integral())
    out_hist.Write()


if __name__ == '__main__':
    regions = ['lo', 'hi']
    draw_string = 'nSelectedJets:abs(dPhiRazor)>>'
    mc_draw_string = draw_string+'+'
    
    fD = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_05Oct2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_Data_NoDuplicates_GoodLumiGolden.root","read")
    tree = fD.Get("RazorInclusive")
    bkg_files = {
            'ttjets':rt.TFile.Open('TTJets.root'),
            'wjets':rt.TFile.Open('WJets.root'),
            'zinv':rt.TFile.Open('ZInv.root'),
            }
    bkg_trees = { proc:f.Get("RazorInclusive") for proc, f in bkg_files.iteritems() } 

    out_f = rt.TFile.Open("qcd_data_njets.root", 'recreate')
    njets = [2, 4, 7, 20]
    dphi = [0.0, 2.8, 3.14159]
    btags = [0, 1, 2, 3, 4]
    for region in regions:
        binning = get_binning('multijet', region=region)
        binning['nSelectedJets'] = njets
        binning['dphi'] = dphi
        binning['nBTaggedJets'] = btags
        for var in ['MR', 'Rsq', 'nBTaggedJets']:
            for var_min, var_max in zip(binning[var][:-1], 
                    binning[var][1:]):
                print "Doing {} {} {}-{}".format(region, var, var_min, var_max)
                hist_name = region+var+str(var_min)
                hist = make_hist(hist_name, binning)
                mc_hist_name = 'mc'+hist_name
                mc_hist = make_hist(mc_hist_name, binning)
                cuts = get_cuts_var(var, var_min, var_max)
                tree.Draw(draw_string+hist_name, cuts)
                mc_cuts = get_cuts_var(var, var_min, var_max, mc=True)
                for bkg, t in bkg_trees.iteritems():
                    t.Draw(mc_draw_string+mc_hist_name, mc_cuts)
                hist.Add(mc_hist, -1)
                if region == 'lo':
                    write_hist(hist, out_f)
                normalize_columns(hist)
                print_hist(hist)
    out_f.Close()
