import array
import ROOT as rt

from getInclusiveTransferFactors import get_cuts, get_binning
from mcStudy import normalize_columns

def make_hist(name, binning):
    hist = rt.TH2F(name, 'QCD',
            len(binning['Rsq'])-1, array.array('d', binning['Rsq']),
            len(binning['nBTaggedJets'])-1, array.array('d', binning['nBTaggedJets'])
            )
    hist.Sumw2()
    return hist

def print_hist(hist):
    name = hist.GetName()
    c = rt.TCanvas('c'+name, 'c', 400, 300)
    hist.SetStats(0)
    hist.GetXaxis().SetTitle('R^{2}')
    hist.GetYaxis().SetTitle('N_{b-tags}')
    hist.SetTitle('')
    hist.Draw("colztexte")
    c.Print('QCD_NBtags_{}.pdf'.format(name))

def write_sys_hist(hist, f):
    """
    Computes a systematic uncertainty based on the spread
    of the data in Rsq
    """
    syshist = rt.TH1F(hist.GetName().replace('lo', 'btags')
            +'sys', '',
            hist.GetNbinsY(), 0, hist.GetNbinsY()-1)
    syshist.SetDirectory(0)
    for btagbin in range(1, hist.GetNbinsY()+1):
        minq = -1
        maxq = -1
        for rsqbin in range(1, hist.GetNbinsX()+1):
            thisq = hist.GetBinContent(rsqbin, btagbin)
            if minq == -1 or thisq < minq:
                minq = thisq
            if maxq == -1 or thisq > maxq:
                maxq = thisq
        syshist.SetBinContent(btagbin, maxq-minq)
    f.cd()
    syshist.Write()

def write_hist(hist, f):
    f.cd()
    name = hist.GetName().replace('lo', 'btags')
    out_hist = hist.Clone(hist.GetName()+'_out')
    out_hist = out_hist.ProjectionY(name)
    out_hist.Scale(1./out_hist.Integral())
    out_hist.Write()


if __name__ == '__main__':
    boxes = ['dijet', 'multijet', 'sevenjet']
    regions = ['lo', 'hi']
    draw_string = 'nBTaggedJets:Rsq>>'
    mc_draw_string = draw_string+'+'
    
    fD = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_27Nov2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root","read")
    tree = fD.Get("RazorInclusive")
    bkg_files = {
            'ttjets':rt.TFile.Open('TTJets.root'),
            'wjets':rt.TFile.Open('WJets.root'),
            'zinv':rt.TFile.Open('ZInv.root'),
            }
    bkg_trees = { proc:f.Get("RazorInclusive") for proc, f in bkg_files.iteritems() } 

    out_f = rt.TFile.Open("qcd_data_btags.root", 'recreate')
    for box in boxes:
        btags = [0, 1, 2, 3, 4]
        if box == 'dijet':
            btags = [0, 1, 2, 3]
        for region in regions:
            print "Doing {} {}".format(box, region)
            binning = get_binning(box, region)
            # I am removing one bin division in Rsq
            # to increase statistics
            try:
                binning['Rsq'].remove(0.28)
            except ValueError:
                pass
            binning['nBTaggedJets'] = btags
            hist_name = box+region
            hist = make_hist(hist_name, binning)
            mc_hist_name = 'mc'+hist_name
            mc_hist = make_hist(mc_hist_name, binning)
            cuts = get_cuts(box)
            tree.Draw(draw_string+hist_name, cuts)
            mc_cuts = get_cuts(box, mc=True)
            for bkg, t in bkg_trees.iteritems():
                t.Draw(mc_draw_string+mc_hist_name, mc_cuts)
            hist.Add(mc_hist, -1)
            if region == 'lo':
                write_hist(hist, out_f)
            normalize_columns(hist)
            if region == 'lo':
                write_sys_hist(hist, out_f)
            print_hist(hist)
    out_f.Close()
