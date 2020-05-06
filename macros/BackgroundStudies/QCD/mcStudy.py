import array
import ROOT as rt
rt.gROOT.SetBatch()
rt.gStyle.SetPaintTextFormat(".3f")

def normalize_columns(hist):
    """
    Input: TH2 
    Output: the same TH2 with each column normalized to 1
    """
    for ix in range(1, hist.GetNbinsX()+1):
        col_total = 0.
        col_total_err = 0.
        for iy in range(1, hist.GetNbinsY()+1):
            col_total += hist.GetBinContent(ix, iy)
            col_total_err += hist.GetBinError(ix, iy)**2
        col_total_err = col_total_err**(0.5)
        for iy in range(1, hist.GetNbinsY()+1):
            cur_content = hist.GetBinContent(ix, iy)
            cur_err = hist.GetBinError(ix, iy)
            err = ((cur_err / col_total)**2 + (cur_content * col_total_err / col_total**2)**2)**(0.5)
            hist.SetBinContent(ix, iy, cur_content / col_total)
            hist.SetBinError(ix, iy, err)


if __name__ == '__main__':
    baseline_cuts = "MR > 650 && MR < 4000 && (Rsq > 0.2 || (MR > 1600 && Rsq > 0.1)) && (box==11||box==12||box==14)"
    filters = [  
            "Flag_HBHENoiseFilter",
            "Flag_HBHEIsoNoiseFilter",
            "Flag_goodVertices",
            "Flag_eeBadScFilter",
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_CSCTightHaloFilter",
            "Flag_badChargedCandidateFilter",
            "Flag_badMuonFilter",
            ]
    for filt in filters:
        baseline_cuts += " && " + filt
    baseline_cuts = 'weight*35800*('+baseline_cuts+')'

    f_qcd = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_05Oct2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_QCD_1pb_weighted.root","read")
    tree = f_qcd.Get("RazorInclusive")

    c = rt.TCanvas('c', 'c', 400, 300)

    variables = ['MR', 'Rsq', 'abs(dPhiRazor)', 'nBTaggedJets', 'nSelectedJets']
    titles = {
            'MR':'M_{R}',
            'Rsq':'R^{2}',
            'abs(dPhiRazor)':'#Delta #phi_{razor}',
            'nBTaggedJets':'N_{b-tags}',
            'nSelectedJets':'N_{jets}',
            }
    binning = {
            'MR':[650, 800, 1000, 1200, 1600, 4000],
            'Rsq':[0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30],
            'abs(dPhiRazor)':[0, 2.8, 3.14159],
            'nBTaggedJets':[0, 1, 2, 3, 4],
            'nSelectedJets':[2, 4, 7, 20],
            }
    log_vars = ['MR', 'Rsq']

    for var in variables:
        var_title = titles[var]
        reduced_var = var.replace('(', '').replace(')', '')
        for frac_var in ['nBTaggedJets', 'nSelectedJets']:
            if var == frac_var: continue
            frac_var_title = titles[frac_var]
            hist = rt.TH2F("h{}{}".format(frac_var, reduced_var), 
                    "{} vs {}".format(frac_var_title, var_title),
                    len(binning[var])-1, 
                    array.array('d', binning[var]), 
                    len(binning[frac_var])-1, 
                    array.array('d', binning[frac_var]))
            tree.Draw('{}:{}>>h{}{}'.format(frac_var, var, frac_var, reduced_var), 
                    baseline_cuts)
            normalize_columns(hist)
            hist.SetStats(0)
            hist.GetXaxis().SetTitle(var_title)
            hist.GetYaxis().SetTitle(frac_var_title)
            hist.Draw("colztexte")
            if var in log_vars:
                c.SetLogx()
            else:
                c.SetLogx(False)
            c.Print("QCDMC_{}_{}.pdf".format(frac_var, reduced_var))

