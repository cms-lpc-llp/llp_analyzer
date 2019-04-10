import array
import ROOT as rt
rt.gROOT.SetBatch()
rt.gStyle.SetPaintTextFormat(".3f")

import QCDFits2D as fits

if __name__ == '__main__':
    rt.gROOT.SetBatch()
    rt.gStyle.SetOptFit(1)

    baseline_cuts = "MR > 550 && MR < 4000 && (Rsq > 0.2 || (MR > 1600 && Rsq > 0.1)) && (box==11||box==12||box==14)"
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

    f_qcd = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_27Nov2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_QCD_1pb_weighted.root","read")
    tree = f_qcd.Get("RazorInclusive")

    c = rt.TCanvas('c', 'c', 400, 300)

    binning = {
            'MR':[550, 650, 800, 1000, 1200, 1400, 1600, 4000],
            'Rsq':[0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30],
            'Rsq_high':[0.3, 0.41, 0.52, 0.64, 1.50],
            }
    dphi_cuts = {
            'lo':'abs(dPhiRazor)<2.8',
            'hi':'abs(dPhiRazor)>2.8',
            }
    hists = {}
    high_hists = {} # these are for plotting the Rsq extrapolation
    for region in ['lo', 'hi']:
        hists[region] = rt.TH2F("h{}".format(region), 'tfs',
                len(binning['MR'])-1, 
                array.array('d', binning['MR']),
                len(binning['Rsq'])-1,
                array.array('d', binning['Rsq']))
        high_hists[region] = rt.TH2F("h{}_high".format(region), 'tfs',
                len(binning['MR'])-1, 
                array.array('d', binning['MR']),
                len(binning['Rsq_high'])-1,
                array.array('d', binning['Rsq_high']))
        cuts = "{} && {}".format(baseline_cuts, dphi_cuts[region])
        tree.Draw('Rsq:MR>>h{}'.format(region), cuts)
        tree.Draw('Rsq:MR>>h{}_high'.format(region), cuts)
    tfs = hists['lo'].Clone("qcd_mc")
    tfs.Divide(hists['hi'])
    tfs_high = high_hists['lo'].Clone("qcd_mc_high")
    tfs_high.Divide(high_hists['hi'])
    fun, n_pars, names, multipliers = fits.make_tf2('mc')
    result = tfs.Fit("qcdmc", "SMF")
    params = result.GetParams()
    errors = result.GetErrors()
    pars = [params[i] for i in range(n_pars)]
    errs = [errors[i] for i in range(n_pars)]

    names += ['Chi Square / NDOF', 'Probability']
    multipliers += [1, 1]
    pars.append(fun.GetChisquare() / fun.GetNDF())
    pars.append(fun.GetProb())
    errs += [0, 0]

    fits.plot_slices(tfs, fun, result)
    fits.plot_slices(tfs_high, fun, result)
    fits.make_pull_distribution(tfs, fun, result)
    for par, err, parname, multiplier in zip(pars, errs,
            names, multipliers):
        print "{}: {:.3f} +/- {:.2f}".format(
                parname, par*multiplier, err*multiplier)
    fits.export_fit(fun, result, 'mc')
