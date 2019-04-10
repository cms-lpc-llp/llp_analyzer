### This script performs two checks to derive uncertainties on the QCD prediction:
# 1) Compares fit predictions to actual transfer factors in the low-MR, high-Rsq
#    bins to assess the goodness of the linearity assumption
# 2) Does the same check in the high-MR, low-Rsq region to assess whether
#    it's appropriate to extrapolate the fit to that region

import numpy as np
import ROOT as rt

from QCDFits3D import plot_pull_distribution

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

def load_hists(mode, region):
    if mode == 'njets':
        in_name = 'qcd_njets_data_inclusive.root'
        hist_names = ['nj{}{}'.format(nj, region) for nj in range(2, 8)]
    else:
        in_name = 'qcd_data.root'
        hist_names = ['{}_{}'.format(mode, region.replace('mrhi', 'hi'))]
    in_file = rt.TFile.Open(in_name)
    hists = []
    for h in hist_names:
        hist = in_file.Get(h)
        hist.SetDirectory(0)
        hists.append(hist)
    return hists

def read_mean_and_cov(in_name):
    """
    Retrieves mean vector and covariance matrix
    from first two lines of text file and returns them
    as numpy arrays.
    """
    with open(in_name, 'r') as in_f:
        arrs = np.load(in_f)
        mean = arrs['arr_0']
        cov = arrs['arr_1']
    return mean, cov
    
def load_function(mode):
    if mode == 'njets':
        in_name = 'qcd_best_fit_3d'
        fun_name = 'qcd'
    else:
        in_name = 'qcd_best_fit_2d_'+mode
        fun_name = 'qcd'+mode
    in_root = rt.TFile.Open(in_name+'.root')
    fun = in_root.Get(fun_name)
    mean, cov = read_mean_and_cov(in_name+'.npz')
    return fun, mean, cov

def get_pred(hist, fun, mean, cov, mode, nj=None):
    pred = hist.ProjectionY(hist.GetName()+'_pred')
    pred.Reset()
    for ibin in range(1, hist.GetNbinsY()+1):
        xmin = hist.GetXaxis().GetBinLowEdge(1)
        xmax = hist.GetXaxis().GetBinUpEdge(1)
        ymin = hist.GetYaxis().GetBinLowEdge(ibin)
        ymax = hist.GetYaxis().GetBinUpEdge(ibin)
        coords = [xmin, xmax, ymin, ymax]
        area = (xmax-xmin)*(ymax-ymin)
        if mode == 'njets':
            zmin = nj-0.5
            zmax = nj+0.5
            coords += [zmin, zmax]
            area *= (zmax-zmin)
        function_val = fun.Integral(*coords)/area
        function_err = integral_error(fun, mean, cov, *coords)/area
        pred.SetBinContent(ibin, function_val)
        pred.SetBinError(ibin, function_err)
    return pred

def make_pull_distribution(hist, fit):
    pull_hist = rt.TH1F("pulls"+hist.GetName(), 
            "Pull distribution", 25, -5, 5)
    for iy in range(1, hist.GetNbinsY()+1):
        fitted = pred.GetBinContent(iy)
        obs = hist.GetBinContent(1, iy)
        obs_err = hist.GetBinError(1, iy)
        if obs_err == 0:
            continue
        pull = (obs-fitted)/obs_err
        pull_hist.Fill(pull)
        if abs(pull) > 1:
            print "Diff {}, Obs err {}, Pred err {}, Pull {}".format(
                    obs-fitted, obs_err, pred.GetBinError(iy), pull)
            print "For agreement within one sigma, would need to multiply pred by a factor of {}".format(
                    (obs-obs_err)/pred.GetBinContent(iy))
    plot_pull_distribution(pull_hist, fit=False)

def make_plot(hist, pred):
    name = hist.GetName()
    c = rt.TCanvas("c"+name, "c", 800, 600)
    hist1d = hist.ProjectionY(name+'proj', 1, 1, "e")
    hist1d.GetXaxis().SetTitle("R^{2}")
    hist1d.GetYaxis().SetTitle("QCD Transfer Factor")
    hist1d.SetTitle("")
    hist1d.SetStats(0)
    hist1d.SetMaximum(
            min(15, 
            max(hist1d.GetBinContent(hist1d.GetMaximumBin())*1.5,
                pred.GetBinContent(pred.GetMaximumBin())*1.5))
            )
    hist1d.SetMinimum(min(hist1d.GetBinContent(hist1d.GetMinimumBin())-0.5,
        pred.GetBinContent(pred.GetMinimumBin())-0.5))
    hist1d.Draw()
    pred.SetStats(0)
    pred.SetLineColor(rt.kRed)
    pred.Draw('same')
    leg = rt.TLegend(0.1, 0.7, 0.5, 0.9)
    leg.AddEntry(hist1d, "QCD Transfer Factors")
    leg.AddEntry(pred, "Fitted prediction")
    leg.Draw()
    c.Print("closure_{}.pdf".format(name))


if __name__ == '__main__':
    rt.gROOT.SetBatch()
    regions = ['mrsideband']
    modes = ['dijet', 'multijet', 'sevenjet']

    # We have histograms binned in njets (for extrapolation to 7jet box)
    # and histograms for the dijet (2-3) and multijet (4-6) jet boxes.
    # For the njets binned case we make predictions using a TF3 in MR, Rsq, Njets
    # (and have one histogram for each njets slice).
    # For the other cases we use a TF2 in MR and Rsq.
    for mode in modes:
        if mode == 'njets':
            njets = range(2, 8)
        else:
            njets = [None for _ in range(2, 8)]
        for region in regions:
            hists = load_hists(mode, region)
            fun, mean, cov = load_function(mode)

            preds = [get_pred(hist, fun, mean, cov, mode, nj) 
                    for hist, nj in zip(hists, njets)]
            for hist, pred in zip(hists, preds):
                make_plot(hist, pred)
                make_pull_distribution(hist, pred)
