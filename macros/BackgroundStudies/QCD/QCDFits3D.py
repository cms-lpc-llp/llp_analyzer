import copy
import argparse
import array
import numpy as np
import ROOT as rt

def get_mean_and_cov(fun, result):
    npars = fun.GetNumberFreeParameters()
    best_pars = result.GetParams()
    cov_from_fit = result.GetCovarianceMatrix()

    mean = np.zeros(npars)
    cov = np.zeros((npars, npars))
    for ix in range(npars):
        mean[ix] = best_pars[ix]
        for iy in range(ix+1):
            cov[ix, iy] = cov_from_fit(ix, iy)
            cov[iy, ix] = cov_from_fit(ix, iy)
    return mean, cov

def integral_error(tf3, result,
        x1, x2, y1, y2, z1, z2):
    n_toys = 1000
    mean, cov = get_mean_and_cov(tf3, result)

    toys = np.zeros(n_toys)
    for itoy in range(n_toys):
        rand = np.random.multivariate_normal(mean, cov)
        tf3.SetParameters(rand)
        integral = tf3.Integral(x1, x2, y1, y2, z1, z2)
        toys[itoy] = integral
    tf3.SetParameters(mean)
    return toys.std()

def fill_fit_slice(hist, tf3, result, xlow, xhigh, njets):
    for ibin in range(1, hist.GetNbinsX()+1):
        ylow = hist.GetXaxis().GetBinLowEdge(ibin)
        yhigh = hist.GetXaxis().GetBinUpEdge(ibin)
        zlow = njets - 0.5
        zhigh = njets + 0.5
        integral = tf3.Integral(xlow, xhigh, ylow, yhigh, zlow, zhigh)
        integral_err = integral_error(tf3, result, xlow, xhigh, ylow, yhigh, zlow, zhigh)
        area = (xhigh-xlow)*(yhigh-ylow)*(zhigh-zlow)
        hist.SetBinContent(ibin, integral/area)
        hist.SetBinError(ibin, integral_err/area)

def plot_slice(obs_hist, fit_hist):
    c = rt.TCanvas(obs_hist.GetName(), "c", 800, 600)

    obs_hist.SetStats(0)
    obs_hist.Draw()

    fit_hist.SetLineColor(rt.kRed)
    fit_hist.Draw("same")

    c.Print('njets3dqcdfit_'+obs_hist.GetName()+'.pdf')

def plot_slices(hist, fit, result, njets):
    for ibin in range(1, hist.GetNbinsX()+1):
        njets_bin = hist.GetZaxis().FindFixBin(njets)
        this_slice = hist.ProjectionY("{}_{}_bin{}".format(
            hist.GetName(), njets, ibin), ibin, ibin, 
            njets_bin, njets_bin, "e")
        this_fit = this_slice.Clone("{}_{}_fit{}".format(
            hist.GetName(), njets, ibin))
        low = hist.GetXaxis().GetBinLowEdge(ibin)
        high = hist.GetXaxis().GetBinUpEdge(ibin)
        this_fit.Reset()
        fill_fit_slice(this_fit, fit, result, low, high, njets)
        plot_slice(this_slice, this_fit)

def combine_th2s(hists, z_bins):
    assert len(z_bins) == len(hists) + 1
    ref_hist = hists[0]
    x_bins = get_bins(ref_hist.GetXaxis())
    y_bins = get_bins(ref_hist.GetYaxis())
    z_bins = array.array('d', z_bins)
    hist3d = rt.TH3D(ref_hist.GetName().replace('nj2','')+"3d", "QCD",
            len(x_bins)-1, x_bins,
            len(y_bins)-1, y_bins,
            len(z_bins)-1, z_bins)
    for z, hist in enumerate(hists):
        iz = z + 1
        for ix in range(1, len(x_bins)):
            for iy in range(1, len(y_bins)):
                content = hist.GetBinContent(ix, iy)
                if content > 10:
                    hist3d.SetBinContent(ix, iy, iz, 0)
                else:
                    hist3d.SetBinContent(ix, iy, iz, 
                            content) 
                hist3d.SetBinError(ix, iy, iz,
                        hist.GetBinError(ix, iy))
    return hist3d

def get_bins(axis):
    return array.array('d',
            [ axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins()+2) ])

def fill_par_numbers(expr, pars):
    n_pars_needed = expr.count('{')
    pars_to_use = [pars.pop() for _ in range(n_pars_needed)]
    return expr.format(*pars_to_use)

def plot_pull_distribution(dist, fit=True):
    c = rt.TCanvas("c"+dist.GetName(), "c", 800, 600)
    dist.SetLineWidth(2)
    if fit:
        dist.Fit('gaus')
    dist.Draw()
    c.Print(dist.GetName()+".pdf")

def make_pull_distributions(hist, fit, result):
    pull_hist_all = rt.TH1F("pulls_all", "pulls_all", 25, -5, 5)
    nj = np.asarray(get_bins(hist.GetZaxis()))
    # center
    nj = (nj[1:] + nj[:-1])/2.0
    pull_hists = []
    for n in nj:
        pull_hists.append(rt.TH1F("pulls_{}j".format(n), "pulls", 25, -5, 5))
        zlow = n - 0.5
        zhigh = n + 0.5
        iz = hist.GetZaxis().FindFixBin(n)
        for ix in range(1, hist.GetNbinsX()+1):
            xlow = hist.GetXaxis().GetBinLowEdge(ix)
            xhigh = hist.GetXaxis().GetBinUpEdge(ix)
            for iy in range(1, hist.GetNbinsY()+1):
                ylow = hist.GetYaxis().GetBinLowEdge(iy)
                yhigh = hist.GetYaxis().GetBinUpEdge(iy)
                area = (xhigh-xlow)*(yhigh-ylow)*(zhigh-zlow)
                fitted = fit.Integral(xlow, xhigh, ylow, yhigh, zlow, zhigh)/area
                obs = hist.GetBinContent(ix, iy, iz)
                obs_err = hist.GetBinError(ix, iy, iz)
                pull = (obs-fitted)/obs_err
                pull_hists[-1].Fill(pull)
                pull_hist_all.Fill(pull)
        plot_pull_distribution(pull_hists[-1])
    plot_pull_distribution(pull_hist_all)

def make_tf3():
    """
    Returns 3D fit function and the number of free fit parameters,
    the list of parameter names, and a list of multipliers 
    for display convenience
    """
    names = ['Intercept', 'Z slope', 'X slope', 
             'XZ slope', 'Y slope', 'YZ slope', 'YZZ slope']
    multipliers = [1, 1, 1000, 1000, 1, 1, 1]
    num_pars = len(names)
    pars = range(num_pars-1, -1, -1)

    par_templates = [
            "[{}] + [{}]*(z-2)",
            "[{}] + [{}]*(z-2)",
            "[{}] + [{}]*(z-2) + [{}]*(z-2)*(z-2)",
            ]
    intercept, xslope, yslope = [
            fill_par_numbers(temp, pars) for temp in par_templates]
    function_expr = "{} + ({})*(x-400) + ({})*(y-0.2)".format(
        intercept, xslope, yslope)
    print "Fit function: "+function_expr
    return rt.TF3("qcd", function_expr), num_pars, names, multipliers

def export_fit(fun, result, exclusive=False):
    out_name = 'qcd_best_fit_3d.root'
    if args.exclusive:
        out_name = out_name.replace('.root', '_exclusive.root')
    out_f = rt.TFile.Open(out_name, "recreate")
    fun.Write()
    out_f.Close()
    mean, cov = get_mean_and_cov(fun, result)
    out_name = out_name.replace('.root', '.npz')
    with open(out_name, 'w') as out_f:
        np.savez(out_f, mean, cov)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exclusive', action='store_true',
            help='use exclusive njets transfer factors')
    args = parser.parse_args()
    rt.gROOT.SetBatch()
    rt.gStyle.SetOptFit(1)

    tag = 'inclusive'
    if args.exclusive:
        tag = 'exclusive'
    inFile = rt.TFile.Open("qcd_njets_data_{}.root".format(tag))
    nj = [2,3,4,5,6,7]
    nj_binning = [n-0.5 for n in nj]+[nj[-1]+0.5]

    inHists = [ inFile.Get("nj{}".format(n)) for n in nj ]
    hist3d = combine_th2s(inHists, nj_binning)
    if args.exclusive:
        hist3d.SetName(hist3d.GetName()+'_exclusive')

    fun, n_pars, names, multipliers = make_tf3()

    result = hist3d.Fit("qcd", "SMF")
    params = result.GetParams()
    errors = result.GetErrors()
    pars = [params[i] for i in range(n_pars)]
    errs = [errors[i] for i in range(n_pars)]

    names += ['Chi Square / NDOF', 'Probability']
    multipliers += [1, 1]
    pars.append(fun.GetChisquare() / fun.GetNDF())
    pars.append(fun.GetProb())
    errs += [0, 0]

    for n in nj:
        plot_slices(hist3d, fun, result, n)
    make_pull_distributions(hist3d, fun, result)
    for par, err, parname, multiplier in zip(pars, errs, 
            names, multipliers):
        print "{}: {:.3f} +/- {:.2f}".format(
                parname, par*multiplier, err*multiplier)
    export_fit(fun, result, args.exclusive)
