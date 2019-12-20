import numpy as np
import ROOT as rt

from QCDFits3D import get_mean_and_cov, get_bins

def combine_hists(hlo, hhi):
    """
    Combines hists hlo and hhi assuming that hlo lives
    at lower MR and higher Rsq than hhi.
    """
    x_binning = np.asarray(get_bins(hlo.GetXaxis())[:-1]
            + get_bins(hhi.GetXaxis())[:])
    y_lo = [y for y in get_bins(hlo.GetYaxis())]
    y_hi = [y for y in get_bins(hhi.GetYaxis()) if y <= min(y_lo)]
    y_binning = np.asarray(y_hi[:-1] + y_lo[:])
    comb_hist = rt.TH2F(hlo.GetName()+hhi.GetName(), 'hist',
            len(x_binning)-1, x_binning,
            len(y_binning)-1, y_binning)
    for ix in range(1, comb_hist.GetNbinsX()+1):
        for iy in range(1, comb_hist.GetNbinsY()+1):
            mr = comb_hist.GetXaxis().GetBinCenter(ix)
            rsq = comb_hist.GetYaxis().GetBinCenter(iy)
            if mr > hlo.GetXaxis().GetXmax():
                this_h = hhi
            else:
                this_h = hlo
            ibin = this_h.FindFixBin(mr, rsq)
            if this_h.GetBinError(ibin) > 0:
                comb_hist.SetBinContent(ix, iy,
                        this_h.GetBinContent(ibin))
                comb_hist.SetBinError(ix, iy,
                        this_h.GetBinError(ibin))
    c = rt.TCanvas("c", "c", 800, 600)
    comb_hist.SetStats(0)
    comb_hist.Draw("colztexte")
    c.Print("test.pdf")
    return comb_hist

def integral_error(tf2, result,
        x1, x2, y1, y2):
    n_toys = 1000
    mean, cov = get_mean_and_cov(tf2, result)

    toys = np.zeros(n_toys)
    for itoy in range(n_toys):
        rand = np.random.multivariate_normal(mean, cov)
        tf2.SetParameters(rand)
        integral = tf2.Integral(x1, x2, y1, y2)
        toys[itoy] = integral
    tf2.SetParameters(mean)
    return toys.std()

def fill_fit_slice(hist, tf2, result, xlow, xhigh):
    for ibin in range(1, hist.GetNbinsX()+1):
        ylow = hist.GetXaxis().GetBinLowEdge(ibin)
        yhigh = hist.GetXaxis().GetBinUpEdge(ibin)
        integral = tf2.Integral(xlow, xhigh, ylow, yhigh)
        integral_err = integral_error(tf2, result, xlow, xhigh, ylow, yhigh)
        area = (xhigh-xlow)*(yhigh-ylow)
        hist.SetBinContent(ibin, integral/area)
        hist.SetBinError(ibin, integral_err/area)

def plot_slice(obs_hist, fit_hist):
    c = rt.TCanvas(obs_hist.GetName(), "c", 800, 600)
    obs_hist.SetStats(0)
    obs_hist.Draw()
    fit_hist.SetLineColor(rt.kRed)
    fit_hist.Draw("same")
    leg = rt.TLegend(0.1, 0.7, 0.5, 0.9)
    leg.AddEntry(obs_hist, "QCD Transfer Factors")
    leg.AddEntry(fit_hist, "Fitted prediction")
    leg.Draw()
    c.Print('njetsqcdfit_'+obs_hist.GetName()+'.pdf')

def plot_slices(hist, fit, result):
    for ibin in range(1, hist.GetNbinsX()+1):
        this_slice = hist.ProjectionY("{}_bin{}".format(
            hist.GetName(), ibin), ibin, ibin, "e")
        this_slice.GetXaxis().SetTitle("R^{2}")
        this_slice.GetYaxis().SetTitle("QCD Transfer factor")
        this_slice.SetTitle("")
        this_fit = this_slice.Clone("{}_fit{}".format(
            hist.GetName(), ibin))
        low = hist.GetXaxis().GetBinLowEdge(ibin)
        high = hist.GetXaxis().GetBinUpEdge(ibin)
        this_fit.Reset()
        fill_fit_slice(this_fit, fit, result, low, high)
        plot_slice(this_slice, this_fit)

def plot_pull_distribution(dist):
    c = rt.TCanvas("c"+dist.GetName(), "c", 800, 600)
    dist.SetLineWidth(2)
    dist.Fit('gaus')
    dist.Draw()
    c.Print(dist.GetName()+".pdf")

def make_pull_distribution(hist, fit, result):
    pull_hist = rt.TH1F("pulls"+hist.GetName(), "Pull distribution", 25, -5, 5)
    for ix in range(1, hist.GetNbinsX()+1):
        xlow = hist.GetXaxis().GetBinLowEdge(ix)
        xhigh = hist.GetXaxis().GetBinUpEdge(ix)
        for iy in range(1, hist.GetNbinsY()+1):
            ylow = hist.GetYaxis().GetBinLowEdge(iy)
            yhigh = hist.GetYaxis().GetBinUpEdge(iy)
            area = (xhigh-xlow)*(yhigh-ylow)
            fitted = fit.Integral(xlow, xhigh, ylow, yhigh)/area
            obs = hist.GetBinContent(ix, iy)
            obs_err = hist.GetBinError(ix, iy)
            if obs_err == 0:
                continue
            pull = (obs-fitted)/obs_err
            pull_hist.Fill(pull)
    plot_pull_distribution(pull_hist)

def make_tf2_sevenjet():
    box = 'SevenJet'
    names = ['Intercept', 'Y slope']
    multipliers = [1, 1]
    num_pars = len(names)
    pars = range(num_pars)
    function_expr = "[{}] + [{}]*y".format(*pars)
    print "Fit function: "+function_expr
    return (rt.TF2("qcd{}".format(box.lower()), function_expr), 
            num_pars, names, multipliers)

def make_tf2(box):
    """
    Returns 2D fit function and the number of free fit parameters,
    the list of parameter names, and a list of multipliers 
    for display convenience
    """
    if box.lower() == 'sevenjet':
        return make_tf2_sevenjet()
    names = ['Intercept', 'X slope', 'Y slope', 'XY slope']
    multipliers = [1, 1000, 1, 1000]
    num_pars = len(names)
    pars = range(num_pars)
    function_expr = "[{}] + [{}]*x + ([{}] + [{}]*x)*y".format(*pars)
    print "Fit function: "+function_expr
    return (rt.TF2("qcd{}".format(box.lower()), function_expr), 
            num_pars, names, multipliers)

def export_fit(fun, result, box):
    out_name = 'qcd_best_fit_2d_{}.root'.format(box.lower())
    out_f = rt.TFile.Open(out_name, 'recreate')
    fun.Write()
    out_f.Close()
    mean, cov = get_mean_and_cov(fun, result)
    out_name = out_name.replace('.root', '.npz')
    with open(out_name, 'w') as out_f:
        np.savez(out_f, mean, cov)


if __name__ == '__main__':
    rt.gROOT.SetBatch()
    rt.gStyle.SetOptFit(1)

    inFile = rt.TFile.Open("qcd_data.root")
    for box in ['SevenJet','MultiJet','DiJet']:
        print "Box {}".format(box)
        boxname = box.lower()
        hist_lo = inFile.Get("{}_lo".format(boxname))
        hist_hi = inFile.Get("{}_hi".format(boxname))
        hist = combine_hists(hist_lo, hist_hi)
        fun, n_pars, names, multipliers = make_tf2(box)

        result = hist.Fit("qcd{}".format(box.lower()), "SMF")
        params = result.GetParams()
        errors = result.GetErrors()
        pars = [params[i] for i in range(n_pars)]
        errs = [errors[i] for i in range(n_pars)]

        names += ['Chi Square / NDOF', 'Probability']
        multipliers += [1, 1]
        pars.append(fun.GetChisquare() / fun.GetNDF())
        pars.append(fun.GetProb())
        errs += [0, 0]

        plot_slices(hist, fun, result)
        make_pull_distribution(hist, fun, result)
        for par, err, parname, multiplier in zip(pars, errs,
                names, multipliers):
            print "{}: {:.3f} +/- {:.2f}".format(
                    parname, par*multiplier, err*multiplier)
        export_fit(fun, result, box)
