import numpy as np
import ROOT as rt

def getXbins(hist_lo, hist_hi):
    """Returns array of bin edges on x-axis for hist"""
    xbins = []
    for i in range(hist_lo.GetNbinsX()):
        xbins.append( hist_lo.GetXaxis().GetXbins()[i] )
    for i in range(hist_hi.GetNbinsX()):
        xbins.append( hist_hi.GetXaxis().GetXbins()[i] )
    xbins.append(4000)
    return np.array(xbins) 

def makeQCDHist(name, xbins):
    """Makes histogram with the given name and bin edges
        on the x-axis"""
    hist = rt.TH1F(name, name, len(xbins)-1, xbins)
    hist.GetXaxis().SetTitle("MR (GeV)")
    return hist


if __name__ == '__main__':
    rt.gROOT.SetBatch()

    inFile = rt.TFile.Open("qcd_data.root")
    outFile = rt.TFile.Open(
            "../../../data/ScaleFactors/RazorMADD2015/RazorQCDScaleFactors_Razor2016_MoriondRereco.root", 
            "recreate")
    c = rt.TCanvas("c","c",800,600)
    for box in ['SevenJet','MultiJet','DiJet']:
        print "Box {}".format(box)
        boxname = box.lower()
        inHists = [ inFile.Get("{}_lo".format(boxname.lower())),
                    inFile.Get("{}_hi".format(boxname).lower())]
        xbins = getXbins(inHists[0], inHists[1])

        slopes = makeQCDHist("QCDSlopes_{}".format(box), xbins)
        inters = makeQCDHist("QCDInters_{}".format(box), xbins)
        covars = makeQCDHist("QCDCovars_{}".format(box), xbins)

        for ibin in range(1, len(xbins)):
            # The high MR region is in a separate histogram
            # with different Rsq binning
            if xbins[ibin-1] >= 1600:
                hist = inHists[1]
                proj_bin = 1
            else:
                hist = inHists[0]
                proj_bin = ibin
            print "{} box, {} < MR < {}".format(
                    box, xbins[ibin-1], xbins[ibin])
            histSlice = hist.ProjectionY("{}_bin{}".format(
                box, ibin), proj_bin, proj_bin)
            if histSlice.Integral() == 0:
                print "Slice has no data"
            else:
                result = histSlice.Fit("pol1", "smf")
                params = result.GetParams()
                errors = result.GetErrors()
                covar   = result.GetCovarianceMatrix()(0,1)
                inters.SetBinContent(ibin, params[0])
                inters.SetBinError(ibin, errors[0])
                slopes.SetBinContent(ibin, params[1])
                slopes.SetBinError(ibin, errors[1])
                covars.SetBinContent(ibin, covar)
                covars.SetBinError(ibin, 0.)
        print "Slopes"
        slopes.Print("all")
        slopes.Draw("colztext")
        c.Print("QCDSlopes_{}.pdf".format(box))
        print "Intercepts"
        inters.Print("all")
        inters.Draw("colztext")
        c.Print("QCDIntercepts_{}.pdf".format(box))
        print "Covariances"
        covars.Print("all")
        covars.Draw("colztext")
        c.Print("QCDCovariances_{}.pdf".format(box))

        outFile.cd()
        slopes.Write()
        inters.Write()
        covars.Write()
