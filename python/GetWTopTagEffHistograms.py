from array import array
import numpy as np
import ROOT as rt

def get_graphs_from_canvas(c):
    graphs = []
    for p in c.GetListOfPrimitives():
        if p.InheritsFrom(rt.TGraphAsymmErrors.Class()):
            graphs.append(p)
    return graphs

def graph_to_hist(g):
    """
    Converts a TGraphAsymmErrors to a TH1F.
    Bin widths are obtained from the x-errors 
    (x error bars are assumed to partition the entire axis range)
    and y-errors are averaged over the upper and lower TGraph errors.
    """
    n = g.GetN()
    # This is a way of accessing the arrays of TGraph points 
    # without segfaults or hanging
    x = np.asarray( [ g.GetX()[i] for i in range(1, n-1) ] )
    y = np.asarray( [ g.GetY()[i] for i in range(1, n-1) ] )
    width = [ g.GetEXhigh()[i] + g.GetEXlow()[i] for i in range(1, n-1) ]
    err = np.asarray( [ (g.GetEYhigh()[i] + g.GetEYlow()[i])/2.0 for i in range(1, n-1) ] )
    start = x[0] - g.GetErrorXlow(1)
    binning = array('d', np.cumsum([start] + width))
    h = rt.TH1F(g.GetName(), g.GetTitle(), n-2, binning)
    for i, (iy, ierr) in enumerate(zip(y, err)):
        h.SetBinContent(i+1, iy)
        h.SetBinError(i+1, ierr)
    return h

if __name__ == '__main__':
    names = ['TagEffFullsim', 'TagEffFastsim', 'TagEffFastsimSF']
    f_w = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/AK8JetTag/2016//WTaggingEfficiency_vs_GenWPtBins_FullFastSim.root")
    c_w = f_w.Get("WTaggingEfficiency_vs_GenWPtBins_FullFastSim")
    g_w = get_graphs_from_canvas(c_w)
    for name, g in zip(names, g_w):
        g.SetName('W'+name)
    h_w = [graph_to_hist(g) for g in g_w]

    f_t = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/AK8JetTag/2016//TopTaggingEfficiency_vs_GenTopPtBins_FullFastSim.root")
    c_t = f_t.Get("TopTaggingEfficiency_vs_GenTopPtBins_FullFastSim")
    g_t = get_graphs_from_canvas(c_t)
    for name, g in zip(names, g_t):
        g.SetName('Top'+name)
    h_t = [graph_to_hist(g) for g in g_t]

    f_out = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/AK8JetTag/2016/AK8WTopTagEff.root", "RECREATE")
    for h in h_w:
        h.Write()
    for h in h_t:
        h.Write()
    f_out.Close()

