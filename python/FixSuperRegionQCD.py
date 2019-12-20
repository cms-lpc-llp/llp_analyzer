import subprocess as sp
import ROOT as rt

import macro.macro as macro

def fix_super_hist(super_hist, standard_hist):
    """
    Creates a new version of super_hist where
    the bin contents are obtained by aggregating the appropriate
    bins in standard_hist, with uncertainties summed 
    in quadrature over the aggregated bins.
    This fixes an issue where a very large QCD transfer factor
    is used in the super regions due to the contribution from
    the highest Rsq bin, which is oversized.
    Assumes 2D histograms.
    """
    contents = [0.] * super_hist.GetSize()
    errors = [0.] * super_hist.GetSize()
    for bx in range(1, standard_hist.GetNbinsX()+1):
        for by in range(1, standard_hist.GetNbinsY()+1):
            this_content = standard_hist.GetBinContent(bx, by)
            this_error = standard_hist.GetBinError(bx, by)

            x = standard_hist.GetXaxis().GetBinCenter(bx)
            y = standard_hist.GetYaxis().GetBinCenter(by)
            super_bin = super_hist.FindBin(x, y)

            err = errors[super_bin-1]
            contents[super_bin-1] += this_content
            errors[super_bin-1] = (err*err + this_error*this_error)**(0.5)
    for b in range(1, super_hist.GetSize()+1):
        super_hist.SetBinContent(b, contents[b-1])
        super_hist.SetBinError(b, errors[b-1])


if __name__ == '__main__':
    boxes = [
        #'DiJet', 
        'MultiJet', 
        'SevenJet'
        ]
    btags = [0, 1, 2, 3]
    tag = 'Razor2016_MoriondRereco'

    for box in boxes:
        for nb in btags:
            if box == 'DiJet' and nb == 3: continue

            hist_dir = 'Plots/{}/{}{}BFineGrainedSuper/'.format(
                    tag, box, nb)
            hist_file = 'razorHistograms{}{}BFineGrainedSuper.root'.format(
                    box, nb)
            super_fname = hist_dir + hist_file
            standard_fname = super_fname.replace('Super', '')

            # Create a backup of the histogram file
            old_super_fname = super_fname.replace('.root', '_BKP.root')
            sp.call(['cp', super_fname, old_super_fname])

            # Get both super and non-super region histograms
            super_hists = macro.importHists(super_fname)
            standard_hists = macro.importHists(standard_fname)

            # First do the nominal QCD histogram
            var = ('MR', 'Rsq')
            proc = 'QCD'
            fix_super_hist(
                    super_hists[proc][var],
                    standard_hists[proc][var])
            # Then do all systematic histograms
            for syst in super_hists['Sys'][proc]:
                fix_super_hist(
                        super_hists['Sys'][proc][syst][var],
                        standard_hists['Sys'][proc][syst][var])

            # Output file
            macro.exportHists(super_hists, super_fname)
