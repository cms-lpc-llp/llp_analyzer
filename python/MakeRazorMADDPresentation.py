### Put all relevant MADD plots into a pdf file for great justice

### Author: 
### Dustin Anderson

import os
import sys
import subprocess
from datetime import date

class Slide(object):
    def __init__(self, slide_title, files):
        """Arguments:
             slide_title: big title for slide
             files: list of file names"""
        self.slide_title = slide_title
        self.files = files
        self.width = 1.0/len(self.files)
        if self.width == 1.0:
            self.width = 0.8 # let's not get too crazy here

    def write(self, f):
        """f: output file object"""
        if not self.check_files: return
        f.write("\\begin{frame} \n")
        f.write("\\begin{center} \n")
        f.write("\\frametitle{%s} \n" % self.slide_title)
        for path in self.files:
            f.write("\\includegraphics[width=%f\\paperwidth]{%s} \n" % (self.width,path))
        f.write("\\end{center} \n")
        f.write("\\end{frame} \n")

    def append_dir_prefix(self, prefix):
        self.files = [ prefix+'/'+path for path in self.files ]

    def check_files(self):
        for f in self.files:
            if not os.path.isfile(f): 
                print "Didn't find file",f
                return False
        return True


def def_slides(plot_dir):
    """Returns: list of slide objects"""
    slides = []
    titles = { "TTJetsSingleLepton":"$t\\bar{t}$+jets control sample",
               "WJetsSingleLepton":"W+jets control sample",
               "WJetsSingleLeptonInv":"W+jets invisible control sample",
               }
    for control_region in ["TTJetsSingleLepton", "WJetsSingleLepton"]:
        slides += [
                Slide(titles[control_region], 
                    ["%s/MR_%s.pdf" % (control_region,control_region),
                     "%s/Rsq_%s.pdf" % (control_region,control_region)]),
                Slide(titles[control_region], 
                    ["%s/MRRsq_%sUnrolledDataMC.pdf" % (control_region,control_region),
                     "%s/%sScaleFactors.pdf" % (control_region,control_region.replace('SingleLepton',''))]),
                Slide(titles[control_region], 
                    ["%s/NBJetsMedium_%s.pdf" % (control_region,control_region),]),
            ]
    control_region = "WJetsSingleLeptonInv"
    slides += [
            Slide(titles[control_region], 
                ["%s/MR_NoW_%s.pdf" % (control_region,control_region),
                 "%s/Rsq_NoW_%s.pdf" % (control_region,control_region)]),
            Slide(titles[control_region], 
                ["%s/MR_NoWRsq_NoW_%sUnrolledDataMC.pdf" % (control_region,control_region),
                 "%s/%sScaleFactors.pdf" % (control_region,control_region.replace('SingleLepton',''))]),
            Slide(titles[control_region], 
                ["%s/NBJetsMedium_%s.pdf" % (control_region,control_region),]),
        ]
    slides += [
            Slide("N$_{jets}$ correction", 
                ["TTJetsForNJets/NJets40_TTJetsSingleLepton.pdf",
                 "WJetsForNJets/NJets40_WJetsSingleLepton.pdf",])
            ]
    for jets in ['DiJet','MultiJet','SevenJet']:
        slides += [
                Slide("1L %s region, after corrections" % jets.lower(), 
                    ["OneLepton%s/MRRsq_SingleLeptonUnrolledDataMC.pdf" % jets,
                     "OneLepton%s/NBJetsMedium_SingleLepton.pdf" % jets,
                        ]),
                ]
        for B in ['0','1','2','3']:
            if jets == 'DiJet' and B == '3': continue
            slides += [
                Slide("1L %s region, after corrections -- %s B" % (jets.lower(),B),
                    ["OneLepton%s%sB/MR_SingleLepton.pdf" % (jets,B),
                     "OneLepton%s%sBMRCorr/Rsq_SingleLepton.pdf" % (jets,B)])
                    ]
    tt2l_dir = "TTJetsDilepton"
    for jets in ['DiJet','MultiJet','SevenJet']:
        suffix = ("MultiJet")*(jets == "MultiJet" or jets == "SevenJet")
        tt2l_title = "t$\\bar{t}$ dilepton %s control region" % jets.lower()
        tt2l_this_dir = tt2l_dir+jets
        slides += [
                Slide(tt2l_title,
                    ["%s/MR_TTJetsDilepton%s.pdf" % (tt2l_this_dir,suffix),
                     "%s/Rsq_TTJetsDilepton%s.pdf" % (tt2l_this_dir,suffix)]),
                Slide(tt2l_title,
                    ["%s/MRRsq_TTJetsDilepton%sUnrolledDataMC.pdf" % (tt2l_this_dir,suffix)]),
                Slide(tt2l_title,
                    ["%s/NJets40_TTJetsDilepton%s.pdf" % (tt2l_this_dir,suffix),
                     "%s/NBJetsMedium_TTJetsDilepton%s.pdf" % (tt2l_this_dir,suffix)]),
                    ]
    photon_dir = "GJetsInv"
    suffixes = ["","ForNJets"]
    title_endings = [""," after razor correction"]
    for suffix,ending in zip(suffixes,title_endings):
        photon_this_dir = photon_dir+suffix
        photon_title = "photon+jets control region"+ending
        if suffix == "":
            slides += [Slide(photon_title,
                ["%s/MR_NoPhoRsq_NoPho_GJetsInvUnrolledDataMC.pdf" % photon_this_dir,
                 "%s/GJetsInvScaleFactors.pdf" % photon_this_dir])]
        else:
            slides += [Slide(photon_title,
                ["%s/MR_NoPhoRsq_NoPho_GJetsInvUnrolledDataMC.pdf" % photon_this_dir])]
        slides += [
                Slide(photon_title,
                    ["%s/MR_NoPho_GJetsInv.pdf" % photon_this_dir,
                     "%s/Rsq_NoPho_GJetsInv.pdf" % photon_this_dir]),
                Slide(photon_title,
                    ["%s/NJets_NoPho_GJetsInv.pdf" % photon_this_dir,
                     "%s/NBJetsMedium_GJetsInv.pdf" % photon_this_dir]),
                    ]
    for jets in ['DiJet','MultiJet','SevenJet']:
        for B in ['0','1','2']:
            if jets == 'DiJet' and B == '3': continue
            slides += [
                Slide("photon+jets %s region, after corrections -- %s B" % (jets.lower(),B),
                    ["GJetsInv%s%sB/MR_NoPho_GJetsInv.pdf" % (jets,B),
                     "GJetsInv%s%sBMRCorr/Rsq_NoPho_GJetsInv.pdf" % (jets,B)])
                    ]
    dy_dir = "DYJetsDileptonInv"
    suffixes = ["Uncorr","","DiJet_Closure","MultiJet_Closure",'SevenJet_Closure']
    title_endings = [" before correction"," after normalization",
            ", dijet closure test",", multijet closure test", ", seven-jet closure test"]
    for suffix,ending in zip(suffixes,title_endings):
        dy_this_dir = dy_dir+suffix
        dy_title = "DY+jets control region"+ending
        dy_suffix = ("MultiJet")*(suffix == "MultiJet_Closure" or suffix == 'SevenJet_Closure')
        if "Jet" in suffix:
            slides += [ Slide(dy_title,
                    ["%s/MR_NoZRsq_NoZ_DYJetsDileptonInv%sUnrolledDataMC.pdf" % (dy_this_dir,dy_suffix)]),
                        Slide(dy_title,
                            ["%s/NJets_NoZ_DYJetsDileptonInv%s.pdf" % (dy_this_dir,dy_suffix),
                             ]),
                            ]
        else:
            slides += [ Slide(dy_title,
                    ["%s/MR_NoZRsq_NoZ_DYJetsDileptonInv%sUnrolledDataMC.pdf" % (dy_this_dir,dy_suffix)
                     ]),
                    Slide(dy_title,
                        ["%s/NJets_NoZ_DYJetsDileptonInv%s.pdf" % (dy_this_dir,dy_suffix),
                         "%s/1_DYJetsDileptonInv%s.pdf" % (dy_this_dir,dy_suffix)]),
                        ]
        slides += [
                Slide(dy_title,
                    ["%s/MR_NoZ_DYJetsDileptonInv%s.pdf" % (dy_this_dir,dy_suffix),
                     "%s/Rsq_NoZ_DYJetsDileptonInv%s.pdf" % (dy_this_dir,dy_suffix)])]
    for leptype in ["Lepton","Tau"]:
        for jets in ["DiJet","MultiJet",'SevenJet']:
            for corr in ["","PtCorr"]:
                veto_dir = "Veto"+leptype+jets+corr
                veto_title = "Veto %s %s control region" % (leptype.lower(),jets.lower())
                veto_region = "Veto"+leptype+"ControlRegion"
                if corr:
                    veto_title += " (after $p_{T}$ correction)"
                slides += [
                        Slide(veto_title,
                            ["%s/lep1_Pt_%s.pdf" % (veto_dir,veto_region),
                             "%s/abslep1_Eta_%s.pdf" % (veto_dir,veto_region)]),
                        Slide(veto_title,
                            ["%s/MRRsq_%sUnrolledDataMC.pdf" % (veto_dir,veto_region)]),
                        Slide(veto_title,
                            ["%s/NJets40_%s.pdf" % (veto_dir,veto_region),
                             "%s/NBJetsMedium_%s.pdf" % (veto_dir,veto_region)]),
                            ]
    for box in ['DiJet','MultiJet','SevenJet','LeptonMultiJet','LeptonSevenJet']:
        for nb in range(4):
            if (box == 'DiJet') and nb == 3: continue
            slides += [
                    Slide("%s (%d b-tags) sideband and signal region"%(box,nb),
                        ["%s/MRRsq_%sUnrolledDataMC.pdf" % (box+str(nb)+"BFineGrained",box)])]

    for slide in slides:
        slide.append_dir_prefix(plot_dir)
    return slides

if __name__ == '__main__':
    if len(sys.argv) < 2:
        in_dir = "Plots/Razor2016_MoriondRereco"
    else:
        in_dir = sys.argv[1]

    fname = "MADD_plots_%s.tex"%date.today()
    with open(fname,'w') as f:
        f.write("\\documentclass{beamer}\n")
        f.write("\\beamertemplatenavigationsymbolsempty\n")
        f.write("\\usepackage{graphicx}\n")
        f.write("\\setbeamersize{text margin left=0pt,text margin right=0pt}\n")
        f.write("\\setbeamertemplate{footline}[frame number]\n")
        f.write("\\begin{document}\n")
        for slide in def_slides(in_dir):
            slide.write(f)
        f.write("\\end{document}\n")

    subprocess.call(['pdflatex', fname])
