import ROOT as rt
import csv
import re
import sys
import collections
from collections import OrderedDict
import uproot
import time
from matplotlib import pyplot as plt
sys.path.append('/nfshome/christiw/llp/delayed_jet_analyzer/lib/')
from histo_utilities import create_TH1D, create_TH2D, create_TGraph,std_color_list

import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 0

wH = 1
Z_MASS = 91.2
lumi = 137000 #in pb-1
lepid = 13 #muon

mindeltaR_cuts =[0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15]
time_cuts = [0.5,0.6,0.7,0.8,0.9]
gamma_cuts = [0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15]






fpath_bkg =OrderedDict()
tree_bkg = OrderedDict()
path = '/mnt/hadoop/store/group/phys_exotica/delayedjets/llp_analyzer/V1p4/MC_Summer16/v1/'
bkg_path = path+"/bkg/wH/normalized/"

fpath_bkg['bkg_QCD'] = bkg_path+"QCD_HT50toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root"
fpath_bkg['bkg_DYJetsToLL'] = bkg_path+"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root"
fpath_bkg['bkg_TTJets_DiLept'] = bkg_path+"TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root"
fpath_bkg['bkg_TTJets_SingleLeptFromTbar'] = bkg_path+"TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root"
fpath_bkg['bkg_TTJets_SingleLeptFromT'] = bkg_path+"TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root"
fpath_bkg['bkg_WJetsToLNu'] = bkg_path+"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root"
fpath_bkg['bkg_ZJetsToNuNu'] = bkg_path+"ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root"

if wH:
    sig_path = path+'/signals/wH/normalized/'
else:    
    sig_path = path+'/signals/zH/normalized/'

fpath_bkg['bbbb'] = sig_path + 'WplusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_1pb_weighted.root'
# fpath_sig['met+bb'] = sig_path + 'ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh125_mx50_pl1000_ev100000_1pb_weighted.root'

for k,v in fpath_bkg.items():
    root_dir = uproot.open(v) 
    tree_bkg[k] = root_dir['vH']
    a = tree_bkg[k]["weight"].array()

legend = {}
legend['bbbb'] = 'Signal (m_{h} ,m_{x})=(125, 50) GeV, c#tau = 1 m'
legend['bbbb_mh2000_ctau1'] = 'Signal (m_{h} ,m_{x})=(2000, 975) GeV, c#tau = 1 m'
legend['bbbb_mh125_ctau10'] = 'Signal (m_{h} ,m_{x})=(125, 50) GeV, c#tau = 10 m'
legend['bbbb_mh2000_ctau10'] = 'Signal (m_{h} ,m_{x})=(2000, 975) GeV, c#tau = 10 m'
legend['bkg_QCD'] = 'QCD'
legend['bkg_DYJetsToLL'] = 'DYJetsToLL'
legend['bkg_TTJets_DiLept'] = 'TTJets_DiLept'
legend['bkg_TTJets_SingleLeptFromTbar'] = 'TTJets_SingleLeptFromTbar'
legend['bkg_TTJets_SingleLeptFromT'] = 'TTJets_SingleLeptFromT'
legend['bkg_WJetsToLNu'] = 'WJetsToLNu'
legend['bkg_ZJetsToNuNu'] = 'ZJetsToNuNu'



trigger_names_file = '/nfshome/christiw/llp/delayed_jet_analyzer/data/trigger_names_llp_v1.dat'
trigger_names = []
with open(trigger_names_file) as f:
    reader = csv.reader(f, delimiter=" ")
    for line in reader:
        trigger_names.append(line[2])
if wH:
    trigger_paths = [87,135] #PFMET120
else:
    trigger_paths = [177,362,87,135] #PFMET120





jetPt = {}
jetEta = {}
jetPhi = {}
jetTime = {}

lepPt = {}
lepEta = {}
lepPhi = {}
start_t = time.time()
jetGammaMax = {}
jetMinDeltaRPVTracks = {}
jetchef = {}
met = {}

mass = {} #mt if wH, Mz if zH
weight = {}
start_t = time.time()

for k in tree_bkg.keys():
    if True:
        T = tree_bkg[k]
        
        #choose only muon or ele
        lepPt[k] = T.array('lepPt')[np.abs(T.array('lepPdgId')) == lepid]
        
        #jet selection
        jet_sel = np.abs(T.array('jetEta'))<1.48
        jetEta[k] =T.array('jetEta')[jet_sel]
        jetPt[k] = T.array('jetPt')[jet_sel]
        jetTime[k] =T.array('jetTime')[jet_sel]
        jetGammaMax[k] = T.array('jetGammaMax_ET')[jet_sel]
        jetMinDeltaRPVTracks[k] = T.array('jetMinDeltaRPVTracks')[jet_sel]
        jetchef[k] = T.array('jetChargedHadronEnergyFraction')[jet_sel]
    
        # event level selection
        print(time.time()-start_t)

        hlt = T['HLTDecision'].array()
        # select only triggered events
        sel_ev= np.zeros(hlt[:,0].shape)
        for tr in trigger_paths:
            sel_ev  = np.logical_or(sel_ev,hlt[:,tr])
        sel_ev = np.logical_and(sel_ev,lepPt[k].counts == 1)
        sel_ev = np.logical_and(sel_ev,jetPt[k].counts >= 1)

        if wH:
            mass_temp = T.array('MT')
        else:
            mass_temp = T.array('ZMass')

        met[k] = T.array('met')[sel_ev]
        weight[k] =  T.array('weight')[sel_ev] * lumi
        mass[k] = mass_temp[sel_ev]
        
        lepPt[k] = lepPt[k][sel_ev]
        jetEta[k] =jetEta[k][sel_ev]
        jetPt[k] = jetPt[k][sel_ev]
        jetTime[k] =jetTime[k][sel_ev]
        jetGammaMax[k] = jetGammaMax[k][sel_ev]
        jetMinDeltaRPVTracks[k] = jetMinDeltaRPVTracks[k][sel_ev]
        jetchef[k] = jetchef[k][sel_ev]
    
# loop over all the cuts
# loop over all the cuts and save the efficiency for each cut and a list for cuts

nTags = {}
nJets = {}
eff_sig = []
eff_bkg = []
purity = []
cuts = []
with open('output.csv','wb') as csvfile:
    filewriter = csv.writer(csvfile,delimiter = ',')

    for gamma_cut in gamma_cuts: #loop over the three variables
        for time_cut in time_cuts:
             for mindeltaR_cut in mindeltaR_cuts:
            #go through each sample
                num_temp = 0.0
                denom_temp = 0.0
                for k,T in tree_bkg.items():
                    tagged_jet = np.logical_and(jetGammaMax[k]<gamma_cut,jetTime[k]>time_cut)
                    tagged_jet = np.logical_and(tagged_jet,np.abs(jetEta[k])<1.48)
                    nTags[k] = jetPt[k][tagged_jet].count()
                    nJets[k] = jetPt[k].count() #njets with preselection
                    start_t = time.time()
                    if k[:4] == 'bbbb':
                        eff_sig.append(1.0*np.sum(nTags[k])/np.sum(nJets[k]))
                    else:
                        num_temp += np.sum(nTags[k])
                        denom_temp += np.sum(nJets[k])
                eff_bkg.append(1.0*num_temp/denom_temp)
                cuts.append([gamma_cut,time_cut])
                filewriter.writerow([gamma_cut, time_cut,eff_bkg[-1],eff_sig[-1]])
                        
#c = rt.TCanvas('c','c', 800, 800)
#leg = rt.TLegend(0.58,0.16,0.92,0.33)

#leg.SetTextSize(0.022)
#leg.SetEntrySeparation(0.01)
#gr = {}
#gr['bbbb'] = create_TGraph(eff_bkg,eff_sig,axis_title = ['#epsilon_{bkg}','#epsilon_{signal}'])
#gr['bbbb'].SetMarkerColor(std_color_list[0])
#gr['bbbb'].SetMarkerSize(1)
#gr['bbbb'].SetMarkerStyle(20)
#gr['bbbb'].GetXaxis().SetRangeUser(0.0, 1.0)
#gr['bbbb'].Draw('ap')
#c.Draw()
#write to csv file
c.SaveAs("../plots/ntags_optimization/mindeltaRcut"+str(mindeltaR_cut)+"_timecut"+str(time_cut)+"gammacut"+str(gamma_cut)+".png")



















