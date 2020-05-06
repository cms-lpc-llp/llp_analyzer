import ROOT as rt
import csv
import re
import sys
import collections
from collections import OrderedDict
import uproot
import numpy as np
import time
import numba
from numba import jit
from matplotlib import pyplot as plt
from numpy import linalg as LA
from ROOT import TLorentzVector
from array import array
# import PyTEX
import os

import math

sys.path.append('/storage/user/jmao/gpu/jmao/cms-llp/delayed_jet_analyzer/lib/')

from histo_utilities import create_TH1D, create_TH2D, create_TGraph,std_color_list

import CMS_lumi, tdrstyle 
tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "     Simulation Preliminary"

lumi = 137000 #in pb-1
     
print(sys.version)

'''
#labels
tags = []
tags = [
        'W(lv)Jets, background',
        'W(lv)H(bb), C1N2, 200 GeV, 150 GeV, 1 m',
        'W(lv)H(bb), C1N2, 200 GeV, 1 GeV, 1 m',
        'W(lv)H(bb), TwinHiggs, 125 GeV, 40 GeV, 1 m',
        ]


# directory
home_dir = '/mnt/hadoop/store/group/phys_exotica/jmao/susy_llp/llp_analyzer/'


# file names

fnames = {}

fnames['W(lv)H(bb), C1N2, 200 GeV, 150 GeV, 1 m'] = 'testrun/RunIIFall17_x1n2-n1-wlv-hbb_mchi200_mlsp150_pl1000_ev100000_fullsim_signal_aod_llp_analyzer_1pb_weighted.root'
fnames['W(lv)H(bb), C1N2, 200 GeV, 1 GeV, 1 m'] =  'testrun/RunIIFall17_x1n2-n1-wlv-hbb_mh200_pl1000_ev100000_fullsim_signal_aod_llp_analyzer_1pb_weighted.root'
fnames['W(lv)H(bb), TwinHiggs, 125 GeV, 40 GeV, 1 m'] = 'testrun/WH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_llp_analyzer_1pb_weighted.root' 
fnames['W(lv)Jets, background'] = 'testrun/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_llp_analyzer_1pb_weighted.root'
# fnames['W(lv)Jets, background'] = 'testrun/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6-v1_v1_v1_Job16_Of_16_1pb_weighted.root'

# output plots directory
out_dir = '/storage/user/jmao/gpu/jmao/cms-llp/plots/20191024/'

 
# colors
cols = {}

cols['W(lv)H(bb), C1N2, 200 GeV, 150 GeV, 1 m'] = 215 
cols['W(lv)H(bb), C1N2, 200 GeV, 1 GeV, 1 m'] = 207  
cols['W(lv)H(bb), TwinHiggs, 125 GeV, 40 GeV, 1 m'] = 221 
cols['W(lv)Jets, background'] = 209 


#styles
stys = {}

stys['W(lv)H(bb), C1N2, 200 GeV, 150 GeV, 1 m'] = 3 
stys['W(lv)H(bb), C1N2, 200 GeV, 1 GeV, 1 m'] = 4  
stys['W(lv)H(bb), TwinHiggs, 125 GeV, 40 GeV, 1 m'] = 2 
stys['W(lv)Jets, background'] = 1

# get sig n bkg the trees
tree = OrderedDict()

# run_sig = 'TwinHiggs'
run_sig = '150 GeV'
# run_sig = '1 GeV'

for tag in tags:

    if 'background' in tag or run_sig in tag:
        
        if 'background' not in tag:
            run_sig_key = tag
            
        print(home_dir+fnames[tag])
        
 
        root_dir = uproot.open(home_dir+fnames[tag]) 
        tree[tag] = root_dir['SusyLLPTree']

        print ('Open ready')
        print(tag)
        print ('Tree ready') 
        print(tree[tag], tree)
        v = tree[tag]
    #     llp_d_pid = v['gLLP_daughter_pid'].array()
    #     print(llp_d_pid[:10])
        weight = v['weight'].array()
        print(weight[:10])
        jet_amax = v['jetGammaMax_ET'].array()
        print(jet_amax[:10])
        if 'back' not in tag:
            calo_jet_amax = v['gLLP0_EB'].array()
            print(calo_jet_amax[:10])

trigger_names_file = '/storage/user/jmao/gpu/jmao/cms-llp/delayed_jet_analyzer/data/trigger_names_llp_v1.dat'
trigger_names = []
with open(trigger_names_file) as f:
    reader = csv.reader(f, delimiter=" ")
    for line in reader:
        trigger_names.append(line[2])

trigger_paths = [87,135]  

# 310 HLT_PFMET120_PFMHT120_IDTight
# 87 HLT_Ele32_WPTight_Gsf
# 135 HLT_IsoMu24

# weight

weight = {}

# lep muon or electron

lep = 13 #muon
# lep = 11 #ele

lep_id = {}
lep_pt = {}

# pf jet

pf_jet_pt = {}
pf_jet_time = {}
pf_jet_nrechits = {}

pf_jet_chef = {}
# pf_jet_h_over_e = {}

pf_jet_tmf_et = {}
pf_jet_delta_r = {}
# pf_jet_pt_trk = {}

# HLT
hlt_trg = {}

i = 0

for k,v in tree.items(): 
    print(k,v)
    
    weight[k] = []
    
    lep_id[k] = []
    lep_pt[k] = []
    
    pf_jet_pt[k] = []
    pf_jet_time[k] = []
    pf_jet_nrechits[k] = []
    pf_jet_chef[k] = []
#     pf_jet_h_over_e[k] = []
    pf_jet_tmf_et[k] = []
    pf_jet_delta_r[k] = []
#     pf_jet_pt_trk[k] = []
    
    hlt_trg[k] = []
 
    #branches for selection
    
    hlt = v['HLTDecision'].array()

    # select only triggered events
    sel_hlt = np.zeros(hlt[:,0].shape)
    for tr in trigger_paths:
        sel_hlt  = np.logical_or(sel_hlt, hlt[:,tr])
    
    lepid = v['lepPdgId'].array()
    leppt = v['lepPt'].array()
    
    sel_lep_id = abs(lepid==lep)
    
    lep_id_sel = lepid[sel_lep_id]
    lep_pt_sel = leppt[sel_lep_id]
    
    sel_lep = lep_pt_sel.counts >= 1
    
    pt = v['jetPt'].array()
    eta = v['jetEta'].array()
    nrechits = v['ecalNRechits'].array()

    sel_jet_eta = np.logical_and(abs(eta)<1.48, nrechits>5)

    j_pt_sel = pt[sel_jet_eta]
    
    sel_jet = j_pt_sel.counts >= 1
    
    sel_obj = np.logical_and(sel_lep, sel_jet)
    sel_evt = np.logical_and(sel_hlt, sel_obj)
    
    # variables
    
    w = v['weight'].array()

    time = v['jetTime'].array()
    
    charged_em = v['jetChargedEMEnergyFraction'].array()
    neutral_em = v['jetNeutralEMEnergyFraction'].array()
        
    charged_had = v['jetChargedHadronEnergyFraction'].array()
    neutral_had = v['jetNeutralHadronEnergyFraction'].array()
    
    tmf_et = v['jetGammaMax_ET'].array()
    delta_r = v['jetMinDeltaRPVTracks'].array()
#     pt_trk = v['jetPtAllPVTracks'].array()
        
    # selection applied to variables

    basic_str = sel_evt
 
    # selected variables        
    
    w_sel = w[basic_str]
    
    j_pt = pt[basic_str]
    
    j_time = time[basic_str]
    j_nrechits = nrechits[basic_str]
    j_charged_em = charged_em[basic_str]
    j_neutral_em = neutral_em[basic_str]
    j_charged_had = charged_had[basic_str]
    j_neutral_had = neutral_had[basic_str]
    
    j_em = [x+y for x,y in zip(j_charged_em, j_neutral_em)]
    j_had = [x+y for x,y in zip(j_charged_had, j_neutral_had)]
    
#     j_h_over_e = []
#     for index1, (x,y) in enumerate(zip(j_had, j_em)):
#         for index2,(x1,y1) in enumerate(zip(x,y)):
#             if y1>0:
#                 j_h_over_e.append(x1/y1)
    
    j_tmf_et = tmf_et[basic_str]
    j_delta_r = delta_r[basic_str]
#     j_pt_trk = pt_trk[basic_str]

    #per event

    
    if 'C1N2' in k:
        # br of w -> l v
        w_all = w_sel * lumi * 0.324
    elif 'TwinHiggs' in k:
        w_all = w_sel * lumi * br
    else:
        w_all = w_sel * lumi
    
    j_pt_all = j_pt 
    j_time_all = j_time 
    j_nrechits_all = j_nrechits
    
    j_chef_all = j_charged_had 
#     j_h_over_e_all = j_h_over_e.flatten()
#     j_h_over_e_all = np.array(j_h_over_e)
    
    j_tmf_et_all = j_tmf_et 
    j_delta_r_all = j_delta_r 
#     j_pt_trk_all = j_pt_trk 

    # print checkout
    
    print(np.unique(w_all))
    
#     print(w_all[:10])
#     print(j_pt_all[:10])
#     print(j_time_all[:10])
#     print(j_nrechits_all[:10])
    
#     print(j_chef_all[:10])
# #     print(j_h_over_e_all[:10])
    
#     print(j_tmf_et_all[:10])
#     print(j_delta_r_all[:10])
# #     print(j_pt_trk_all[:10])

    # print len
    print(len(w_all))
    
    print(len(j_pt_all))
    print(len(j_nrechits_all))
    
    print(len(j_time_all))
    print(len(j_chef_all))
    print(len(j_tmf_et_all))
    print(len(j_delta_r_all))    
    
    # print len of 5
#     print(len(j_pt_all[5]))
#     print(len(j_nrechits_all[5]))
    
 
    
    # assign 
    
#     weight[k] = np.array(w_all)
#     pf_jet_pt[k] = np.array(j_pt_all)
#     pf_jet_time[k] = np.array(j_time_all)
#     pf_jet_nrechits[k] = np.array(j_nrechits_all)
#     pf_jet_chef[k] = np.array(j_chef_all)
#     pf_jet_tmf_et[k] = np.array(j_tmf_et_all)
# #     pf_jet_h_over_e[k] = np.array(j_h_over_e_all)
#     pf_jet_delta_r[k] = np.array(j_delta_r_all)
# #     pf_jet_pt_trk[k] = np.array(j_pt_trk_all)

    weight[k] = w_all
    pf_jet_pt[k] = j_pt_all
    pf_jet_time[k] = j_time_all
    pf_jet_nrechits[k] = j_nrechits_all
    pf_jet_chef[k] = j_chef_all
    pf_jet_tmf_et[k] = j_tmf_et_all
    pf_jet_delta_r[k] = j_delta_r_all

 
    i += 1

# optimization
from skopt import gp_minimize
from skopt.space import Real, Integer
from skopt.utils import use_named_args

sig_bins = [0,1,2]



space = [Real(-5,5, name='jet_time'),
          Real(0.0,1.5, name='jet_tmf'),
          Real(0.0,1.5,name='jetMinDeltaRPVTracks'),
          Real(0.0,1.5,name = 'jet_chef')
          ]

def figure_of_merit(time_cut = -50, tmf_et_cut = 50, delta_r_cut = -50, chef_cut = -50):
    nTags = {}
    sig_count = [] 
    bkg_count = np.zeros((3,))
    
    for k,T in pf_jet_tmf_et.items():
#         start_t = time.time()
        nTags[k] = []
        tagged_jet = np.logical_and(pf_jet_tmf_et[k]<tmf_et_cut, pf_jet_delta_r[k]>delta_r_cut)
        tagged_jet = np.logical_and(tagged_jet, pf_jet_time[k]>time_cut)
        tagged_jet = np.logical_and(tagged_jet, pf_jet_chef[k]<chef_cut)
#         nTags[k] = [len(x) for x in pf_jet_pt[k][tagged_jet]] #event level variable, number of tags jet per event, no weights added
#         nTags[k] = np.array(nTags[k])
        nTags[k] = pf_jet_pt[k][tagged_jet].count() #event level variable, number of tags jet per event, no weights added
#         nTags[k] = np.sum(w[tagged_jet])
#         nJets[k] = np.sum(w[w<8000000])
#         print(k, nTags[k][:10])
        print(k, np.unique(nTags[k]))
        nTags[k] = np.array(nTags[k])
            
        occurCount = []
        for i in sig_bins:
            if i == 2:
#               #overflow bin
                tag_cut = nTags[k] >= i 
            else:
                tag_cut = nTags[k] == i
            w_sum = sum(weight[k][tag_cut])
            occurCount.append(w_sum)
        print(occurCount)
        occurCount = np.array(occurCount)
        
        if run_sig in k :
            sig_count = list(occurCount[:2])
            sig_count.append(sum(occurCount[2:]))
            sig_count = np.array(sig_count)
        elif 'background' in k:
            bkg_count = bkg_count + occurCount

#     cond = np.logical_not(np.logical_or(sig_count < 5,bkg_count < 5))
    sig_count = sig_count
    bkg_count = bkg_count
    print("sig count", sig_count)
    print("bkg count", bkg_count)
    cond = sig_count > 5
    if not np.count_nonzero(cond) == len(sig_bins):
        return 1

#     return -1.0*np.sum(2*((s+b)*math.log(s/b+1)-s))**0.5

    # loss function: - significance  = - (sqrt(S/(S+B)))
    return -1.0 * np.sum(sig_count**2/(sig_count+bkg_count))**0.5

@use_named_args(space)
def objective(**X):
    print("New configuration: {}".format(X))
#     begt = time.time()
    fom = figure_of_merit(time_cut=X['jet_time'], tmf_et_cut=X['jet_tmf'],delta_r_cut=X['jetMinDeltaRPVTracks'],chef_cut = X['jet_chef'])
#     fom = train(model, learning_rate=X['learning_rate'])
    return fom


# res_gp = gp_minimize(objective, space, n_calls=5, n_random_starts=1, random_state=123, verbose=True)
res_gp = gp_minimize(objective, space, n_calls=500, n_random_starts=100, random_state=123, verbose=True)

from skopt.plots import plot_convergence
print(run_sig_key)
 plot_convergence(res_gp)
 print("Best parameters: \
 \n jet_time_cuts = [{}] \
 \n jet_tmf_cuts = [{}] \
 \n jet_delta_r_cuts = [{}]\
  \n jet_chef_cuts = [{}]".format(res_gp.x[0],
                                  res_gp.x[1],
                                  res_gp.x[2],
                                  res_gp.x[3]))
'''
