import ROOT as rt
import os
rt.gROOT.SetBatch()

#myfile = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p8_24Apr2017_Trigger/OneLeptonFull/RazorDM_OneLeptonFull_SingleLeptonSkim_Razor2016_MoriondRereco_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root")
myfile = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p8_18Jul2017_Trigger_SingleElectron2016/OneLeptonFull/RazorDM_OneLeptonFull_SingleLeptonSkim_Razor2016_MoriondRereco_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root")

assert(myfile)
mytree = myfile.Get("ControlSampleEvent")

MR_eff = rt.TEfficiency("MR_eff",";M_{R}",60,0,800)
HT_eff = rt.TEfficiency("HT_eff",";HT",60,0,2500)
Rsq_eff = rt.TEfficiency("Rsq_eff",";R^{2}",60,0,1)
MRRsq_eff = rt.TEfficiency("MRRssq_eff",";M_{R};R^{2}",60,0,800,60,0,1.)
MR_Rsq_gt_0p5_eff = rt.TEfficiency("MR_Rsq_gt_0p03_eff",";M_{R}",60,0,800)
Rsq_MR_gt_200_eff = rt.TEfficiency("Rsq_MR_gt_200_eff",";R^{2}",60,0,1.)

MR_eff.SetStatisticOption(rt.TEfficiency.kFCP)
MR_eff.SetConfidenceLevel(0.68)

Rsq_eff.SetStatisticOption(rt.TEfficiency.kFCP)
Rsq_eff.SetConfidenceLevel(0.68)

MRRsq_eff.SetStatisticOption(rt.TEfficiency.kFCP)
MRRsq_eff.SetConfidenceLevel(0.68)

Rsq_MR_gt_200_eff.SetStatisticOption(rt.TEfficiency.kFCP)
Rsq_MR_gt_200_eff.SetConfidenceLevel(0.68)

MR_Rsq_gt_0p5_eff.SetStatisticOption(rt.TEfficiency.kFCP)
MR_Rsq_gt_0p5_eff.SetConfidenceLevel(0.68)

print "NEntries = %d" %(mytree.GetEntries())
for i in range(mytree.GetEntries()/10):
    mytree.GetEntry(i)
    if (i%10000==0): print("Get entry %d" %(i))
    #if ((mytree.HLTDecision[150] or mytree.HLTDecision[151] or mytree.HLTDecision[152] or mytree.HLTDecision[153] or mytree.HLTDecision[154])):
    if (mytree.HLTDecision[27] or mytree.HLTDecision[29] or mytree.HLTDecision[34] or mytree.HLTDecision[36] or mytree.HLTDecision[37] or mytree.HLTDecision[38] or mytree.HLTDecision[39] or mytree.HLTDecision[42]) and mytree.MR > 0 and mytree.Rsq > 0: # Ele
    #if (mytree.HLTDecision[4] or mytree.HLTDecision[13] or mytree.HLTDecision[18] or mytree.HLTDecision[20] or mytree.HLTDecision[24] or mytree.HLTDecision[29] or mytree.HLTDecision[29] or mytree.HLTDecision[34] or mytree.HLTDecision[36] or mytree.HLTDecision[37] or mytree.HLTDecision[38] or mytree.HLTDecision[39] or mytree.HLTDecision[42]) and mytree.NJets80 > 1: # 1L and 2L
    #if (mytree.HLTDecision[150] or mytree.HLTDecision[151]): #PFHT125 || PFHT200
        if mytree.MR > 200:
            passSel = mytree.HLTDecision[172] or mytree.HLTDecision[166]
            Rsq_MR_gt_200_eff.Fill(passSel, mytree.Rsq)

        if mytree.Rsq > 0.5:
            passSel = mytree.HLTDecision[172] or mytree.HLTDecision[166]
            MR_Rsq_gt_0p5_eff.Fill(passSel, mytree.MR)
            
        passSel = mytree.HLTDecision[172] or mytree.HLTDecision[166]
        MRRsq_eff.Fill(passSel, mytree.MR, mytree.Rsq)
        MR_eff.Fill(passSel, mytree.MR)
        Rsq_eff.Fill(passSel, mytree.Rsq)
        HT_eff.Fill(passSel, mytree.HT)
        
rt.gStyle.SetOptStat(0)
c1 = rt.TCanvas("c1","",600,600)

Rsq_MR_gt_200_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq_MR_gt_200.png")

MR_Rsq_gt_0p5_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_MR_Rsq_gt_0p5.png")

Rsq_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq.png")

MR_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_MR.png")

HT_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_HT.png")

MRRsq_eff.Draw("COLZ")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_MRRsq.png")
