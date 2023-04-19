import numpy as np
import sys
import math
import ROOT as rt
import uproot
sys.path.append('/storage/user/christiw/gpu/christiw/llp/delayed_jet_analyzer/lib/')
from histo_utilities import create_TH2D
def reweightXsec_tau(ctau,mass):
    #ctau in mm
    #mass in GeV
    #xsec in pb
    cof={
        "1":10.80079931781248,
        "2":7.338888133838879,
        "4":3.351656936102019,
        "7":-2.2837633387534577,
        "10":-4.251700372998307,
    }
    return np.exp(-1*np.log(ctau)+cof[mass])
def weight_calc(llp_ct, new_ctau, old_ctau):
    if new_ctau>10000 and old_ctau ==10000: #100m, the largest ctau sample
        #assert(old_ctau == 10000)
        # set the bins
        bins = [0,10,20,30,50,75,100,200,300,400,500,600,700,800,900,1000] #10m-->100
        bins = bins *2
        bins+=[int(len(bins)/2),int(len(bins)/2)]
        h_weight = create_TH2D( np.column_stack((np.array([]),np.array([]))), name='weight', binning=bins)

        for x in range(h_weight.GetNbinsX()):
            for y in range(h_weight.GetNbinsY()):
                x1,x2 = h_weight.GetXaxis().GetBinLowEdge(x+1), h_weight.GetXaxis().GetBinLowEdge(x+1)+h_weight.GetXaxis().GetBinWidth(x+1)
                y1,y2 = h_weight.GetYaxis().GetBinLowEdge(y+1), (h_weight.GetYaxis().GetBinLowEdge(y+1)+h_weight.GetYaxis().GetBinWidth(y+1))
                if x == h_weight.GetNbinsX()-1:
                    bin_content = (math.exp(-x1/new_ctau))*(math.exp(-y1/new_ctau)-math.exp(-y2/new_ctau))
                    bin_content /= (math.exp(-x1/old_ctau))*(math.exp(-y1/old_ctau)-math.exp(-y2/old_ctau))
                elif y == h_weight.GetNbinsY()-1:
                    bin_content = (math.exp(-x1/new_ctau)-math.exp(-x2/new_ctau))*(math.exp(-y1/new_ctau))
                    bin_content /= (math.exp(-x1/old_ctau)-math.exp(-x2/old_ctau))*(math.exp(-y1/old_ctau))
                else:
                    bin_content = (math.exp(-x1/new_ctau)-math.exp(-x2/new_ctau))*(math.exp(-y1/new_ctau)-math.exp(-y2/new_ctau))
                    bin_content /=  (math.exp(-x1/old_ctau)-math.exp(-x2/old_ctau))*(math.exp(-y1/old_ctau)-math.exp(-y2/old_ctau))

                h_weight.SetBinContent(x+1,y+1, bin_content)

        outFile = rt.TFile('test.root', 'RECREATE')
        outFile.WriteTObject(h_weight, 'h_weight', "WriteDelete");
        outFile.Close();

        root_dir = uproot.open('test.root')
        h_reweight = root_dir['h_weight']
        return h_reweight.values[np.argmax(h_reweight.edges[0]>llp_ct[:,0][:,None],axis=1)-1, np.argmax(h_reweight.edges[1]>llp_ct[:,1][:,None],axis=1)-1]



    else:
        if llp_ct.ndim > 1:llp_ct = np.array(np.sum(llp_ct,axis=1))
        source = np.exp(-1.0*llp_ct/old_ctau)/old_ctau**2
        weight = 1.0/new_ctau**2 * np.exp(-1.0*llp_ct/new_ctau)/source
        return weight


# def weight_calc_1llp(llp_ct, new_ctau, old_ctau):
#     llp_ct = llp_ct[:,0]
#     source = np.exp(-1.0*llp_ct/old_ctau)/old_ctau**2
#     weight = 1.0/new_ctau**2 * np.exp(-1.0*llp_ct/new_ctau)/source
#     return weight
def make_datacard_2sig(outDataCardsDir,modelName,  signal_rate, norm, bkg_rate, observation, bkg_unc, bkg_unc_name, sig_unc, sig_unc_name, prefix):
    a,b,c,d = bkg_rate[0], bkg_rate[1], bkg_rate[2], bkg_rate[3]
    c1 = a/b
    c2 = c/b
    nSig = len(signal_rate.keys())
    text_file = open(outDataCardsDir+modelName+".txt", "w")
    text_file.write('# signal norm {0} \n'.format(norm))

    text_file.write('imax {0} \n'.format(4))
    text_file.write('jmax {0} \n'.format(nSig))
    text_file.write('kmax * \n')
    text_file.write('shapes * * FAKE \n')



    text_file.write('--------------- \n')
    text_file.write('--------------- \n')
    text_file.write('bin \t chA \t chB \t chC \t chD \n')
    text_file.write('observation \t {0:6.2f} \t {1:6.2f} \t {2:6.2f} \t {3:6.2f} \n'.format(observation[0],observation[1],observation[2],observation[3]))
    text_file.write('------------------------------ \n')
    text_file.write('bin '+'\t chA ' * (1+nSig) + '\t chB ' * (1+nSig) +'\t chC '*(1+nSig) +'\t chD '*(1+nSig) +'\n')
    process_name = '\t '+ (' \t ').join(list(signal_rate.keys())) + '\t bkg '
    text_file.write('process ' + process_name * 4 + '\n')
    process_number = '\t '+ (' \t ').join(list((np.arange(nSig)*-1).astype(str))) + '\t 1'
    text_file.write('process ' + process_number * 4 + '\n')
    rate_string = 'rate'
    for i in range(4):# 4 bins
        for k,v in signal_rate.items():
            rate_string +='\t {0:e} '.format(v[i])
        rate_string += '\t 1 '
    text_file.write(rate_string+'\n')
    text_file.write('------------------------------ \n')

    text_file.write(prefix+'A	 rateParam	 chA	 bkg 	  (@0*@2/@1)			'+prefix+'B,'+prefix+'C,'+prefix+'D \n')
    text_file.write(prefix+'B	 rateParam	 chB	 bkg 	 {0:.2f}        [0,{1:.2f}] \n'.format(b, b*7)) 
    text_file.write(prefix+'C	 rateParam	 chC	 bkg 	 {0:.2f}        [0,{1:.2f}] \n'.format(c, c*7)) 
    text_file.write(prefix+'D	 rateParam	 chD	 bkg 	 {0:.2f}        [0,{1:.2f}] \n'.format(d, d*7)) 
    for k,v in signal_rate.items():
        text_file.write('norm rateParam * {0} 1  \n'.format(k))


  #### uncertainties ####
    for k,v in sig_unc.items():assert(len(sig_unc_name)==len(v))
    for i in range(len(sig_unc_name)):
        if 'mc_stats' in sig_unc_name[i]:
            for j, bin in enumerate(['A', 'B', 'C', 'D']):#bin
                    for l, k in enumerate(sig_unc.keys()): #channels
                        before = (len(sig_unc.keys())+1)*j+l
                        after = (len(sig_unc.keys())+1)*4-before-1
                        if sig_unc[k][i][j] > 0.0: text_file.write(sig_unc_name[i]+'_'+k+'_'+bin+' \t gmN ' +str(int(sig_unc[k][i][j]))+ '  '+'\t -  '*before + str(signal_rate[k][j]/int(sig_unc[k][i][j])) + '\t - '*after +'\n')

        else:

            unc_text = sig_unc_name[i]+' \t lnN'
            if len(sig_unc[list(sig_unc.keys())[0]][i])==4:#symmetric uncertainties
                for j in range(4):#bin
                    for k,v in sig_unc.items():
                        if v[i][j] == 0.0:unc_text += ' \t -'
                        else: unc_text += ' \t '+str(v[i][j]+1)
                    unc_text += '\t - '
            else:#asymmetric
                for j in range(4):#bin A, B, C, D
                    for k,v in sig_unc.items():
                        if  v[i][j] == 0.0 and v[i][j+4] == 0.0: unc_text += ' \t -'
                        else:unc_text += ' \t {0}/{1}'.format(1-v[i][j],1+v[i][j+4])
                    unc_text += '\t -'
            text_file.write(unc_text + ' \n')
    for i in range(len(bkg_unc_name)):
        bkg_unc_text = bkg_unc_name[i] + ' \t lnN ' + '\t - '*(4*nSig+3) + '\t ' + str(1+bkg_unc[i]) + ' \n'
        text_file.write(bkg_unc_text)

    text_file.close()

def make_datacard_2tag(outDataCardsDir,modelName,  signal_rate, norm, bkg_rate, observation, bkg_unc, bkg_unc_name, sig_unc, sig_unc_name,signal_region, prefix):
    a,b,c,d = bkg_rate[0], bkg_rate[1], bkg_rate[2], bkg_rate[3]
    nSig = len(signal_rate.keys())
    text_file = open(outDataCardsDir+modelName+".txt", "w")
    text_file.write('# signal norm {0} \n'.format(norm))

    text_file.write('imax {0} \n'.format(4))
    text_file.write('jmax {0} \n'.format(nSig))
    text_file.write('kmax * \n')
    text_file.write('shapes * * FAKE \n')



    text_file.write('--------------- \n')
    text_file.write('--------------- \n')
    text_file.write('bin \t chA \t chB \t chC \t chD \n')
    text_file.write('observation \t {0:6.2f} \t {1:6.2f} \t {2:6.2f} \t {3:6.2f} \n'.format(observation[0],observation[1],observation[2],observation[3]))
    text_file.write('------------------------------ \n')
    text_file.write('bin '+'\t chA ' * (1+nSig) + '\t chB ' * (1+nSig) +'\t chC '*(1+nSig) +'\t chD '*(1+nSig) +'\n')
    process_name = '\t '+ (' \t ').join(list(signal_rate.keys())) + '\t bkg '
    text_file.write('process ' + process_name * 4 + '\n')
    process_number = '\t '+ (' \t ').join(list((np.arange(nSig)*-1).astype(str))) + '\t 1'
    text_file.write('process ' + process_number * 4 + '\n')
    rate_string = 'rate'
    for i in range(4):# 4 bins
        for k,v in signal_rate.items():
            rate_string +='\t {0:e} '.format(v[i])
        rate_string += '\t 1 '
    text_file.write(rate_string+'\n')
    text_file.write('------------------------------ \n')

    text_file.write(prefix+'A   rateParam       chA     bkg      (@0*@2/@1)                    '+prefix+'B,'+prefix+'C,'+prefix+'D \n')
    if b == 0: text_file.write(prefix+'B   rateParam       chB     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(b, c*7))
    else: text_file.write(prefix+'B   rateParam       chB     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(b, b*7))
    text_file.write(prefix+'C   rateParam       chC     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(c, c*7))
    if d == 0:text_file.write(prefix+'D   rateParam       chD     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(d, c*7))
    else: text_file.write(prefix+'D   rateParam       chD     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(d, d*7))


    for k,v in signal_rate.items():
        text_file.write('norm rateParam * {0} 1  \n'.format(k))


  #### uncertainties ####
    for k,v in sig_unc.items():assert(len(sig_unc_name)==len(v))
    for i in range(len(sig_unc_name)):
        if 'mc_stats' in sig_unc_name[i]:
            for j, bin in enumerate(['A', 'B', 'C', 'D']):#bin
                    for l, k in enumerate(sig_unc.keys()): #channels
                        before = (len(sig_unc.keys())+1)*j+l
                        after = (len(sig_unc.keys())+1)*4-before-1
                        if sig_unc[k][i][j] > 0.0: text_file.write(sig_unc_name[i]+'_'+k+'_'+bin+' \t gmN ' +str(int(sig_unc[k][i][j]))+ '  '+'\t -  '*before + str(signal_rate[k][j]/int(sig_unc[k][i][j])) + '\t - '*after +'\n')

        else:

            unc_text = sig_unc_name[i]+' \t lnN'
            if len(sig_unc[list(sig_unc.keys())[0]][i])==4:#symmetric uncertainties
                for j in range(4):#bin
                    for k,v in sig_unc.items():
                        if v[i][j] == 0.0:unc_text += ' \t -'
                        else: unc_text += ' \t '+str(v[i][j]+1)
                    unc_text += '\t - '
            else:#asymmetric
                for j in range(4):#bin A, B, C, D
                    for k,v in sig_unc.items():
                        if  v[i][j] == 0.0 and v[i][j+4] == 0.0: unc_text += ' \t -'
                        else:unc_text += ' \t {0}/{1}'.format(1-v[i][j],1+v[i][j+4])
                    unc_text += '\t -'
            text_file.write(unc_text + ' \n')
    for i in range(len(bkg_unc_name)):
        bkg_unc_text = bkg_unc_name[i] + ' \t lnN ' + '\t - '*(4*nSig+3) + '\t ' + str(1+bkg_unc[i]) + ' \n'
        text_file.write(bkg_unc_text)

    text_file.close()



def make_datacard_2tag_3bins(outDataCardsDir,modelName,  signal_rate, norm, bkg_rate, observation, bkg_unc, bkg_unc_name, sig_unc, sig_unc_name,signal_region, prefix):
    #only works if expected background is the same in B and D, so the selection is the same for both clusters
    nBins = 3

    assert(len(bkg_rate) == nBins)
    assert(len(observation) == nBins)

    a,b,c = bkg_rate[0], bkg_rate[1], bkg_rate[2]

    nSig = len(signal_rate.keys()) # number of production modes
    text_file = open(outDataCardsDir+modelName+".txt", "w")
    text_file.write('# signal norm {0} \n'.format(norm))

    text_file.write('imax {0} \n'.format(nBins))# 3 bins
    text_file.write('jmax {0} \n'.format(nSig))
    text_file.write('kmax * \n')
    text_file.write('shapes * * FAKE \n')

    text_file.write('--------------- \n')
    text_file.write('--------------- \n')
    text_file.write('bin \t chA \t chB \t chC \n')
    text_file.write('observation \t {0:6.2f} \t {1:6.2f} \t {2:6.2f} \n'.format(observation[0],observation[1],observation[2]))
    text_file.write('------------------------------ \n')
    text_file.write('bin '+'\t chA ' * (1+nSig) + '\t chB ' * (1+nSig) +'\t chC '*(1+nSig)  +'\n')
    process_name = '\t '+ (' \t ').join(list(signal_rate.keys())) + '\t bkg '
    text_file.write('process ' + process_name * nBins + '\n')
    process_number = '\t '+ (' \t ').join(list((np.arange(nSig)*-1).astype(str))) + '\t 1'
    text_file.write('process ' + process_number * nBins + '\n')
    rate_string = 'rate'
    for i in range(nBins):# 3 bins
        for k,v in signal_rate.items():
            rate_string +='\t {0:e} '.format(v[i])
        rate_string += '\t 1 '
    text_file.write(rate_string+'\n')
    text_file.write('------------------------------ \n')


    text_file.write(prefix+'A	 rateParam	 chA	 bkg 	  (@0/(2*@1))**2*@1			'+prefix+'B,'+prefix+'C \n')
    if b == 0: text_file.write(prefix+'B   rateParam       chB     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(b, b*7+c*7))
    else: text_file.write(prefix+'B   rateParam       chB     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(b, b*7))
    text_file.write(prefix+'C   rateParam       chC     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(c, c*7))

    for k,v in signal_rate.items():
        text_file.write('norm rateParam * {0} 1  \n'.format(k))


  #### uncertainties ####
    for k,v in sig_unc.items():assert(len(sig_unc_name)==len(v))
    for i in range(len(sig_unc_name)): # loop over different signal efficiencies
        if 'mc_stats' in sig_unc_name[i]:
            for j, bin in enumerate(['A', 'B', 'C']):#bin
                    for l, k in enumerate(sig_unc.keys()): #channels
                        before = (len(sig_unc.keys())+1)*j+l
                        after = (len(sig_unc.keys())+1)*nBins-before-1
                        if sig_unc[k][i][j] > 0.0: text_file.write(sig_unc_name[i]+'_'+k+'_'+bin+' \t gmN ' +str(int(sig_unc[k][i][j]))+ '  '+'\t -  '*before + str(signal_rate[k][j]/int(sig_unc[k][i][j])) + '\t - '*after +'\n')

        else:
            unc_text = sig_unc_name[i]+' \t lnN'
            if len(sig_unc[list(sig_unc.keys())[0]][i])==nBins:#symmetric uncertainties
                for j in range(nBins):#bin
                    for k,v in sig_unc.items():
                        if v[i][j] == 0.0:unc_text += ' \t -'
                        else: unc_text += ' \t '+str(v[i][j]+1)
                    unc_text += '\t - '
            else:#asymmetric
                for j in range(nBins):#bin A, B, C
                    for k,v in sig_unc.items():
                        if  v[i][j] == 0.0 and v[i][j+nBins] == 0.0: unc_text += ' \t -'
                        else:unc_text += ' \t {0}/{1}'.format(1-v[i][j],1+v[i][j+nBins])
                    unc_text += '\t -'
            text_file.write(unc_text + ' \n')
    for i in range(len(bkg_unc_name)):
        bkg_unc_text = bkg_unc_name[i] + ' \t lnN ' + '\t - '*(nBins*nSig+nBins-1) + '\t ' + str(1+bkg_unc[i]) + ' \n'
        text_file.write(bkg_unc_text)

    text_file.close()
