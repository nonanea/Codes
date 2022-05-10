import ROOT
import os
import time
from math import sqrt,pow,fabs
from array import array

op = time.time()

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat("")

c1 = ROOT.TCanvas("c1", "canvas", 1000, 800)

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.05)
c1.SetTopMargin(0.05)
c1.SetBottomMargin(0.15)

#   Reader initialization

colorList = [
    ROOT.kGreen,
    ROOT.kBlue,
    ROOT.kTeal,
    ROOT.kYellow,
    ROOT.kGray,
    ROOT.kCyan,

    ROOT.kViolet,
    ROOT.kRed,
    ROOT.kOrange,
    ROOT.kGreen + 2,
    ROOT.kAzure,
    ROOT.kPink + 4,
    
    ROOT.kGreen + 1,
    ROOT.kAzure + 1,
    ROOT.kViolet + 1,
    ROOT.kOrange + 1,
    ROOT.kCyan + 1,
    ROOT.kYellow + 1,
    ROOT.kTeal + 1,
    ROOT.kMagenta + 1,
    ROOT.kRed + 1,
    ROOT.kBlue + 1,
    ROOT.kGreen + 2,
    ROOT.kAzure + 2,
    ROOT.kViolet + 2,
    ROOT.kOrange + 2,
    ROOT.kCyan + 2,
    ROOT.kYellow + 2,
    ROOT.kTeal + 2,
    ROOT.kMagenta + 2,
    ROOT.kRed + 2,
    ROOT.kBlue + 2
]

f = open("Plots/results.txt","a")
f.seek(0)
f.truncate()

f1 = open("Plots/cutflow.txt","a")
f1.seek(0)
f1.truncate()

lep = [ROOT.TLorentzVector(1,1,0,5) for i in range(3)]

lep_id = [array('f',[0]) for i in range(3)]
lep_pt = [array('f',[0]) for i in range(3)]
lep_eta = [array('f',[0]) for i in range(3)]
lep_phi = [array('f',[0]) for i in range(3)]
lep_E = [array('f',[0]) for i in range(3)]

met_met = array('f',[0.])
met_px = array('f',[0.])
met_py = array('f',[0.])

mll_1 = array('f',[0.])
mll_2 = array('f',[0.])

weight = array('f',[0.])
xsec = array('f',[0.])

WZ_list = [364253,364284,363358]
ZZ_list = [345666,345723,364285,364250,364283,345283,245706,363356,364254]
WW_list = [345718]
VVV_list = [364242,364253,364244,364245,364246,364247,364248,364249]
Zjets_list = [361107,361666,361667]
t_list = [410644,410645,410658,410659]
Wt_list = [410646,410647]
tt_list = [410472]
ttV_list = [410155,410156,410157,410081,410218,410219]
Zgamma_list = [364505,364506,364507,364508,364509]
# extra_list = [345705,345706,345714,364283,364284,364286,364287,364288,364289,364290,361108]
# extra_list = [364505,364506,364507,364508,364509]

plot_name = [
            "PT1", #0
            "PT2", #1
            "PT3", #2
            "MET", #3
            "M_ll_1", #4
            "M_ll_2", #5
            "N_jets", #6
            "N_bjets", #7
            "mT_met_l1", #8
            "mT_met_l2", #9
            "mT_met_l3", #10
            "mT_sum", #11
            "dR_1", #12
            "dR_2", #13
            "Delta_Phi", #14
            "M_ll_Z1", #15
            "M_ll_Z2", #16
            "METS", #17
            "dPhi_1", #18
            "dPhi_2", #19
            "VT", #20
            "HT", #21
            "LT", #22
            ]

h_name = [
            "P_{T,1} [GeV]", #0
            "P_{T,2} [GeV]", #1
            "P_{T,3} [GeV]", #2
            "MET [GeV]", #3
            "M_{ll,1} [GeV]", #4
            "M_{ll,2} [GeV]", #5
            # "M_{Z'} [GeV]", #4
            # "M_{W} [GeV]", #5
            "N_{jets}", #6
            "N_{bjets}", #7
            # "m_{T,met&l1} [GeV]", #8
            # "m_{T,met&l2} [GeV]", #9
            "m_{T,met&l} [GeV]", #8
            "m_{T,met&l_{W}} [GeV]", #9
            "m_{T,met&l3} [GeV]", #10
            "m_{T,sum} [GeV]", #11
            "dR_{ll,1}", #12
            "dR_{ll,2}", #13
            # "#Delta#Phi_{3l&MET}", #12
            # "weight", #14
            "DNN", #14
            # "M_{ll,Z1} [GeV]", #15
            "M_{lll} [GeV]", #15
            "M_{ll,Z2} [GeV]", #16
            "METS", #17
            "dPhi_{1}", #18
            "dPhi_{2}", #19
            "HVT [GeV]", #20
            "HT [GeV]", #21
            "LT [GeV]", #22
            ]

xrange = [
    (0,120), #0
    (0,80), #1
    (0,40), #2
    (0,80), #3
    (0,100), #4
    (0,10), #5
    # (0,300), #0
    # (0,200), #1
    # (0,100), #2
    # (0,300), #3
    # (0,600), #4
    # (0,600), #5
    (0,7), #6
    (0,3), #7
    (0,100), #8
    (0,80), #9
    (20,200), #10
    (40,200), #11
    # (0,300), #8
    # (0,200), #9
    # (0,100), #10
    # (0,800), #11
    (0,7), #12
    (0,7), #13
    # (-0.5,0.5), #14
    (0.,1.), #14
    (0,200), #15
    (0,100), #16
    (0,10), #17
    (-7,7), #18
    (-7,7), #19
    (0,70), #20
    (40,200), #21
    (0,40), #22
]

nbin = [
    20, #0
    20, #1
    20, #2
    20, #3
    200, #4
    50, #5
    8, #6
    4, #7
    20, #8
    20, #9
    20, #10
    20, #11
    20, #12
    20, #13
    20, #14
    20, #15
    20, #16
    20, #17
    20, #18
    20, #19
    20, #20
    20, #21
    20, #22
]

process = [
    "WZ",
    "ZZ",
    # "WW",
    "Z+jets",
    # "top(s)",
    # "t",
    # "Wt",
    "tt",
    "ttV",
    # "VVV",
    # "Offshell H",
    "Zgamma",
    # "Others",
    # "Total background", # N
    # "Loose",
    # "LowPt",
    # "Medium",
    "5 GeV", # N
    "9 GeV",
    "15 GeV",
    # "19 GeV",
    # "23 GeV",
    # "27 GeV",
    # "31 GeV",
    # "35 GeV",
    "39 GeV",
    # "45 GeV",
    # "51 GeV",
    # "54 GeV",
    # "60 GeV",
    # "66 GeV",
    "69 GeV",
    "75 GeV",
    # "200 GeV",

    # "200 GeV-g=0.1",
    # "300 GeV-g=0.1",
    # "400 GeV-g=0.1",
    # "500 GeV-g=0.1",
    # "200 GeV-g=10e-4",

    # "95 GeV",
    # "150 GeV",
    # "250 GeV",
    # "350 GeV",
    # "500 GeV",
    
    "Data" # N + N_sig
]

N_sig = 6
N = len(process) - N_sig - 1

amp = 400.0
# lumi = 36.1
lumi = 139

# print N

hs = [ROOT.THStack() for i in range(len(plot_name))]
# ratio = [ROOT.TH1F() for i in range(len(plot_name))]
h_bkg = [ROOT.TH1F() for i in range(len(plot_name))]
h = [[ROOT.TH1F() for j in range(len(plot_name))] for i in range(len(process))]
ratio = [[ROOT.TH1F() for j in range(len(plot_name))] for i in range(N_sig)]

# print h 

for i in range(len(plot_name)):
    hs[i] = ROOT.THStack(plot_name[i],plot_name[i])
    h_bkg[i] = ROOT.TH1F(plot_name[i],plot_name[i], nbin[i], xrange[i][0], xrange[i][1])
    
    for j in range(N+N_sig+1):
        h[j][i] = ROOT.TH1F(process[j] + plot_name[i], process[j] + plot_name[i], nbin[i], xrange[i][0], xrange[i][1])
        if j < N:
            h[j][i].SetFillColor(int(colorList[j]))
            # h[j][i].SetLineColor(int(colorList[j]))
            # h[j][i].SetLineWidth(2)

        elif j != N+N_sig:
            ratio[j-N][i] = ROOT.TH1F(process[j] + plot_name[i], process[j] + plot_name[i], nbin[i], xrange[i][0], xrange[i][1])
            ratio[j-N][i].SetLineColor(int(colorList[j]))
            # ratio[j-N][i].SetMarkerColor(int(colorList[N+j]))
            ratio[j-N][i].SetLineWidth(2)
            # ratio[j-N][i].SetMarkerSize(2)

            h[j][i].SetLineColor(int(colorList[j]))
            h[j][i].SetLineWidth(2)
        # h[j][i].SetLineColor(ROOT.kRed)

        else:
            h[j][i].SetLineColor(ROOT.kBlack)
            h[j][i].SetLineWidth(2)

# len(cut) = 5
cut = [
    "Total",
    "Isolation tight",
    "Lep pt",
    "Njets",
    "Mll",
]

Events = [array('f',(N+N_sig+1)*[0.]) for i in range(len(cut))]
Err = [array('f',(N+N_sig+1)*[0.]) for i in range(len(cut))]
Raw = [array('i',(N+N_sig+1)*[0]) for i in range(len(cut))]

MVA_flag = 1
input_dir = "."
if MVA_flag: input_dir = "MVA/39GeV/"

# if not MVA_flag:   sample_list  = os.listdir(".")
# else:   sample_list = os.listdir("MVA/")
sample_list = os.listdir(input_dir)

# print(sample_list)
for sample in sample_list:
    if "root" not in sample:    continue
        
    # if "364254" in sample or "364285" in sample or "364250" in sample:
    # if "364" in sample or "410" in sample:    continue
    # if "005" not in sample:    continue
    if "data" in sample:    continue

    # if "364253_mc" not in sample and "mass" not in sample:    continue
    # if "_mc16" not in sample or int(sample[:6]) not in Zjets_list:   continue
    
    # if "364253" not in sample and "01_mass_500" not in sample:
    #     continue
    # if "_mc16a" in sample or "_mc16d" in sample or "_mc16e" in sample:    continue
    if "_mc16" in sample:
        # print int(sample[:6])
        if int(sample[:6]) in WZ_list:    k = 0
        elif int(sample[:6]) in ZZ_list:    k = 1
        # elif int(sample[:6]) in WW_list:    k = 2
        elif int(sample[:6]) in Zjets_list:    k = 2
        # elif int(sample[:6]) in t_list:     k = 5
        # elif int(sample[:6]) in Wt_list:     k = 6
        elif int(sample[:6]) in tt_list:    k = 3
        elif int(sample[:6]) in ttV_list:    k = 4
        elif int(sample[:6]) in Zgamma_list:    k = 5
        elif int(sample[:6]) in extra_list:    k = 6
        else:   continue

    # if "Z+jets" in sample:     k = 4
    # elif "361" in sample:
    # elif "top" in sample:     k = 5
    # elif "_WZ_" in sample:
    # elif "364253_16a.root" == sample:     k = 0
    # elif "ZZ" in sample:     k = 1
    # elif "WW" in sample:     k = 2
    # elif "_VVV_" in sample:     k = 3
    # elif "_ttV_" in sample:     k = 6
    # elif "3641" in sample:     k = 2
    elif "Zp005" in sample:    k = N
    elif "Zp009" in sample:     k = N + 1
    elif "Zp035" in sample:    k = N + 2
    elif "Zp039" in sample:    k = N + 3
    # elif "Zp045" in sample:    k = N + 3
    # elif "Zp051" in sample:    k = N + 4
    # elif "Zp054" in sample:     k = N + 4
    # elif "Zp060" in sample:     k = N + 5
    # elif "Zp066" in sample:     k = N + 4
    elif "Zp069" in sample:     k = N + 4
    elif "Zp075" in sample:    k = N + 5
    # elif "350.root" in sample:    k = N + N_sig - 1
    elif "data" in sample:    k = N + N_sig
    elif MVA_flag:    k = 0
    else:    continue
    
    sample = input_dir + sample
    # if "Loose" in sample:   k = N
    # elif "LowPt" in sample:   k = N + 1
    # elif "Med" in sample:   k = N + 2
    # else: continue
    # if MVA_flag:
    #     if "bkg" in sample: k = 0

    events = array('f',len(cut)*[0.])
    err = array('f',len(cut)*[0.])
    raw = array('i',len(cut)*[0])

    print k,sample

    # sample = "NoCut/" + sample
    f_in = ROOT.TFile(sample)
    t = f_in.Get("nominal")

    for i in range(t.GetEntries()):
        # print "Processing:{}%".format(round((i + 1) * 100 / t.GetEntries(),2)), "\r"
        
        # if i > 50:
        #     break

        t.GetEntry(i)

        if MVA_flag:
            if "bkg" in sample: k = t.bkg
            if "sig" in sample:
                # if t.mass == 5:    k = N
                # elif "Zp019" in sample:     k = N + 1
                # elif t.mass == 9:    k = N + 1
                # elif t.mass == 35:    k = N + 2
                if t.mass == 39:    k = N + 3
                # elif t.mass == 69:    k = N + 4
                # elif t.mass == 75:    k = N + 5
                else:    continue

        if k == N+N_sig:
            # weight = t.weight*lumi
            weight = 1.0
        # elif k >= N:
        #     weight = t.weight_g*lumi
        else:
            weight = t.weight*lumi
            # weight = 0.
        # print weight
        
        # if t.weight > 0.08:
        #     continue


        events[0] += weight
        err[0] += weight*weight
        raw[0] += 1

        Events[0][k] += weight
        Err[0][k] += weight*weight
        Raw[0][k] += 1

        # lep[0].SetPtEtaPhiE(t.lep0_pt,t.lep0_eta,t.lep0_phi,t.lep0_E)
        # lep[1].SetPtEtaPhiE(t.lep1_pt,t.lep1_eta,t.lep1_phi,t.lep1_E)
        # lep[2].SetPtEtaPhiE(t.lep2_pt,t.lep2_eta,t.lep2_phi,t.lep2_E)

        # if fabs(lep[0].Eta()) > 2.5 or fabs(lep[1].Eta()) > 2.5 or fabs(lep[2].Eta()) > 2.5 or t.lep2_pt < 3:
        #     continue

        # if abs(t.lep0_id) != 13 or abs(t.lep1_id) != 13 or abs(t.lep2_id) != 13:
        #     continue

        # if abs(t.lep0_id) + abs(t.lep1_id) + abs(t.lep2_id) != 37:
        #     continue
            
        # if abs(t.lep0_id + t.lep1_id + t.lep2_id) != 11:
        #     continue

        # if t.lep2_isoTight != 0 or t.lep1_isoTight != 0:
        #     continue
        # if t.lep2_isoTight != 0:
        #     continue
        events[1] += weight
        err[1] += weight*weight
        raw[1] += 1

        Events[1][k] += weight
        Err[1][k] += weight*weight
        Raw[1][k] += 1


        if t.lep1_pt < 10:    continue
        if t.lep2_pt < 6:    continue
        events[2] += weight
        err[2] += weight*weight
        raw[2] += 1

        Events[2][k] += weight
        Err[2][k] += weight*weight
        Raw[2][k] += 1


        if t.n_jets > 0:    continue
        # if t.n_bjets > 0:
        #     continue
        events[3] += weight
        err[3] += weight*weight
        raw[3] += 1

        Events[3][k] += weight
        Err[3][k] += weight*weight
        Raw[3][k] += 1

        # if t.mll_2 < 85 or t.mll_1 > 100 :
        if t.mll_1 > 80:    continue
        # if t.mll_1 < 80 or t.mll_1 > 100:     continue
        # if t.met_met < 30:    continue
        # if t.LT > 8:
        #     continue
        # if MVA_flag:
        # #     # if t.SRBDT < -0.2: continue
        #     if t.DNN_5 < 0.7: continue

        events[4] += weight
        err[4] += weight*weight
        raw[4] += 1

        Events[4][k] += weight
        Err[4][k] += weight*weight
        Raw[4][k] += 1
        

        # if t.weight < -0.1:
        #     print t.weight

        # mT_vl1 = pow(lep[0].Et() + t.met_met, 2) - pow( lep[0].Px() + t.met_px , 2) - pow ( lep[0].Py() + t.met_py , 2)
        # mT_vl2 = pow(lep[1].Et() + t.met_met, 2) - pow( lep[1].Px() + t.met_px , 2) - pow ( lep[1].Py() + t.met_py , 2)
        # mT_vl3 = pow(lep[2].Et() + t.met_met, 2) - pow( lep[2].Px() + t.met_px , 2) - pow ( lep[2].Py() + t.met_py , 2)

        # mT = pow(sqrt(lep[0].Pt()*lep[0].Pt() + 0.105*0.105) + sqrt(lep[1].Pt()*lep[1].Pt() + 0.105*0.105) + sqrt(lep[2].Pt()*lep[2].Pt() + 0.105*0.105) + t.met_met, 2) - pow( lep[0].Px() + lep[1].Px() + lep[2].Px() + t.met_px , 2) - pow ( lep[0].Py() + lep[1].Py() + lep[2].Py() + t.met_py , 2)

        # mT = pow(lep[0].Et() + lep[1].Et() + lep[2].Et() + t.met_met, 2) - pow( lep[0].Px() + lep[1].Px() + lep[2].Px() + t.met_px , 2) - pow ( lep[0].Py() + lep[1].Py() + lep[2].Py() + t.met_py , 2)

        # if mT_vl1 < 0:  mT_vl1 = 0.0
        # if mT_vl2 < 0:  mT_vl2 = 0.0
        # if mT_vl3 < 0:  mT_vl3 = 0.0
        # if mT < 0:
        #     mT = 0.0
        
        # mT_vl1 = sqrt(mT_vl1)
        # mT_vl2 = sqrt(mT_vl2)
        # mT_vl3 = sqrt(mT_vl3)
        # mT_vl1 = 0.0
        # mT_vl2 = 0.0
        # mT_vl3 = 0.0
        # mT = sqrt(mT)

        # delta_phi = ((lep[0]+lep[1]+lep[2]).Px()*t.met_px + (lep[0]+lep[1]+lep[2]).Py()*t.met_py)/sqrt((lep[0]+lep[1]+lep[2]).Pt()*(lep[0]+lep[1]+lep[2]).Pt()*t.met_met*t.met_met)

        variables = [
            # lep_pt[0][0], #0
            # lep_pt[1][0], #1
            # lep_pt[2][0], #2
            # met_met[0], #3
            # mll_1[0], #4
            # mll_2[0], #5
            t.lep0_pt, #0
            t.lep1_pt, #1
            t.lep2_pt, #2
            t.met_met, #3
            t.mll_1, #4
            t.mll_2, #5
            # t.m_Zp, #4
            # t.m_W, #5
            t.n_jets, #6
            t.n_bjets, #7
            # mT_vl1, #8
            # mT_vl2, #9
            t.mT_vl, #8
            t.mT_Wvl, #9
            # mT_vl3, #10
            # t.mT_WZ, #10
            t.mT_Wvl, #9
            t.mT, #11
            t.dR_1, #12
            t.dR_2, #13
            # 1.0,
            # 1.0,
            # 1.0,
            # t.weight, #14
            t.DNN_39, #14
            # t.mll_Z1, #15
            t.mlll, #15
            # t.mll_Z2, #16
            t.mll_1 - t.mll_2, #16
            t.met_signif, #17
            t.dPhi_1, #18
            t.dPhi_2, #19
            t.VT, #20
            t.HT, #21
            t.LT, #22
        ]

        # print weight[0],xsec[0],lumi
        # if k > N-1:
        #     print t.weight*lumi
        for j in range(len(plot_name)):
            if k != N+N_sig:
                if k < N:
                    h[k][j].Fill(variables[j],t.weight*lumi)
                    h_bkg[j].Fill(variables[j],t.weight*lumi)
                    # h[k][j].Fill(variables[j],1)
                else:
                    h[k][j].Fill(variables[j],t.weight*lumi*amp)
                # print Events[k]
            else:
                h[N+N_sig][j].Fill(variables[j],1.)

    # for j in range(len(cut)):
    #     Events[j][k] += events[j]
    #     Err[j][k] += err[j]
    #     Raw[j][k] += raw[j]

    #     print(k,t.mass,Events[j][k])

    f.write("\n"+sample)
    f.write("\n"+"Total:"+str(events[0])+"+-"+str(sqrt(err[0])))
    f.write("\n"+"Isolation tight:"+str(events[1])+"+-"+str(sqrt(err[1])))
    f.write("\n"+"Lep pt cut:"+str(events[2])+"+-"+str(sqrt(err[2])))
    f.write("\n"+"Njets cut:"+str(events[3])+"+-"+str(sqrt(err[3])))
    f.write("\n"+"Mll_1 cut:"+str(events[4])+"+-"+str(sqrt(err[4]))+"\n")

# print "Signal:",Events[N-1],"Background:",Events[N],"Data:",Events[N+1]
# flow = ["Process","Total","Isolation tight","l3 pt cut","Mll_1 cut"]
flow_bkg = ["Process"]
flow_sig = ["Mass"]
flow_bkg = flow_bkg + cut
flow_sig = flow_sig + cut

for i in range(N+N_sig+1):
    f.write("\n"+process[i])
    if i < N:   flow_bkg[0] = flow_bkg[0] + "&" + process[i]
    else:   flow_sig[0] = flow_sig[0] + "&" + process[i]

    for j in range(len(cut)):
        f.write("\n"+cut[j]+":"+str(Events[j][i])+"+-"+str(sqrt(Err[j][i]))+",Raw number:"+str(Raw[j][i]))
        if i < N:   flow_bkg[j+1] = flow_bkg[j+1] + "&" + str(round(Events[j][i],2)) + "$\\pm$" + str(round(sqrt(Err[j][i]),2))
        else:   flow_sig[j+1] = flow_sig[j+1] + "&" + str(round(Events[j][i],2)) + "$\\pm$" + str(round(sqrt(Err[j][i]),2))

    f.write("\n"+"Final:"+str(Events[4][i])+"+-"+str(sqrt(Err[4][i]))+",Raw number:"+str(Raw[4][i])+"\n")


######################  ratio calculation
rf_max = -9999
for j in range(len(plot_name)):
    rf_max = -9999
    for i in range(N+N_sig+1):
        # print process[i],Events[i]
        # f.write(process[i]+" "+str(Events[i])+"\n")
        if i < N:
            hs[j].Add(h[i][j])
        elif i != N+N_sig:
            signif = 0
            for k in range(h[i][j].GetSize()):
                # if h_bkg[j].GetBinContent(k) > 0:  signif = sqrt( signif*signif + pow( h[i][j].GetBinContent(k)/amp/sqrt(h_bkg[j].GetBinContent(k)), 2) )

                # ratio[i-N][j].SetBinContent(k, signif)
                
                # if h_bkg[j].Integral(0,k) <= 0:    ratio[i-N][j].SetBinContent(k,0)
                # else:   ratio[i-N][j].SetBinContent(k,h[i][j].Integral(0,k)/amp/sqrt(h_bkg[j].Integral(0,k)))

                # if h_bkg[j].Integral(k,h[i][j].GetSize()) <= 0:    ratio[i-N][j].SetBinContent(k,0)
                # else:   ratio[i-N][j].SetBinContent(k,h[i][j].Integral(k,h[i][j].GetSize())/amp/sqrt(h_bkg[j].Integral(k,h[i][j].GetSize())))
                bc_1 = h_bkg[j].GetBinContent(k)
                bc_2 = h[N+N_sig][j].GetBinContent(k)
                e_1 = h_bkg[j].GetBinError(k)
                e_2 = h[N+N_sig][j].GetBinError(k)
                # print(bc_1,bc_2,e_1,e_2)

                if bc_2 <= 0:    ratio[0][j].SetBinContent(k,0)
                else:   
                    ratio[0][j].SetBinContent(k,bc_1/bc_2)
                    ratio[0][j].SetBinError( k, sqrt( pow(e_1/bc_2,2) + pow(bc_1/bc_2/bc_2*e_2,2) ) )

            f.write(process[i]+":"+str(ratio[i-N][j].GetXaxis().GetBinCenter(ratio[i-N][j].GetMaximumBin()))+"\n")
            if rf_max < ratio[i-N][j].GetMaximum(): rf_max = ratio[i-N][j].GetMaximum()

    f.write(plot_name[j]+"\n")
    # print rf_max
    ratio[0][j].GetYaxis().SetRangeUser(1e-5,rf_max*1.1)
    ratio[0][j].GetYaxis().SetNdivisions(505, ROOT.kTRUE)
    ratio[0][j].GetYaxis().SetTitle ("S/#sqrt{B}")
    ratio[0][j].GetXaxis().SetLabelSize(0.1)
    ratio[0][j].GetYaxis().SetLabelSize(0.1)
    ratio[0][j].GetYaxis().SetTitleSize(0.1)
    ratio[0][j].GetYaxis().SetTitleOffset(0.3)

f.close()

for i in range(len(flow_bkg)):
    f1.write(flow_bkg[i]+"\\\\\n\\midrule\n")
for i in range(len(flow_sig)):
    f1.write(flow_sig[i]+"\\\\\n\\midrule\n")

# c1.SetLogy()
# c1.SetLogx()

c2 = ROOT.TCanvas("c2", "canvas", 1200, 800)

c2.SetLeftMargin(0.25)
c2.SetRightMargin(0.25)
c2.SetTopMargin(0.05)
# c2.SetBottomMargin(0.15)
c2.Divide(3,2,0.01,0.01)

leg_cache = []

for i in range(len(plot_name)):
    # pad = ROOT.TPad(plot_name[i],plot_name[i],0 ,0.3 ,1 ,1)
    # frame = ROOT.gPad.DrawFrame(xrange[i][0],1e-5,xrange[i][1],scale)

    pad1 = ROOT.TPad("pad1","pad1",0,0.,1,1)
    pad2 = ROOT.TPad("pad2","pad2",0,0.05,1,0.3)
    pad1.SetBottomMargin(0.1)
    pad1.SetLeftMargin(0.15)
    pad2.SetTopMargin (0.)
    pad2.SetBottomMargin(0.25)
    pad2.SetLeftMargin(0.15)

    c1.cd()
    pad1.Draw()
    pad1.cd()
    # pad1.SetLogy()
    # ROOT.gPad.Update()

    # h_ref = ROOT.TH1F("ref","ref",nbin[i],xrange[i][0],xrange[i][1])
    h_bkg[i].GetXaxis().SetTitle(h_name[i])
    h_bkg[i].GetYaxis().SetTitle("Events")
    if h[N+N_sig][i].GetMaximum() <= 0: h_bkg[i].GetYaxis().SetRangeUser(1e-5,hs[i].GetMaximum()*3.0/2)
    else:   h_bkg[i].GetYaxis().SetRangeUser(1e-5,h[N][i].GetMaximum()*200)
    # else:   h_bkg[i].GetYaxis().SetRangeUser(1e-5,h[N+N_sig][i].GetMaximum()*3.0/2)
    # h_bkg[i].GetYaxis().SetRangeUser(1e-5,h[N+1][i].GetMaximum()*3.0/2)
    # h_bkg[i].GetYaxis().SetRangeUser(1,hs[i].GetMaximum()*200)
    h_bkg[i].Draw("hist")

    pt = ROOT.TPaveText(0.33,0.92,0.23,0.75,"NDC")
    pt.SetBorderSize(0)
    pt.SetFillColor(0)
    # pt.SetFillStyle(0)
    pt.SetTextFont(42)
    pt.SetTextAlign(12)
    pt.SetTextSize(0.04)

    text = ROOT.TText()
    text = pt.AddText("#scale[1.3]{#it{ATLAS} Internal}")
    text = pt.AddText("Low mass Z")
    text = pt.AddText("13 TeV, 139 fb^{-1}")
    # text=pt.AddText(categories[type])

    # pt.Draw("same")

    leg = ROOT.TLegend(0.90,0.88,0.75,0.88-len(process)*0.03)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.03)

    # leg.AddEntry(h[0][i],process[0],"f")

    hs[i].Draw("same hist")

    for j in range(N+N_sig):
        if j < N:
            # continue
            leg.AddEntry(h[j][i],process[j],"f")
            # h[j][i].Draw("same hist")
        # leg.AddEntry(h[j+1][i],process[j+1],"f")
        else:
            leg.AddEntry(h[j][i],process[j],"f")
            # print(h_name[i],j,h_bkg[i].GetMaximum(),h[j][i].Integral())
            
            # if h[j][i].Integral() > 0: h[j][i].Scale(h_bkg[i].Integral()/h[j][i].Integral())
            h[j][i].Draw("same hist")
    
            # c1.cd()
            # pad2.Draw()
            # pad2.cd()
            # if j - N == 0: ratio[0][i].Draw("L")
            # else: ratio[j-N][i].Draw("L same")

            # c1.cd()
            # pad1.Draw()
            # pad1.cd()

        # h[j+1][i].Draw("same hist")

    h[N+N_sig][i].Draw("same E1")
    leg.AddEntry(h[N+N_sig][i],"Data","lep")
    leg.Draw("same")
    leg_cache.append(leg)

    # c1.cd()
    # pad2.Draw()
    # pad2.cd()
    # ratio[0][i].GetYaxis().SetRangeUser(0.,1.2)
    # ratio[0][i].Draw("E1")
    
    c1.SaveAs("Plots/" + plot_name[i] + ".png")
    c1.SaveAs("Plots/" + plot_name[i] + ".svg")
    # c1.SaveAs("Plots/" + plot_name[i] + ".pdf")
    c1.Clear()

    c2.cd(i+1-6*int(i/6))
    # if i+1-6*int(i/6) == 1:  h_ref.Draw("hist")
    # else: h_ref.Draw("same hist")
    h_bkg[i].Draw("hist")
    hs[i].Draw("same hist")
    for j in range(N+N_sig):
        if j >= N:
            # continue
        #     leg.AddEntry(h[j][i],process[j],"f")
        #     # h[j][i].Draw("same hist")
        # # leg.AddEntry(h[j+1][i],process[j+1],"f")
        # else:
        #     leg.AddEntry(h[j][i],process[j],"f")
            h[j][i].Draw("same hist")

    h[N+N_sig][i].Draw("same E1")
    leg.Draw("same")
    # pt.Draw("same")

    if (i+1)%6 == 0:
        c2.SaveAs("Plots/Multi/" + str((i+1)/6) +"_dist.png")
        # c2.SaveAs("Plots/Multi/" + str((i+1)/6) +"_dist.svg")
        # c2.SaveAs("Plots/Multi/" + str((i+1)/6) +"_dist.pdf")
        c2.Clear()
        c2.Divide(3,2,0.01,0.01)
        
c2.SaveAs("Plots/Multi/" + str((i+1)/6+1) +"_dist.png")
# c2.SaveAs("Plots/Multi/" + str((i+1)/6+1) +"_dist.svg")
# c2.SaveAs("Plots/Multi/" + str((i+1)/6+1) +"_dist.pdf")
c2.Clear()
# c2.Divide(3,2,0.01,0.01)

ed = time.time()
print "time consumption:",ed - op,"s"