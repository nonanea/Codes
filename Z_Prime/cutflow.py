import ROOT
import os
import time
from math import sqrt,pow,fabs
from array import array

op = time.time()

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptTitle(0)

c1 = ROOT.TCanvas("c1", "canvas", 1000, 800)

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.05)
c1.SetTopMargin(0.05)
c1.SetBottomMargin(0.15)

#   Reader initialization

colorList = [
    ROOT.kRed,
    ROOT.kBlue,
    ROOT.kGreen,
    ROOT.kAzure,
    ROOT.kYellow,
    ROOT.kViolet,
    ROOT.kOrange,
    ROOT.kCyan,
    ROOT.kTeal,
    ROOT.kMagenta,
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
]

f = open("Plots/results.txt","a")
f.seek(0)
f.truncate()

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
            "Delta_Phi", #12
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
            "mT_{met&l1} [GeV]", #8
            "mT_{met&l2} [GeV]", #9
            "mT_{met&l3} [GeV]", #10
            "mT_{sum} [GeV]", #11
            # "#Delta#Phi_{3l&MET}", #12
            "weight", #12
            ]

xrange = [
    (0,300), #0
    (0,200), #1
    (0,100), #2
    (0,300), #3
    # (0,100), #4
    # (0,100), #5
    (0,300), #4
    (0,600), #5
    (0,7), #6
    (0,3), #7
    (0,300), #8
    (0,200), #9
    (0,100), #10
    (0,400), #11
    (-0.5,0.5), #12
]

nbin = [
    20, #0
    20, #1
    20, #2
    20, #3
    40, #4
    40, #5
    8, #6
    4, #7
    20, #8
    20, #9
    20, #10
    20, #11
    50, #12
]

process = [
    # "WZ",
    # "ZZ",
    # "VV",
    # "Top(s)",
    # # "Z+jets Powheg",
    # "Z+jets",
    # "VVV",
    # "Offshell H",
    # "Others",
    # "Total background", # N
    # "5 GeV", # N
    # "17 GeV",
    # "27 GeV",
    # "31 GeV",
    # "35 GeV",
    # "39 GeV-DAI",
    # "42 GeV",
    # "45 GeV",
    # "48 GeV",
    # "51 GeV",
    # "54 GeV",
    # "57 GeV",
    # "60 GeV",

    "200 GeV-g=0.5",
    "200 GeV-g=0.1",
    # "200 GeV-g=10e-4",
    
    # "20 GeV",
    # "75 GeV",
    # "200 GeV",

    # "95 GeV",
    # "150 GeV",
    # "250 GeV",
    # "350 GeV",
    # "500 GeV",
    
    "Data" # N + N_sig
]

N_sig = 2
N = len(process) - N_sig - 1

# print N

hs = [ROOT.THStack() for i in range(len(plot_name))]
h = [[ROOT.TH1F() for j in range(len(plot_name))] for i in range(len(process))]

# print h 

for i in range(len(plot_name)):
    hs[i] = ROOT.THStack(plot_name[i],plot_name[i])
    for j in range(N+N_sig+1):
        h[j][i] = ROOT.TH1F(process[j] + plot_name[i], process[j] + plot_name[i], nbin[i], xrange[i][0], xrange[i][1])
        if j < N:
            h[j][i].SetFillColor(int(colorList[j]))
            # h[j][i].SetLineColor(int(colorList[j]))
            # h[j][i].SetLineWidth(2)

        elif j != N+N_sig:
            h[j][i].SetLineColor(int(colorList[N+j]))
            h[j][i].SetLineWidth(2)
        # h[j][i].SetLineColor(ROOT.kRed)

        else:
            h[j][i].SetLineColor(ROOT.kBlack)
            h[j][i].SetLineWidth(2)

Events = [array('f',(N+N_sig+1)*[0.]) for i in range(4)]
Err = [array('f',(N+N_sig+1)*[0.]) for i in range(4)]
Raw = [array('i',(N+N_sig+1)*[0]) for i in range(4)]

sample_list  = os.listdir(".")

for sample in sample_list:
    if "root" not in sample:
        continue
        
    # if "364254" in sample or "364285" in sample or "364250" in sample:
    #     continue
    # if "muvZp" in sample:
    #     continue
    if "data" in sample:
        continue

    # if "361" not in sample and "3641" not in sample:
    #     continue
    
    # if "200G" not in sample:
    #     continue

    # if "075" not in sample and "250.root" not in sample and "150.root" not in sample and "350.root" not in sample :
    #     continue

    # if "Z+jets" in sample:
    #     k = 2
    # # elif "361" in sample:
    # elif "top" in sample:
    #     k = 1
    # elif "_VV_" in sample:
    #     k = 0
    # elif "_VVV_" in sample:
    #     k = 3
    # elif "3641" in sample:
    #     k = 2
    # elif "muvZp005" in sample:
    #     k = N
    # elif "muvZp017" in sample:
    #     k = N + 1
    # elif "muvZp027" in sample:
    #     k = N + 2
    # elif "muvZp039" in sample:
    #     k = N
    # elif "muvZp051" in sample:
    #     k = N + 2
    # elif "muvZp060" in sample:
    #     k = N + 3
    if "200GeV_5e-1.root" in sample:
    # elif "g1-75GeV" in sample:
        k = N
    # elif "5GeV_10e-3" in sample:
    #     k = N + 1
    elif "200GeV_10e-2" in sample:
        # k = N + 1
    # elif "150.root" in sample:
    #     k = N + 1
    # elif "250.root" in sample:
    #     k = N + 2
    # elif "500.root" in sample:
        # k = N + 10
    # elif "muvZp048" in sample:
    #     k = N + 6
    # elif "muvZp051" in sample:
    #     k = N + 7
    # elif "muvZp054" in sample:
    #     k = N + 8
    # elif "muvZp057" in sample:
    #     k = N + 9
    # elif "muvZp060" in sample:
    #     k = N + 10
    # elif "350.root" in sample:
        k = N + N_sig - 1
    elif "data" in sample:
        k = N + N_sig
    else:
        # k = -1
        continue

    events = array('f',4*[0.])
    err = array('f',4*[0.])
    raw = array('i',4*[0])

    lumi = 139
    print k,sample

    f_in = ROOT.TFile(sample)
    t = f_in.Get("nominal")

    for i in range(t.GetEntries()):
        # print "Processing:{}%".format(round((i + 1) * 100 / t.GetEntries(),2)), "\r"
        
        t.GetEntry(i)

        if k != N+N_sig:
            # weight = t.weight*lumi
            weight = 1.
        else:
            weight = 0.
        
        # if t.weight > 0.08:
        #     continue


        events[0] += weight
        err[0] += weight*weight
        raw[0] += 1

        lep[0].SetPtEtaPhiE(t.lep0_pt,t.lep0_eta,t.lep0_phi,t.lep0_E)
        lep[1].SetPtEtaPhiE(t.lep1_pt,t.lep1_eta,t.lep1_phi,t.lep1_E)
        lep[2].SetPtEtaPhiE(t.lep2_pt,t.lep2_eta,t.lep2_phi,t.lep2_E)

        # if i > 50:
        #     break


        if fabs(lep[0].Eta()) > 3 or fabs(lep[1].Eta()) > 3 or fabs(lep[2].Eta()) > 3 or t.lep2_pt < 3:
            continue

        # if t.lep2_pt < 3:
        #     print t.lep2_pt
        # if t.lep0_isoTight == 0 or t.lep1_isoTight == 0 or t.lep2_isoTight == 0:
        #     continue
        # if t.lep0_isoLoose == 0 or t.lep1_isoLoose == 0 or t.lep2_isoLoose == 0:
        #     continue
        events[1] += weight
        err[1] += weight*weight
        raw[1] += 1

        # if t.lep2_pt < 7:
        #     continue
        events[2] += weight
        err[2] += weight*weight
        raw[2] += 1

        # if t.mll_2 < 85 or  t.mll_2 > 100 :
        # if t.mll_1 < 60:
        #     continue
        # if t.n_jets > 0:
        #     continue
        events[3] += weight
        err[3] += weight*weight
        raw[3] += 1

        # if t.weight < -0.1:
        #     print t.weight

        mT_vl1 = pow(lep[0].Et() + t.met_met, 2) - pow( lep[0].Px() + t.met_px , 2) - pow ( lep[0].Py() + t.met_py , 2)
        mT_vl2 = pow(lep[1].Et() + t.met_met, 2) - pow( lep[1].Px() + t.met_px , 2) - pow ( lep[1].Py() + t.met_py , 2)
        mT_vl3 = pow(lep[2].Et() + t.met_met, 2) - pow( lep[2].Px() + t.met_px , 2) - pow ( lep[2].Py() + t.met_py , 2)

        # mT = pow(sqrt(lep[0].Pt()*lep[0].Pt() + 0.105*0.105) + sqrt(lep[1].Pt()*lep[1].Pt() + 0.105*0.105) + sqrt(lep[2].Pt()*lep[2].Pt() + 0.105*0.105) + t.met_met, 2) - pow( lep[0].Px() + lep[1].Px() + lep[2].Px() + t.met_px , 2) - pow ( lep[0].Py() + lep[1].Py() + lep[2].Py() + t.met_py , 2)

        # mT = pow(lep[0].Et() + lep[1].Et() + lep[2].Et() + t.met_met, 2) - pow( lep[0].Px() + lep[1].Px() + lep[2].Px() + t.met_px , 2) - pow ( lep[0].Py() + lep[1].Py() + lep[2].Py() + t.met_py , 2)

        if mT_vl1 < 0:
            mT_vl1 = 0.0
        if mT_vl2 < 0:
            mT_vl2 = 0.0
        if mT_vl3 < 0:
            mT_vl3 = 0.0
        # if mT < 0:
        #     mT = 0.0
        
        mT_vl1 = sqrt(mT_vl1)
        mT_vl2 = sqrt(mT_vl2)
        mT_vl3 = sqrt(mT_vl3)
        # mT_vl1 = 0.0
        # mT_vl2 = 0.0
        # mT_vl3 = 0.0
        # mT = sqrt(mT)

        delta_phi = ((lep[0]+lep[1]+lep[2]).Px()*t.met_px + (lep[0]+lep[1]+lep[2]).Py()*t.met_py)/sqrt((lep[0]+lep[1]+lep[2]).Pt()*(lep[0]+lep[1]+lep[2]).Pt()*t.met_met*t.met_met)

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
            mT_vl1, #8
            mT_vl2, #9
            mT_vl3, #10
            t.mT, #11
            t.weight, #11
        ]

        # print weight[0],xsec[0],lumi
        # if k > N-1:
        #     print t.weight*lumi
        for j in range(len(plot_name)):
            if k != N+N_sig:
                if k < N:
                    h[k][j].Fill(variables[j],t.weight*lumi)
                    # h[k][j].Fill(variables[j],1)
                else:
                    h[k][j].Fill(variables[j],t.weight*lumi*100)
                # print Events[k]
            else:
                h[N+N_sig][j].Fill(variables[j],1.)

    for j in range(4):
        Events[j][k] += events[j]
        Err[j][k] += err[j]
        Raw[j][k] += raw[j]

    f.write("\n"+sample)
    f.write("\n"+"Total:"+str(events[0])+"+-"+str(sqrt(err[0])))
    f.write("\n"+"Isolation tight:"+str(events[1])+"+-"+str(sqrt(err[1])))
    f.write("\n"+"l3 pt cut:"+str(events[2])+"+-"+str(sqrt(err[2])))
    f.write("\n"+"Mll_1 cut:"+str(events[3])+"+-"+str(sqrt(err[3]))+"\n")

# print "Signal:",Events[N-1],"Background:",Events[N],"Data:",Events[N+1]
for i in range(N+N_sig+1):
    # print process[i],Events[i]
    # f.write(process[i]+" "+str(Events[i])+"\n")
    f.write("\n"+process[i])
    f.write("\n"+"Total:"+str(Events[0][i])+"+-"+str(sqrt(Err[0][i]))+",Raw number:"+str(Raw[0][i]))
    f.write("\n"+"Isolation tight:"+str(Events[1][i])+"+-"+str(sqrt(Err[1][i])))
    f.write("\n"+"l3 pt cut:"+str(Events[2][i])+"+-"+str(sqrt(Err[2][i])))
    f.write("\n"+"Mll_1 cut:"+str(Events[3][i])+"+-"+str(sqrt(Err[3][i])))
    f.write("\n"+"Final:"+str(Events[3][i])+"+-"+str(sqrt(Err[3][i]))+",Raw number:"+str(Raw[3][i])+"\n")
    if i < N:
        for j in range(len(plot_name)):
            hs[j].Add(h[i][j])

ROOT.gStyle.SetOptStat("")
# c1.SetLogy()
# c1.SetLogx()

for i in range(len(plot_name)):
    # pad = ROOT.TPad(plot_name[i],plot_name[i],0 ,0.3 ,1 ,1)
    # frame = ROOT.gPad.DrawFrame(xrange[i][0],1e-5,xrange[i][1],scale)

    # pad.Draw()
    # ROOT.gPad.Update()
    
    # frame.GetXaxis().SetNdivisions(505, kTRUE)
    # frame.GetYaxis().SetTitle("Events")
    # frame.GetXaxis().SetTitle(h_name[i])


    h_ref = ROOT.TH1F("ref","ref",nbin[i],xrange[i][0],xrange[i][1])
    h_ref.GetXaxis().SetTitle(h_name[i])
    h_ref.GetYaxis().SetTitle("Events")
    # h_ref.GetYaxis().SetRangeUser(1e-5,hs[i].GetMaximum()*3.0/2)
    h_ref.GetYaxis().SetRangeUser(1e-5,h[N][i].GetMaximum()*3.0/2)
    # h_ref.GetYaxis().SetRangeUser(1e-1,hs[i].GetMaximum()*200)
    # h_ref.GetYaxis().SetRangeUser(1e-1,h[N][i].GetMaximum()*200)
    h_ref.Draw("hist")

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

    leg = ROOT.TLegend(0.92,0.92,0.7,0.92-len(process)*0.03)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.03)

    # leg.AddEntry(h[0][i],process[0],"f")

    # hs[i].Draw("same hist")

    for j in range(N+N_sig):
        if j < N:
            # continue
            leg.AddEntry(h[j][i],process[j],"f")
            # h[j][i].Draw("same hist")
        # leg.AddEntry(h[j+1][i],process[j+1],"f")
        else:
            leg.AddEntry(h[j][i],process[j],"f")
            h[j][i].Draw("same hist")

        # h[j+1][i].Draw("same hist")

    h[N+N_sig][i].Draw("same E1")
    # leg.AddEntry(h[N+N_sig][i],"Data","lep")
    leg.Draw("same")

    c1.SaveAs("Plots/" + plot_name[i] + ".png")
    c1.SaveAs("Plots/" + plot_name[i] + ".svg")
    c1.SaveAs("Plots/" + plot_name[i] + ".pdf")
    c1.Clear()

ed = time.time()
print "time consumption:",ed - op,"s"