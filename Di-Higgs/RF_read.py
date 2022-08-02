import ROOT
import os
import sys

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptTitle(0)

colorList = [
    ROOT.kRed,
    ROOT.kBlue,
    ROOT.kGreen,
    ROOT.kMagenta,
    ROOT.kViolet,
    ROOT.kOrange,
    ROOT.kCyan,
    ROOT.kYellow,
    ROOT.kTeal,
    ROOT.kAzure,
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

c1 = ROOT.TCanvas("c1", "hist", 1000, 800)

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.05)
c1.SetTopMargin(0.05)
c1.SetBottomMargin(0.15)

leg = ROOT.TLegend(0.92,0.92,0.7,0.68)
# frame = ROOT.TH1()

leg.SetBorderSize(0)
leg.SetTextSize(0.05)

dataPath = "Data/"
plotPath = "Plots/"
# obj_name = "MUON_ID__1up"
obj_name = "CategoryReduction_JET_JER_EffectiveNP_9__1up"
obj_list = []

fileList = os.listdir(dataPath)
# print(fileList)
flag = 0

h = [ROOT.TH1F() for j in range(7)]
h_ref = ROOT.TH1F("ref", "ref", 20, -1., 1.)
h_ref.SetLineColor(ROOT.kBlack)

i = -1
ROOT.gStyle.SetOptStat("")

leg = ROOT.TLegend(0.92,0.92,0.7,0.92-7*0.03)
leg.SetBorderSize(0)
leg.SetTextSize(0.03)

for f_name in fileList:
    if "ttV.root" not in f_name:    continue
    # if "ZH.root" not in f_name and "ggF_HH.root" not in f_name:    continue
    i += 1
    
    process = f_name.replace(".root","")
    h[i] = ROOT.TH1F(process, process, 20, -1., 1.)
    h[i].SetLineColor(colorList[i])
    # h[i].SetLineStyle(ROOT.kDashed)
    h[0].GetXaxis().SetTitle("BDT")
    h[0].GetYaxis().SetTitle("Events")
    h[0].GetYaxis().SetRangeUser(0,25)
    print h[i]

    f = ROOT.TFile.Open(dataPath + f_name,"read")
    t = f.Get(obj_name)
    
    for event in t:
        if "ggF_HH.root" in f_name:
            h[i].Fill(event.SRBDT,event.weight*100)
        else:
            h[i].Fill(event.SRBDT,event.weight)

    for j in range(h[i].GetSize()):
        # h_ref.SetBinContent(j,h[i].GetBinContent(j)*1.1)
        # h_ref.SetBinError(j,h[i].GetBinError(j))
        # h_ref.Fill(-1.05+0.2*j,h[1].GetBinContent(j)*1.1)
        print j,h[i].GetBinContent(j),h[i].GetBinError(j)

    if i == 0:
        h[i].Draw("hist")
    else:
        h[i].Draw("same hist")

    # h_ref.Draw("same hist")

    leg.AddEntry(h[i],process,"l")

    t.Delete("")
    f.Close()

leg.Draw("same")
c1.SaveAs(plotPath+"_test.png")
c1.SaveAs(plotPath+"_test.svg")
# for i in range(len(
#     if "root" in fileList[i]:
#         # print(fileList[i],i)
#         f = ROOT.TFile.Open(dataPath + fileList[i],"read")

#         # print(f)
        
#         obj = f.Get(obj_name)

#         # h = ROOT.TH1F("BDT","BDT",15,-0.6,0.4)

        
#         if i == 0:
#             h.Draw("hist")
#         else:
#             h.Draw("same hist")
        
#         # c1.Update()
#         # h.Delete("")

        # objj.SetDirectory(0)
        # objj = f.Get(obj_name).Clone()
        # f.Close()
        # print(objj)
        # f.GetObject(obj_name,obj)

        # obj_list.append(obj)
        # print(obj_list)
        # obj.SetLineColor(i-2)