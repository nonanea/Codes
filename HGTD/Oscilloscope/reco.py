import ROOT
import os

# f_list = ["fra_830.csv","fra_835.csv","fra_840.csv","fra_845.csv"]
# dir = "./temp/"
dir = "./"
f_list = os.listdir(dir)

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)

c1 = ROOT.TCanvas("c1", "canvas", 1000, 800)

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.05)
c1.SetTopMargin(0.05)
c1.SetBottomMargin(0.15)

h = ROOT.TH1F("h","h",20,-25,25)

reco_f = ROOT.TF1("f","[0]*ROOT::Math::erfc((x-[1])*[2])+[3]",0,400)
reco_f.SetParameters(-0.778208,205.205,-0.0148531,0.762492)

# reco_f = ROOT.TF1("f","[0]*(TMath::Gaus(x,[1]+[4],[2])-TMath::Gaus(x,-[1]+[4],[2]))/(TMath::Gaus(x,[4],[2])+TMath::Gaus(x,[1]+[4],[2])+TMath::Gaus(x,-[1]+[4],[2]))+[3]",0,400)
# reco_f.SetParameters(0.775334,-63.9891,34.8793,-0.0150533,205.111)

for f in f_list:
    if "fra" not in f or "csv" not in f: continue
    f = dir + f
    print f
    with open(f,"r") as text:
        content = text.read().splitlines()
        
        for i,y in enumerate(content):
            # if i < 60 or i > 100: continue
            y = y.split(",")
            if int(y[0]) < 820 or int(y[0]) > 860: continue

            # print i*2.5,reco_f.GetX(float(y[3]),-10,410)
            # h.Fill(i*2.5-reco_f.GetX(float(y[3]),-10,410))
            
            # print (int(y[0])-760)*2.5,float(y[3])*75
            # h.Fill(i*2.5-203+float(y[3])*75)
            h.Fill((int(y[0])-760)*2.5-197+float(y[3])*75)
            # print i*2.5,float(y[3])*75.0


h.GetXaxis().SetTitle("x_{reco}-x_{truth} [um]")
h.GetYaxis().SetTitle("Events")
h.Draw("hist")

gaus = ROOT.TF1("g","gaus(0)",-30,30)
gaus.SetLineColor(ROOT.kRed)

r = h.Fit("gaus","QS","C",-10,10)
gaus.SetParameters(r.Value(0),r.Value(1),r.Value(2))
gaus.Draw("same")

c1.SaveAs("reco.png")
c1.SaveAs("reco.svg")
c1.Clear()