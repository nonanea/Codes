import ROOT
import os
from array import array
from math import fabs,sqrt,pow

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptTitle(0)

c1 = ROOT.TCanvas("c1", "canvas", 1000, 800)

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.05)
c1.SetTopMargin(0.05)
c1.SetBottomMargin(0.15)

# f_l = "mmm/data_eta_L.txt"
# f_t = "mmm/data_eta_T.txt"

f_l = "loose_eta.txt"
f_t = "tight_eta.txt"

# bins =  array("f",[0,8,10,12,14,16,18,20,80])
bins =  array("f",[-2.5,-1.5,-0.5,0.5,1.5,2.5])
# bins =  array("f",[0.,0.5,1.5,2.5,3])
# x = array("f",[5.,15.,25.,35.,45.,60.])
# ex = array("f",[5.,5.,5.,5.,5.,5.])

x = [array("f",(len(bins)-1)*[0.]) for i in range(2)]
ex = [array("f",(len(bins)-1)*[0.]) for i in range(2)]
y = [array("f",(len(bins)-1)*[0.]) for i in range(2)]
ey = [array("f",(len(bins)-1)*[0.]) for i in range(2)]

l_y = [array("f",(len(bins)-1)*[0.]) for i in range(2)]
l_ey = [array("f",(len(bins)-1)*[0.]) for i in range(2)]
t_y = [array("f",(len(bins)-1)*[0.]) for i in range(2)]
t_ey = [array("f",(len(bins)-1)*[0.]) for i in range(2)]

# h = ROOT.TH1F("h","h",7,bins)
# data = 1
# if data:    j = 6
# else:   j = 2

for k in range(2):
    j = 4*k + 2
    with open(f_l,"r") as f:
        content = f.read().splitlines()

        for i in range(len(l_y[0])):
            # print var[0],var[1],var[2]
            # for i,bin in enumerate(bins):
            for var in content:
                var = var.split(",")
                # print(var)
                if float(var[1]) < bins[i+1] and float(var[1]) >= bins[i] :
                    l_y[k][i] += float(var[j])
                    l_ey[k][i] += float(var[j+1])*float(var[j+1])
                    # break
            
            x[k][i] = (bins[i+1]+bins[i])/2
            ex[k][i] = (bins[i+1]-bins[i])/2
            l_ey[k][i] = sqrt(l_ey[k][i])

    with open(f_t,"r") as f:
        content = f.read().splitlines()

        for i in range(len(t_y[0])):
            for var in content:
                var = var.split(",")
                if float(var[1]) < bins[i+1] and float(var[1]) >= bins[i] :
                    t_y[k][i] += float(var[j])
                    t_ey[k][i] += float(var[j+1])*float(var[j+1])
                    # break
            
            t_ey[k][i] = sqrt(t_ey[k][i])
            
            # print t_y[i],l_y[i]

            y[k][i] = t_y[k][i]/l_y[k][i]
            ey[k][i] = sqrt( pow(t_ey[k][i]/l_y[k][i],2) + pow(t_y[k][i]/l_y[k][i]/l_y[k][i]*l_ey[k][i],2) )

            print "[",bins[i],",",y[k][i],",",ey[k][i],"],"

# print bin_content
# print x,y,ey

gr1 = ROOT.TGraphAsymmErrors(len(x[0]),x[0],y[0],ex[0],ex[0],ey[0],ey[0])
gr2 = ROOT.TGraphAsymmErrors(len(x[0]),x[1],y[1],ex[1],ex[1],ey[1],ey[1])

gr1.SetLineColor(ROOT.kRed)
gr2.SetLineColor(ROOT.kBlue)

# gr1.GetXaxis().SetLimits(0,60)
gr1.GetXaxis().SetLimits(bins[0],bins[-1])
gr1.GetYaxis().SetRangeUser(min(y[0])/2,max(y[0])*1.1)
gr1.GetXaxis().SetTitle("|#eta|")
# gr1.GetXaxis().SetTitle("P_{T} [GeV]")
gr1.GetYaxis().SetTitle("Fake Factor")

leg = ROOT.TLegend(0.88,0.88,0.75,0.80)
leg.SetBorderSize(0)
leg.SetTextSize(0.03)

leg.AddEntry(gr1,"MC Z+jets","l")
leg.AddEntry(gr2,"Data","l")

gr1.Draw("AP")
gr2.Draw("P")
leg.Draw("same")

c1.SaveAs("graph.png")
c1.SaveAs("graph.svg")