import ROOT
import os
import sys
from array import array

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

c1 = ROOT.TCanvas("c1", "canvas", 1000, 800)

c1.SetLeftMargin(0.15)
# c1.SetRightMargin(0.05)
c1.SetTopMargin(0.05)
c1.SetBottomMargin(0.15)
c1.SetRightMargin(0.20)


indir = "Data/"
plotPath = "Plots/"

x = sys.argv[1]
f_list = os.listdir(indir)

# x = []
# y = []
i = 0
step = 2.5  #   2.5 um

for f in f_list:
    A = [[],[],[]]
    # if not "csv" in f and not "txt" in f:    continue
    if not "txt" in f:    continue
    if "charge" in f:    continue

    out = f.replace(".csv","").replace(".txt","")

    # f_out = open(plotPath + out + "_" + x + ".csv","a")
    f_out = open(plotPath + out + ".csv","a")
    f_out.seek(0)
    f_out.truncate()

    f = os.path.join(indir,f)
    print f
    with open(f,"r") as text:
        content = text.read().splitlines()
    
    x0 = -1
    y0 = -1
    # i = 0
    for p in content:
        # i+=1
        # if i > 10: break
        p = p.split(",")
        # if int(p[0]) > x0:    x0 = int(p[0])
        # if int(p[1]) > y0:    y0 = int(p[1])
        # A.append( float(p[8]) )
        # if int(p[0][3:6]) > x0:    x0 = int(p[0][3:6])
        # if int(p[0][7:]) > y0:    y0 = int(p[0][7:])
        A[0].append( int(p[0][3:6]) )
        A[1].append( int(p[0][7:]) )
        A[2].append( float(p[1]) )

        if int(p[0][3:6]) == int(x):
        # if int(p[0][7:]) == int(x):
        # if int(p[0][7:]) == 2000:
            # if int(p[0][7:]) < 2079 and int(p[0][7:]) > 2040:
            f_out.write( p[0][3:6] + "," + p[0][7:] + "," + p[1] + "\n")

    
    # x0 = x0 - int(content[0][3:6])
    # y0 = y0 - int(content[0][7:11])

    print(min(A[0]),min(A[1]),len(A[2]))

    h_2d = ROOT.TH2D("2d","2d", max(A[0])-min(A[0]), 0, (max(A[0])-min(A[0]))*step, max(A[1])-min(A[1]), 0, (max(A[1])-min(A[1]))*step )

    for i in range(len(A[0])):
        # if i > 10: break
        # print A[0][i]-min(A[0]), A[1][i]-min(A[1]), A[2][i]
        h_2d.SetBinContent( A[0][i]-min(A[0]), A[1][i]-min(A[1]), A[2][i] )

    # for i in range(x0):
    #     for j in range(y0):
    #         h_2d.SetBinContent(i,j,A[i+j/y0*x0])

    h_2d.GetXaxis().SetTitle("X [um]")
    h_2d.GetYaxis().SetTitle("Y [um]")
    h_2d.GetZaxis().SetTitle("Ampl [a.u.]")
    # h_2d.GetXaxis().SetNdivisions(505, ROOT.kTRUE)
    # h_2d.SetMarkerSize(0.2)
    h_2d.Draw("COLZ")

    # c1.SaveAs(plotPath + out + "_2d.png")
    # c1.SaveAs(plotPath + out + "_2d.svg")
    c1.Clear()
    
####    fraction calculation
f_out = open("fra_"+ str(x) +".csv","a")
f_out.seek(0)
f_out.truncate()

Pos = [[],[]]
Var = [[],[],[],[]]

for i in range(4):
    f = plotPath + "amplitude_" + str(i+1) + ".csv"
    print f
    with open(f,"r") as text:
        content = text.read().splitlines()

        for var in content:
            var = var.split(",")
            Var[i].append(float(var[2]))
            if i == 0:
                Pos[0].append(var[0])
                Pos[1].append(var[1])

for i,var in enumerate(Var[0]):
    # print Pos[0][i],Pos[1][i],Var[0][i]
    fra1 = (Var[0][i]-Var[1][i])/(Var[0][i]+Var[1][i])
    # fra2 = (Var[2][i]-Var[3][i])/(Var[2][i]+Var[3][i])
    fra2 = (Var[0][i]+Var[2][i]-Var[1][i]-Var[3][i])/(Var[0][i]+Var[2][i]+Var[1][i]+Var[3][i])
    # fra2 = (Var[0][i]-Var[2][i]+Var[1][i]-Var[3][i])/(Var[0][i]+Var[2][i]+Var[1][i]+Var[3][i])
    f_out.write( Pos[0][i] + "," + Pos[1][i] + "," + str(fra1) + "," + str(fra2) + "\n" )

f_out.close()