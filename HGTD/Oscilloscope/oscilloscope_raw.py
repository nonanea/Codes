import lecroyparser
import argparse
import ROOT
import os
import sys
import numpy as np
import csv
from array import array
import time

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# colorList = [ kRed, kBlue, kGreen, kViolet, kMagenta, kTeal, kPink, kAzure, kCyan, kYellow, kBlue + 1, kRed + 1, kGreen + 1, kTeal + 1, kPink + 1, kMagenta + 1, kViolet + 1, kAzure + 1, kCyan + 1, kYellow + 1]
plotStorePath = "Plots/"

c1 = ROOT.TCanvas("c1", "waveform", 1000, 800)
c1.SetLeftMargin(0.25)
c1.SetRightMargin(0.2)
c1.SetTopMargin(0.2)

def GetOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-dir", "--DataDir", help="Dir of data")
    parser.add_argument("-averageOnly", "--AverageOnly", dest = 'AverageOnly', action = 'store_true', help="Only draw the average of waveform")
    parser.add_argument("-singleDraw", "--SingleDraw", dest = 'SingleDraw', action = 'store_true', help="Draw waveforms one by one")
    parser.add_argument("-drawNum", "--DrawNum", type = int, default = 50, help="Total number of the drawn waveforms")
    parser.add_argument("-multi", "--Multiple", dest = 'Multiple', action = 'store_true', help="Draw multiple waveforms together")
    parser.add_argument("-multiNum", "--MultiNum", type = int, default = 50, help="Total number of the waveforms per plot")
    options = parser.parse_args(args)
    return options

def PedestalAndNoise(y):
    # global ped,noise,noise_amp

    noise_list = []
    noise_amp = 0

    length = len(y)
    for i in range(int(length*0.2)):
        noise_list.append(y[i])
        if y[i] > noise_amp:
            noise_amp = y[i]
    
    ped = np.mean(noise_list)
    noise = np.std(noise_list,ddof = 1)
    noise_amp = noise_amp - ped

    return(ped, noise, noise_amp)

def Peak(x,A,threshold,init_f,fin_f):
    # global Tmax,Amax,Max

    Amax = 0
    Tmax = 0
    Max = 0

    length = len(x)
    for i in range(int(length*init_f),int(length*fin_f)):
        if A[i] > Amax:
            Amax = A[i]
            Tmax = x[i]
            Max = i
    
    # print(Amax,Tmax,Max)
    return(Amax, Tmax, Max)

def LeftTime(x,A,threshold,Max):
    if x[0] == x[1]:
        return(0)

    i = Max
    while A[i] > threshold and i > 0:
        i-=1
    
    a=(A[i+1]-A[i])/(x[i+1]-x[i])
    if not a == 0:
        b=A[i]-a*x[i]
        lefttime=(threshold-b)/a
    else:
        lefttime = x[i]

    return(lefttime)

def RightTime(x,A,threshold,Max):
    if x[0] == x[1]:
        return(0)
    
    i = Max
    while A[i] > threshold and i < len(A) - 2:
        i+=1
    
    a=(A[i]-A[i-1])/(x[i]-x[i-1])
    b=A[i]-a*x[i]
    righttime=(threshold-b)/a

    return(righttime)

def Charge(x,A,Max,left,right,threshold):
    width = int(right - left)
    deltaT = x[1] - x [0]
    charge = 0

    for i in range(int(Max-0.5*width),int(Max+width)):
        charge+=A[i]

    return(charge*deltaT)

#########Files processing#########
op = time.time()

options = GetOptions(sys.argv[1:])

filePath = options.DataDir

plotType = filePath

if not os.path.exists(plotStorePath + plotType):
    os.mkdir(plotStorePath + plotType)

results_file = plotStorePath + plotType + "Average.txt"
f = open(results_file,"a")
f.seek(0)
f.truncate()

plotname = filePath
plotname = plotname.replace("/","")

filePath = "Data/" + filePath
print("\nRead files from",filePath)

fileList = os.listdir(filePath)
fileList = sorted(fileList)

ped = array('f',4*[0.])
noise = array('f',4*[0.])
noise_amp = array('f',4*[0.])
charge = array('f',4*[0.])
charge_hw = array('f',4*[0.])
charge_dw = array('f',4*[0.])
Amax = array('f',4*[0.])
Tmax = array('f',4*[0.])
Max = array('I',4*[0])
CTD = array('f',4*[0.])
CFD = array('f',4*[0.])
TOA = array('f',4*[0.])
TOT = array('f',4*[0.])
risetime = array('f',4*[0.])
risetime_ave = array('f',4*[0.])
ave_risetime = array('f',4*[0.])
# slope = array('f',4*[0.])
id = array('f',4*[0.])

# ped = 4*[0]
# noise = 4*[0]
# noise_amp = 4*[0]
# charge = 4*[0]
# charge_hw = 4*[0]
# charge_dw = 4*[0]
# Amax = 4*[0]
# Tmax = 4*[0]
# Max = 4*[0]
# CTD = 4*[0]
# CFD = 4*[0]
# TOA = 4*[0]
# TOT = 4*[0]
# risetime = 4*[0]
# risetime_ave = 4*[0]
# ave_risetime = 4*[0]
# slope = 4*[0]
# id = 4*[0]

t = np.zeros([4,9999])
sum = np.zeros([4,9999])

i = 0
length = 9999
ch_list = []
ch_dut = 1
ch_ref = 2
Ns = 500
threshold = 0.0
init_f = 0.4
fin_f = 0.6
fra = 0.5
th = 10

# options = sys.argv[2:]
# print(options)
# average_only = False
# average_only = True
# draw = False
# draw = True
# multi = True

list_ch = [[] for i in range(4)]

t_ref = []

for file in fileList:
    if ".trc" in file:
        if "C1" in file:
            list_ch[0].append(filePath+file)
        elif "C2" in file:
            list_ch[1].append(filePath+file)
        elif "C3" in file:
            list_ch[2].append(filePath+file)
        else:
            list_ch[3].append(filePath+file)

for i in range(3):
    if list_ch[i]:
        ch_list.append(i)

# Ref processing
for file in list_ch[2]:
    if len(list_ch[2]) < 1:
        print("No reference!")
        break

    x_ave = []
    y_ave = []
    waveform = lecroyparser.ScopeData(file)
    for n in range(len(waveform.x)):
        x_ave.append(waveform.x[n]*1e9)
        y_ave.append(-waveform.y[n])
    
    if length > len(x_ave):
        length = len(x_ave)

    pulse_length = length/5 # half pulse length
    
    (ped[ch_ref],noise[ch_ref],noise_amp[ch_ref]) = PedestalAndNoise(y_ave)
    (Amax[ch_ref],Tmax[ch_ref],Max[ch_ref]) = Peak(x_ave,y_ave,threshold,init_f,fin_f)

    if Amax[ch_ref]*1000 < 10:
        continue
    
    n = Max[ch_ref]
    while y_ave[n] > (Amax[ch_ref]+ped[ch_ref])*0.5:  
        n-=1
    Max[ch_ref] = n

    t_ref.append(Max[ch_ref])

    # print(int(Max[ch_ref]-length/5), int(Max[ch_ref]+length/5))
    for n in range(int(Max[ch_ref]-pulse_length), int(Max[ch_ref]+pulse_length)):
        t[ch_ref][int(n-Max[ch_ref]+pulse_length)]=x_ave[n]-x_ave[Max[ch_ref]]
        sum[ch_ref][int(n-Max[ch_ref]+pulse_length)]=sum[ch_ref][int(n-Max[ch_ref]+pulse_length)]+y_ave[n]-ped[ch_ref]

    i+=1
    if i > Ns:
        break

print(Ns,"signals of the reference")

i = 0
# DUT processing
for j in range(len(t_ref)):
    if len(list_ch[1]) < 1:
        print("No DUT!")
        break

    x_ave = []
    y_ave = []
    waveform = lecroyparser.ScopeData(list_ch[1][j])
    for n in range(len(waveform.x)):
        x_ave.append(waveform.x[n]*1e9)
        y_ave.append(-waveform.y[n])
    
    if length > len(x_ave):
        length = len(x_ave)
    
    (ped[ch_dut],noise[ch_dut],noise_amp[ch_dut]) = PedestalAndNoise(y_ave)
    (Amax[ch_dut],Tmax[ch_dut],Max[ch_dut]) = Peak(x_ave,y_ave,threshold,init_f,fin_f)

    if Amax[ch_dut]*1000 < 10:
        continue
    
    for n in range(int(t_ref[j]-pulse_length), int(t_ref[j]+pulse_length)):
        t[ch_dut][int(n-t_ref[j]+pulse_length)]=x_ave[n]-x_ave[t_ref[j]]
        sum[ch_dut][int(n-t_ref[j]+pulse_length)]=sum[ch_dut][int(n-t_ref[j]+pulse_length)]+y_ave[n]-ped[ch_dut]

    i+=1
    
print(i,"signals of the DUT")

gr_ave = [ROOT.TGraph() for i in range(4)]
T_ave = [[] for i in range(4)]
A_ave = [[] for i in range(4)]

for n in range(4):
    for i in range(int(2*pulse_length)):
        T_ave[n].append(t[n][i])
        A_ave[n].append(sum[n][i]/Ns)
        gr_ave[n].SetPoint(i,t[n][i],sum[n][i]/Ns)

P_left = 4*[0]
P_right = 4*[0]

for i in range(4):
    (Amax[i],Tmax[i],Max[i]) = Peak(T_ave[i],A_ave[i],threshold,0.2,0.8)
    
    n = Max[i]
    while A_ave[i][n] > Amax[i]*0.1:
        n-=1
    P_left[i] = n
    
    n = Max[i]
    while A_ave[i][n] > Amax[i]*0.1:
        # print(n)
        n+=1
        if n == len(A_ave[i]):
            break
    P_right[i] = n

    print("Length:",len(A_ave[i]),"M:",Max[i],"L:",P_left[i],"R:",P_right[i])
    f.write("Length:"+str(len(A_ave[i]))+",M:"+str(Max[i])+",L:"+str(P_left[i])+",R:"+str(P_right[i])+"\n")
    
    y_max = ROOT.TMath.MaxElement(gr_ave[i].GetN(),gr_ave[i].GetY())

    gr_ave[i].SetMarkerStyle(21)
    gr_ave[i].SetMarkerSize(0.4)
    gr_ave[i].SetMarkerColor(2)
    gr_ave[i].SetLineColor(2)

    gr_ave[i].GetXaxis().SetTitle("Time [ns]")
    gr_ave[i].GetYaxis().SetTitle("Ampl [V]")

    gr_ave[i].GetXaxis().SetNdivisions(505, ROOT.kTRUE)
    gr_ave[i].GetXaxis().SetLimits(-7, 7)
    # gr_ave[i].GetYaxis().SetRangeUser(0.0025, 0.0045)
    gr_ave[i].GetYaxis().SetRangeUser(-0.08*y_max, 1.1*y_max)

    gr_ave[i].Draw("ALP")
    
    rise_start = LeftTime(T_ave[i],A_ave[i],Amax[i]*0.4,Max[i])
    rise_end = LeftTime(T_ave[i],A_ave[i],Amax[i]*0.6,Max[i])
    risetime_ave[i] = rise_end - rise_start
    
    print("Amax_20:",format((Amax[i]*0.2*1000),'.6f'),"Rise time:",format((rise_end-rise_start),'.6f'))

    pt = ROOT.TLatex()
    pt.SetNDC()
    pt.SetTextAlign(12)
    pt.SetTextSize(0.02)
    
    # text = "Amax_20:" + str(Amax[i]*0.2*1000)
    pt.DrawText(0.6,0.74,"Amax_20:" + format((Amax[i]*0.2*1000),'.6f'))
    # text = "Rise time_4060:" + str(rise_end-rise_start)
    pt.DrawText(0.6,0.7,"Rise time_4060:" + format((rise_end-rise_start),'.6f'))

    f.write("Amax_20:"+format((Amax[i]*0.2*1000),'.6f')+",Rise time_4060:"+format((rise_end-rise_start),'.6f')+"\n")
    
    c1.Update()
    c1.SaveAs(plotStorePath + plotType + str(i) + "average.png")
    c1.SaveAs(plotStorePath + plotType + str(i) + "average.svg")
    c1.Clear()

if options.AverageOnly:
    print("Only Average Drawn")
    sys.exit(0)

gr = ROOT.TGraph()
f = ROOT.TFile(plotStorePath + plotType + plotname + ".root","recreate")
# tt = [ROOT.TTree() for i in range(4)]
# tt[0] = ROOT.TTree("t1","channel 1")
# tt[1] = ROOT.TTree("t2","channel 2")
# tt[2] = ROOT.TTree("t3","channel 3")
# tt[3] = ROOT.TTree("t4","channel 4")

tt = ROOT.TTree("tt","All 4 Channels")

nCh = array( 'I', [4])

# for n in range(4):
    # tt[n].Branch("Amax",Amax[n],"Amax[{}]/F".format(n))
    # tt[n].Branch("Tmax",Tmax[n],"Tmax[{}]/F".format(n))
    # tt[n].Branch("Charge",charge[n],"charge[{}]/F".format(n))
    # tt[n].Branch("Charge_HW",charge_hw[n],"charge_hw[{}]/F".format(n))
    # tt[n].Branch("Charge_DW",charge_dw[n],"charge_dw[{}]/F".format(n))
    # tt[n].Branch("CTD",CTD[n],"CTD[{}]/F".format(n))
    # tt[n].Branch("CFD",CFD[n],"CFD[{}]/F".format(n))
    # # tt[n].Branch("ZCD",ZCD[n]) 
    # tt[n].Branch("TOT",TOT[n],"TOT[{}]/F".format(n))
    # tt[n].Branch("risetime_ave",ave_risetime[n],"ave_risetime[{}]/F".format(n))
    # tt[n].Branch("risetime",risetime[n],"risetime[{}]/F".format(n))
    # tt[n].Branch("slope",slope[n],"slope[{}]/F".format(n))
    # tt[n].Branch("pedestal",ped[n],"ped[{}]/F".format(n))
    # # tt[n].Branch("ped_1",ped_1[n])
    # # tt[n].Branch("ped_2",ped_2[n])
    # tt[n].Branch("Noise",noise[n],"noise[{}]/F".format(n))
    # tt[n].Branch("Noise_Amp",noise_amp[n],"noise_amp[{}]/F".format(n))
    # tt[n].Branch("ID",id[n],"id[{}]/F".format(n))
    
    # tt[n].Branch('nCh',nCh,'nCh/I')
    # tt[n].Branch("Amax",Amax[n],"Amax[nCh]/F")
    # tt[n].Branch("Tmax",Tmax[n],"Tmax[nCh]/F")
    # tt[n].Branch("Charge",charge[n],"charge[nCh]/F")
    # tt[n].Branch("Charge_HW",charge_hw[n],"charge_hw[nCh]/F")
    # tt[n].Branch("Charge_DW",charge_dw[n],"charge_dw[nCh]/F")
    # tt[n].Branch("CTD",CTD[n],"CTD[nCh]/F")
    # tt[n].Branch("CFD",CFD[n],"CFD[nCh]/F")
    # # tt[n].Branch("ZCD",ZCD[n]) 
    # tt[n].Branch("TOT",TOT[n],"TOT[nCh]/F")
    # tt[n].Branch("risetime_ave",ave_risetime[n],"ave_risetime[nCh]/F")
    # tt[n].Branch("risetime",risetime[n],"risetime[nCh]/F")
    # tt[n].Branch("slope",slope[n],"slope[nCh]/F")
    # tt[n].Branch("pedestal",ped[n],"ped[nCh]/F")
    # # tt[n].Branch("ped_1",ped_1[n])
    # # tt[n].Branch("ped_2",ped_2[n])
    # tt[n].Branch("Noise",noise[n],"noise[nCh]/F")
    # tt[n].Branch("Noise_Amp",noise_amp[n],"noise_amp[nCh]/F")
    # tt[n].Branch("ID",id[n],"id[nCh]/F")
    
tt.Branch('nCh',nCh,'nCh/I')
tt.Branch("Amax",Amax,"Amax[nCh]/F")
tt.Branch("Tmax",Tmax,"Tmax[nCh]/F")
tt.Branch("Charge",charge,"charge[nCh]/F")
tt.Branch("Charge_HW",charge_hw,"charge_hw[nCh]/F")
tt.Branch("Charge_DW",charge_dw,"charge_dw[nCh]/F")
tt.Branch("CTD",CTD,"CTD[nCh]/F")
tt.Branch("CFD",CFD,"CFD[nCh]/F")
# tt[n].Branch("ZCD",ZCD[n]) 
tt.Branch("TOT",TOT,"TOT[nCh]/F")
tt.Branch("risetime_ave",ave_risetime,"ave_risetime[nCh]/F")
tt.Branch("risetime",risetime,"risetime[nCh]/F")
# tt[n].Branch("slope",slope,"slope[nCh]/F")
tt.Branch("pedestal",ped,"ped[nCh]/F")
# tt[n].Branch("ped_1",ped_1[n])
# tt[n].Branch("ped_2",ped_2[n])
tt.Branch("Noise",noise,"noise[nCh]/F")
tt.Branch("Noise_Amp",noise_amp,"noise_amp[nCh]/F")
tt.Branch("ID",id,"id[nCh]/F")

gr = ROOT.TGraph()
hA_T = [ROOT.TH2D() for i in range(4)]
t_range = [-10,10]

for n in range(4):
    hA_T[n] = ROOT.TH2D("Amax vs. Tmax","Amax vs. Tmax",100,0,200,40,t_range[0],t_range[1])

for i in range(len(list_ch[1])):
    # if i > 50:
    #     break
    # if not "trc" in file:
    #     continue
    for j in ch_list:
        x = []
        y = []
        T = []
        A = []

        plotInfo = list_ch[j][i]
        plotInfo = plotInfo[plotInfo.find("/C")+1:plotInfo.rfind(".trc")]

        # print(i)
        waveform = lecroyparser.ScopeData(list_ch[j][i])
        for n in range(len(waveform.x)):
            x.append(waveform.x[n]*1e9)
            y.append(-waveform.y[n])
        
        if length > len(x):
            length = len(x)

        (ped[j],noise[j],noise_amp[j]) = PedestalAndNoise(y)
        for n in range(length):
            A.append(y[n]-ped[j])
            if options.SingleDraw or options.Multiple:
                gr.SetPoint(n,x[n],A[n])

        inif_f = 0.1
        fin_f = 0.9
        if j == ch_dut or j == ch_ref:
            init_f = 0.4
            fin_f = 0.6

        (Amax[j],Tmax[j],Max[j]) = Peak(x,A,threshold,init_f,fin_f)

        charge[j] = Charge(x,A,Max[j],P_left[j],P_right[j],threshold)
        charge_hw[j] = Charge(x,A,Max[j],0.8*P_left[j],0.8*P_right[j],threshold)
        # charge_dw[j] = Charge(x,A,Max[j],3*P_left[j],3*P_right[j],threshold)
        charge_dw[j] = Charge(x,A,Max[j],2*P_left[j],2*P_right[j],threshold)
        
        CTD[j] = LeftTime(x,A,th,Max[j])
        # CTD_C = rightTime(xlist,Alist,th,Max)
        # TOT = CTD_C - CTD; 
        CFD[j] = LeftTime(x,A,Amax[j]*fra,Max[j])

        rise_start = LeftTime(x,A,Amax[j]*0.4,Max[j])
        rise_end = LeftTime(x,A,Amax[j]*0.6,Max[j])
        risetime[j] = rise_end - rise_start
        # slope[j] = Amax[j]*0.2/risetime[j]
        TOT[j] = RightTime(x,A,Amax[j]*0.1,Max[j]) - LeftTime(x,A,Amax[j]*0.1,Max[j])

        ave_risetime[j] = risetime_ave[j]
        
        id[j] = float(file[-9:-4])

        hA_T[j].Fill(Amax[j]*1000,Tmax[j])

        if (options.SingleDraw or options.Multiple) and i < options.DrawNum:
            gr.SetMarkerStyle(21)
            gr.SetMarkerSize(0.4)
            gr.SetMarkerColor(2)
            gr.SetLineColor(2)

            gr.GetXaxis().SetTitle("Time [ns]")
            gr.GetYaxis().SetTitle("Ampl [V]")

            gr.GetXaxis().SetNdivisions(505, ROOT.kTRUE)
            gr.GetXaxis().SetLimits(-20, 20)
            # gr.GetYaxis().SetRangeUser(0.0025, 0.0045)
            gr.GetYaxis().SetRangeUser(-0.8*Amax[j], 1.1*Amax[j])

            # gr.Draw("ALP")
            
            if not options.Multiple:
                gr.Draw("ALP")
                c1.SaveAs(plotStorePath + plotType + plotInfo + ".png")
                c1.Clear()
            
            else:
                # print(i)
                if i != 0:
                    gr.Draw("LP")
                else:
                    gr.Draw("ALP")

                if (i+1)%options.MultiNum == 0:
                    # print(i+1)
                    c1.SaveAs(plotStorePath + plotType + str(j) + ".png")
                    c1.Clear()

    tt.Fill()

    # i+=1
    print("Processing:{}%".format(round((i + 1) * 100 / len(list_ch[1]),2)), end="\r")

print("File number per Ch:",i)

tt.Write()
for n in range(4):

    hA_T[n].GetXaxis().SetTitle("Amax [mV]")
    hA_T[n].GetYaxis().SetTitle("Tmax [ns]")
    hA_T[n].GetXaxis().SetNdivisions(505, ROOT.kTRUE)
    hA_T[n].SetMarkerSize(0.2)
    hA_T[n].Draw("COLZ")

    c1.SaveAs(plotStorePath + plotType + str(n) + "Amax_vs_Tmax.png")
    c1.Clear()

ed = time.time()
print("time consumption:",ed - op,"s")