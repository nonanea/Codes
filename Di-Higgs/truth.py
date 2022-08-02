import ROOT
from operator import itemgetter

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat("")

c1 = ROOT.TCanvas("c1", "canvas", 1000, 800)

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.05)
c1.SetTopMargin(0.05)
c1.SetBottomMargin(0.15)
f = ROOT.TFile("sample.root")
t = f.Get("tree_NOMINAL")

lep = [ROOT.TLorentzVector(1,1,0,5) for i in range(3)]

h_0 = ROOT.TH1F("0","0",20,110,140)
h_1 = ROOT.TH1F("1","1",20,110,140)
h_2 = ROOT.TH1F("2","2",20,110,140)

h_0.SetLineColor(ROOT.kBlue)
h_1.SetLineColor(ROOT.kRed)
h_2.SetLineColor(ROOT.kGreen)

for event in t:
    # vl = t.v_tlv_L
    # vp_PID = t.vp_PID
    if len(t.vp_PID) != 4: continue
    # if len(t.vp_PID) != 4 or t.nJet < 2: continue
    flag = 0
    l_order = sorted( enumerate( [t.v_tlv_L[0].Pt(),t.v_tlv_L[1].Pt(),t.v_tlv_L[2].Pt(),t.v_tlv_L[2].Pt()] ) ,key=itemgetter(1))
    # print l_order
    # for j, vp in enumerate(t.vp_PID):
    for vp in t.vp_PID:
        # print vp,t.event
        if "24" in str(vp):
            flag = 0
            break
        elif "25" not in str(vp):
            flag = 1
            break
        else: flag = 2

    quad_mass = (t.v_tlv_L[0] + t.v_tlv_L[1] + t.v_tlv_L[2] + t.v_tlv_L[3] ).M()/1000
    if quad_mass > 135 or quad_mass < 115: continue
    # print((t.v_tlv_L[0] + t.v_tlv_L[1] + t.v_tlv_L[2] + t.v_tlv_L[3] ).M())

    if flag == 0: h_0.Fill(quad_mass)
    elif flag == 1: h_1.Fill(quad_mass)
    else: h_2.Fill(quad_mass)

    # if flag == 0: h_0.Fill( t.v_tlv_L[l_order[2][0]].Pt()/1000 )
    # elif flag == 1: h_1.Fill( t.v_tlv_L[l_order[2][0]].Pt()/1000 )
    # else: h_2.Fill( t.v_tlv_L[l_order[2][0]].Pt()/1000 )


h_0.GetXaxis().SetTitle("M_{4l} [GeV]")
# h_0.GetXaxis().SetTitle("p_{T,2} [GeV]")
h_0.GetXaxis().SetNdivisions(505, ROOT.kTRUE)
h_0.GetYaxis().SetTitle("Events")
h_0.GetYaxis().SetRangeUser(1e-5,20)
h_0.Draw("hist")
h_1.Draw("same hist")
h_2.Draw("same hist")

print h_0.GetEntries(),h_1.GetEntries(),h_2.GetEntries()

leg = ROOT.TLegend(0.45,0.90,0.25,0.80)
leg.SetBorderSize(0)
leg.SetTextSize(0.03)
leg.AddEntry(h_0,"H->WW","f")
leg.AddEntry(h_1,"H->ZZ->2l2x","f")
leg.AddEntry(h_2,"H->ZZ->4l","f")
leg.Draw("same")

c1.SaveAs("truth.png")
c1.SaveAs("truth.svg")
# c1.SaveAs("BDT.pdf")
c1.Clear()
    #     if "25" not in str(vp):
    #         flag += 1
    #         wl_id = j
    # if flag > 1: continue

    # ll = [0,1,2].remove(wl_id)
    # lep[0] = t.v_tlv_L[ll[0]]
    # lep[1] = t.v_tlv_L[ll[1]]
    # lep[2] = t.v_tlv_L[wl_id]

    # print wl_id,t.event
            # print vp