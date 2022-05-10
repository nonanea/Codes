import os
import sys
import ROOT
from math import fabs,pow,sqrt
from array import array
from operator import itemgetter

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

c1 = ROOT.TCanvas("c1", "canvas", 1000, 800)

c1.SetLeftMargin(0.15)
# c1.SetRightMargin(0.05)
c1.SetRightMargin(0.15)
c1.SetTopMargin(0.05)
c1.SetBottomMargin(0.15)

Z_mass = 91.1876

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

# input variables initialization

# p4 = ROOT.TLorentzVector(1,1,0,5)
f_xsec = "Xsec/3muv.txt"
with open(f_xsec, "r") as f:
    xsec_list = f.readlines()

W = ROOT.TLorentzVector(1,1,0,5)
Zp = ROOT.TLorentzVector(1,1,0,5)
lep = [ROOT.TLorentzVector(1,1,0,5) for i in range(3)]
nv = ROOT.TLorentzVector(1,1,0,5)
jet = [ROOT.TLorentzVector(1,1,0,5) for i in range(3)]

lep_id = [array('i',[0]) for i in range(3)]
Zpp_id = array('i',[0])
vp_id = array('i',[0])
# lep3_id = array('i',[0])

dir = "./High_Mass/"
# dir = "./g_test/"

input_f_list = os.listdir(dir)

plot_name = [
            "M_W", #0
            "M_Zp", #1
            "M_ll1", #2
            "M_ll2", #3
            "MET", #4
            "Pt_1", #5
            "Pt_2", #6
            "Pt_3", #7
            "Eta_1", #8
            "Eta_2", #9
            "Eta_3", #10
            "Phi_1", #11
            "Phi_2", #12
            "Phi_3", #13
            "ID_1", #14
            "ID_2", #15
            "ID_3", #16
            "ID_W", #17
            "JID_1", #18
            "JID_2", #19
]

h_name = [
            "M_{W} [GeV]", #0
            "M_{Z'} [GeV]", #1
            "M_{ll,1} [GeV]", #2
            "M_{ll,2} [GeV]", #3
            "MET [GeV]", #4
            "P_{T,1} [GeV]", #5
            "P_{T,2} [GeV]", #6
            "P_{T,3} [GeV]", #7
            "Eta_{1}", #8
            "Eta_{2}", #9
            "Eta_{3}", #10
            "Phi_{1}", #11
            "Phi_{2}", #12
            "Phi_{3}", #13
            "PID_{l,1}", #14
            "PID_{l,2}", #15
            "PID_{l,3}", #16
            "PID_{W}", #17
            "PID_{j,1}", #18
            "PID_{j,2}", #19
            ]

xrange = [
    (0,300), #0
    (0,120), #1
    (0,120), #2
    (0,120), #3
    (0,120), #4
    (0,200), #5
    (0,120), #6
    (0,80), #7
    
    # (0,550), #1
    # (0,550), #2
    # (0,550), #3
    # (0,250), #4
    # (0,300), #5
    # (0,200), #6
    # (0,100), #7

    (-4,4), #8
    (-4,4), #9
    (-4,4), #10
    (-4,4), #11
    (-4,4), #12
    (-4,4), #13
    (-15,15), #14
    (-15,15), #15
    (-15,15), #16
    (-30,30), #17
    (-6,6), #18
    (-6,6), #19
]

nbin = [
    40, #0
    40, #1
    40, #2
    40, #3
    40, #4
    40, #5
    40, #6
    40, #7
    40, #8
    40, #9
    40, #10
    40, #11
    40, #12
    40, #13
    30, #14
    30, #15
    30, #16
    30, #17
    12, #17
    12, #17
]

h = [ROOT.TH1F() for i in range(len(plot_name))]

for i in range(len(h_name)):
    # h_name = "hist" + str(i)
    h[i] = ROOT.TH1F(plot_name[i],plot_name[i],nbin[i],xrange[i][0],xrange[i][1])
    # h[i].SetLineColor(colorList[i])
    h[i].SetLineColor(ROOT.kRed)

h_2d = ROOT.TH2D("2d","2d",5,0,5,30,0,300)

for input_f in input_f_list:
    # if "muvZp005" not in input_f and "muvZp017" not in input_f and "muvZp027" not in input_f and "muvZp039" not in input_f and "muvZp051" not in input_f and "muvZp060" not in input_f and "muvZp075" not in input_f:
    #     continue

    # if "mass_400.root" not in input_f:
    if "05_" not in input_f:
        continue

    # if "100.root" not in input_f and "150.root" not in input_f and "250.root" not in input_f and "350.root" not in input_f and "500.root" not in input_f:
    # if not input_f == "01_mass_200.root":
    #     continue

    # if "muv" in input_f or "data" in input_f or "16a" in input_f:
        # continue

    # print input_f

    f_in = ROOT.TFile(dir+input_f)
    t = f_in.Get("t")

    t.SetBranchAddress("lep1_id",lep_id[0])
    t.SetBranchAddress("lep2_id",lep_id[1])
    t.SetBranchAddress("lep3_id",lep_id[2])
    
    # t.SetBranchAddress("Zpp_id",Zpp_id)
    # t.SetBranchAddress("vp_id",vp_id)

    # t.SetBranchAddress("lep1_id",lep1_id)
    # t.SetBranchAddress("lep2_id",lep2_id)
    # t.SetBranchAddress("lep3_id",lep3_id)

    t.SetBranchAddress("W_p4",W)
    t.SetBranchAddress("Zp_p4",Zp)
    t.SetBranchAddress("lep1_p4",lep[0])
    t.SetBranchAddress("lep2_p4",lep[1])
    t.SetBranchAddress("lep3_p4",lep[2])
    t.SetBranchAddress("vl_p4",nv)
    t.SetBranchAddress("jet1_p4",jet[0])
    t.SetBranchAddress("jet2_p4",jet[1])
    t.SetBranchAddress("jet3_p4",jet[2])

    # output variables initialization

    lep0_pt = array('f',[0.])
    lep1_pt = array('f',[0.])
    lep2_pt = array('f',[0.])
    lep0_eta = array('f',[0.])
    lep1_eta = array('f',[0.])
    lep2_eta = array('f',[0.])
    lep0_phi = array('f',[0.])
    lep1_phi = array('f',[0.])
    lep2_phi = array('f',[0.])
    lep0_E = array('f',[0.])
    lep1_E = array('f',[0.])
    lep2_E = array('f',[0.])
    lep0_id = array('i',[0])
    lep1_id = array('i',[0])
    lep2_id = array('i',[0])

    lep0_isoLoose = array('i',[1])
    lep1_isoLoose = array('i',[1])
    lep2_isoLoose = array('i',[1])

    lep0_isoTight = array('i',[1])
    lep1_isoTight = array('i',[1])
    lep2_isoTight = array('i',[1])

    met_met = array('f',[0.])
    met_px = array('f',[0.])
    met_py = array('f',[0.])

    mll_1 = array('f',[0.])
    mll_2 = array('f',[0.])
    mll_Z1 = array('f',[0.])
    mll_Z2 = array('f',[0.])

    mT = array('f',[0.])
    mT_vl = array('f',[0.])
    mT_Wvl = array('f',[0.])
    
    dR_1 = array('f',[0.])
    dR_2 = array('f',[0.])
    
    m_Zp = array('f',[0.])
    m_W = array('f',[0.])

    n_jets = array('i',[0])
    n_bjets = array('i',[0])

    weight = array('f',[0.])
    weight_g = array('f',[0.])
    xsec = array('f',[0.])

    f_out = ROOT.TFile("Outputs/"+input_f,"recreate")
    tt = ROOT.TTree("nominal","nominal")

    tt.Branch('lep0_pt',lep0_pt,'lep0_pt/F')
    tt.Branch('lep1_pt',lep1_pt,'lep1_pt/F')
    tt.Branch('lep2_pt',lep2_pt,'lep2_pt/F')
    tt.Branch('lep0_eta',lep0_eta,'lep0_eta/F')
    tt.Branch('lep1_eta',lep1_eta,'lep1_eta/F')
    tt.Branch('lep2_eta',lep2_eta,'lep2_eta/F')
    tt.Branch('lep0_phi',lep0_phi,'lep0_phi/F')
    tt.Branch('lep1_phi',lep1_phi,'lep1_phi/F')
    tt.Branch('lep2_phi',lep2_phi,'lep2_phi/F')
    tt.Branch('lep0_E',lep0_E,'lep0_E/F')
    tt.Branch('lep1_E',lep1_E,'lep1_E/F')
    tt.Branch('lep2_E',lep2_E,'lep2_E/F')
    tt.Branch('lep0_id',lep0_id,'lep0_id/I')
    tt.Branch('lep1_id',lep1_id,'lep1_id/I')
    tt.Branch('lep2_id',lep2_id,'lep2_id/I')
    
    tt.Branch('lep0_isoLoose',lep0_isoLoose,'lep0_isoLoose/I')
    tt.Branch('lep1_isoLoose',lep1_isoLoose,'lep1_isoLoose/I')
    tt.Branch('lep2_isoLoose',lep2_isoLoose,'lep2_isoLoose/I')
    
    tt.Branch('lep0_isoTight',lep0_isoTight,'lep0_isoTight/I')
    tt.Branch('lep1_isoTight',lep1_isoTight,'lep1_isoTight/I')
    tt.Branch('lep2_isoTight',lep2_isoTight,'lep2_isoTight/I')

    tt.Branch('met_met',met_met,'met_met_f/F')
    tt.Branch('met_px',met_px,'met_px_f/F')
    tt.Branch('met_py',met_py,'met_py_f/F')

    tt.Branch('mll_1',mll_1,'mll_1/F')
    tt.Branch('mll_2',mll_2,'mll_2/F')
    tt.Branch('mll_Z1',mll_Z1,'mll_Z1/F')
    tt.Branch('mll_Z2',mll_Z2,'mll_Z2/F')

    tt.Branch('mT',mT,'mT/F')
    tt.Branch('mT_vl',mT_vl,'mT_vl/F')
    tt.Branch('mT_Wvl',mT_Wvl,'mT_Wvl/F')

    tt.Branch('dR_1',dR_1,'dR_1/F')
    tt.Branch('dR_2',dR_2,'dR_2/F')
    
    tt.Branch('m_Zp',m_Zp,'m_Zp/F')
    tt.Branch('m_W',m_W,'m_W/F')

    tt.Branch('n_jets',n_jets,'n_jets/I')
    tt.Branch('n_bjets',n_bjets,'n_bjets/I')

    tt.Branch('weight',weight,'weight/F')
    tt.Branch('weight_g',weight_g,'weight_g/F')
    tt.Branch('xsection',xsec,'xsec/F')

    xsec_flag = 0
    weight[0] = 1.0/t.GetEntries()

    for xsec_value in xsec_list:
        # print float(input_f[input_f.find("mass")+5:-5])
        if float(input_f[input_f.find("mass")+5:-5]) == float(xsec_value.split()[0]):
            print "xsec:",xsec_value.split()[2],",coupling:",xsec_value.split()[1]
            xsec[0] = float(xsec_value.split()[2])*1e6
            weight_g[0] = xsec[0]/float(xsec_value.split()[1])/float(xsec_value.split()[1])/t.GetEntries()
            # print weight_g[0]
            xsec_flag = 1

    # if xsec_flag == 0:
    #     print "Xsec not found!"
    #     xsec[0] = 1.0
            
    # f_out.Close()

    # weight[0] = 1.0*xsec[0]/t.GetEntries()

    j = 0
    k = 0

    # print "Total number:",t.GetEntries()
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        # if i > 50:
        #     break

        l_order = sorted(enumerate([lep[0].Pt(),lep[1].Pt(),lep[2].Pt()]),key=itemgetter(1))

        lep0 = lep[l_order[2][0]]
        lep1 = lep[l_order[1][0]]
        lep2 = lep[l_order[0][0]]

        if fabs(lep0.Eta()) > 2.5 or fabs(lep1.Eta()) > 2.5 or fabs(lep2.Eta()) > 2.5:
            continue

        if lep2.Pt() < 3:
            continue

        # if lep0.Pt() < 20:
        #     continue

        # if abs(lep_id[0][0]+lep_id[1][0]+lep_id[2][0]) != 13:
        #     continue

        # print(lep_id[0][0],lep_id[1][0],lep_id[2][0])

        # print(lep[l_order[0][0]].Pt(),lep[l_order[1][0]].Pt(),lep[l_order[2][0]].Pt())

        # if lep_id[0][0] != lep_id[1][0]:
        #     if lep_id[0][0] != lep_id[2][0]:
        #         if (lep[0]+lep[1]).M() > (lep[0]+lep[2]).M():
        #             ll_1 = lep[0] + lep[1]
        #             ll_2 = lep[0] + lep[2]
        #         else:
        #             ll_2 = lep[0] + lep[1]
        #             ll_1 = lep[0] + lep[2]

        #     else:
        #         if (lep[0]+lep[1]).M() > (lep[1]+lep[2]).M():
        #             ll_1 = lep[0] + lep[1]
        #             ll_2 = lep[1] + lep[2]
        #         else:
        #             ll_2 = lep[0] + lep[1]
        #             ll_1 = lep[1] + lep[2]
        # else:
        #     if (lep[0]+lep[2]).M() > (lep[1]+lep[2]).M():
        #         ll_1 = lep[0] + lep[2]
        #         ll_2 = lep[1] + lep[2]
        #     else:
        #         ll_2 = lep[0] + lep[2]
        #         ll_1 = lep[1] + lep[2]
        
        if lep_id[0][0]*lep_id[1][0] == lep_id[0][0]*lep_id[2][0]:
            os = 0
        else:
            os = 2 if lep_id[0][0]*lep_id[1][0] > lep_id[0][0]*lep_id[2][0] else 1

        ss_1 = 0.0
        ss_2 = -1.0

        for k in range(3):
            if k != os:
                ss_1 = k
                if ss_2 < 0:
                    ss_2 = ss_1

        ll_1 = lep[os] + lep[ss_1]
        ll_2 = lep[os] + lep[ss_2]
        
        if fabs(ll_1.M()-Z_mass) < fabs(ll_2.M()-Z_mass):
        # if sqrt(pow(lep[ss_1].Phi()-lep[os].Phi(),2)+pow(lep[ss_1].Eta()-lep[os].Eta(),2)) > sqrt(pow(lep[ss_2].Phi()-lep[os].Phi(),2)+pow(lep[ss_2].Eta()-lep[os].Eta(),2)):
            mll_Z1[0] = ll_1.M()
            mll_Z2[0] = ll_2.M()

            mT_Wvl[0] = pow(lep[ss_2].Et() + nv.Et(), 2) - pow( lep[ss_2].Px() + nv.Px() , 2) - pow ( lep[ss_2].Py() + nv.Py() , 2)
            if mT_Wvl[0] > 0:   mT_Wvl[0] = sqrt(mT_Wvl[0])
            else:   mT_Wvl[0] = 0
            # dR_1[0] = sqrt(pow(lep[ss_1].Phi()-lep[os].Phi(),2)+pow(lep[ss_1].Eta()-lep[os].Eta(),2))
            # dR_2[0] = sqrt(pow(lep[ss_2].Phi()-lep[os].Phi(),2)+pow(lep[ss_2].Eta()-lep[os].Eta(),2))
        else:
            mll_Z2[0] = ll_1.M()
            mll_Z1[0] = ll_2.M()

            mT_Wvl[0] = pow(lep[ss_1].Et() + nv.Et(), 2) - pow( lep[ss_1].Px() + nv.Px() , 2) - pow ( lep[ss_1].Py() + nv.Py() , 2)
            if mT_Wvl[0] > 0:   mT_Wvl[0] = sqrt(mT_Wvl[0])
            else:   mT_Wvl[0] = 0
        
        if ll_1.M() > ll_2.M():
        # if sqrt(pow(lep[ss_1].Phi()-lep[os].Phi(),2)+pow(lep[ss_1].Eta()-lep[os].Eta(),2)) > sqrt(pow(lep[ss_2].Phi()-lep[os].Phi(),2)+pow(lep[ss_2].Eta()-lep[os].Eta(),2)):
            mll_1[0] = ll_1.M()
            mll_2[0] = ll_2.M()

            mT_vl[0] = pow(lep[ss_2].Et() + nv.Et(), 2) - pow( lep[ss_2].Px() + nv.Px() , 2) - pow ( lep[ss_2].Py() + nv.Py() , 2)
            # mT_vl[0] = sqrt(mT_vl[0])
            if mT_vl[0] > 0:   mT_vl[0] = sqrt(mT_vl[0])
            else:   mT_vl[0] = 0
            # dR_1[0] = sqrt(pow(lep[ss_1].Phi()-lep[os].Phi(),2)+pow(lep[ss_1].Eta()-lep[os].Eta(),2))
            # dR_2[0] = sqrt(pow(lep[ss_2].Phi()-lep[os].Phi(),2)+pow(lep[ss_2].Eta()-lep[os].Eta(),2))
        else:
            mll_2[0] = ll_1.M()
            mll_1[0] = ll_2.M()

            mT_vl[0] = pow(lep[ss_1].Et() + nv.Et(), 2) - pow( lep[ss_1].Px() + nv.Px() , 2) - pow ( lep[ss_1].Py() + nv.Py() , 2)
            if mT_vl[0] > 0:   mT_vl[0] = sqrt(mT_vl[0])
            else:   mT_vl[0] = 0
            # dR_2[0] = sqrt(pow(lep[ss_1].Phi()-lep[os].Phi(),2)+pow(lep[ss_1].Eta()-lep[os].Eta(),2))
            # dR_1[0] = sqrt(pow(lep[ss_2].Phi()-lep[os].Phi(),2)+pow(lep[ss_2].Eta()-lep[os].Eta(),2))

        # if ll_2.M() < 4:
        #     continue
            
        # if mll_2[0] < 4:
        #     continue

        dR = 0
        dRR = -1.0

        for k in range(3):
            if k != os:
                # dR = sqrt(pow(lep[k].Phi()*lep[k].Phi()-lep[os].Phi()*lep[os].Phi(),2)+pow(lep[k].Eta()*lep[k].Eta()-lep[os].Eta()*lep[os].Eta(),2))
                dR = sqrt(pow(lep[k].Phi()-lep[os].Phi(),2)+pow(lep[k].Eta()-lep[os].Eta(),2))
                if dRR < 0:
                    dRR = dR
        
        dR_1[0] = dR if dR >= dRR else dRR
        dR_2[0] = dR if dRR >= dR else dRR

        # h[0].Fill(ll_1.M())
        # h[1].Fill(ll_2.M())
        # h[2].Fill(lep1.Pt())
        # h[3].Fill(lep2.Pt())
        # h[4].Fill(lep3.Pt())

        # fill branch

        lep0_pt[0] = lep0.Pt()
        lep1_pt[0] = lep1.Pt()
        lep2_pt[0] = lep2.Pt()
        lep0_eta[0] = lep0.Eta()
        lep1_eta[0] = lep1.Eta()
        lep2_eta[0] = lep2.Eta()
        lep0_phi[0] = lep0.Phi()
        lep1_phi[0] = lep1.Phi()
        lep2_phi[0] = lep2.Phi()
        lep0_E[0] = lep0.E()
        lep1_E[0] = lep1.E()
        lep2_E[0] = lep2.E()
        lep0_id[0] = lep_id[l_order[2][0]][0]
        lep1_id[0] = lep_id[l_order[1][0]][0]
        lep2_id[0] = lep_id[l_order[0][0]][0]

        met_met[0] = nv.Et()
        met_px[0] = nv.Px()
        met_py[0] = nv.Py()

        mts = pow(lep[0].Et() + lep[1].Et() + lep[2].Et() + nv.Et(), 2) - pow( lep[0].Px() + lep[1].Px() + lep[2].Px() + nv.Px() , 2) - pow ( lep[0].Py() + lep[1].Py() + lep[2].Py() + nv.Py() , 2)
        mT[0] = sqrt(mts)

        # mll_1[0] = ll_1.M()
        # mll_2[0] = ll_2.M()

        m_Zp[0] = Zp.M()
        m_W[0] = (lep0+lep1+lep2+nv).M()

        n_jets[0] = 0
        n_bjets[0] = 0

        # xsec[0]

        # if Zpp_id[0] > 2.5:
        #     continue

        tt.Fill()
        h[0].Fill((lep0+lep1+lep2+nv).M())
        h[1].Fill(Zp.M())
        h[2].Fill(ll_1.M())
        h[3].Fill(ll_2.M())
        h[4].Fill(nv.Et())
        h[5].Fill(lep0.Pt())
        h[6].Fill(lep1.Pt())
        h[7].Fill(lep2.Pt())
        h[8].Fill(lep0.Eta())
        h[9].Fill(lep1.Eta())
        h[10].Fill(lep2.Eta())
        h[11].Fill(lep0.Phi())
        h[12].Fill(lep1.Phi())
        h[13].Fill(lep2.Phi())
        h[14].Fill(lep_id[0][0])
        h[15].Fill(lep_id[1][0])
        h[16].Fill(lep_id[2][0])
        # h[14].Fill(lep_id[l_order[2][0]][0])
        # h[15].Fill(lep_id[l_order[1][0]][0])
        # h[16].Fill(lep_id[l_order[0][0]][0])
        h[17].Fill(t.W_id)
        h[18].Fill(t.jet1_id)
        h[19].Fill(t.jet2_id)

        h_2d.Fill(Zpp_id[0],Zp.M())
        j+=1

    tt.Write()
    # f_out.Close()
    # print "Final number:",j
    print input_f,t.GetEntries(),j,j*1.0/t.GetEntries()

    # c1.SetLogy()
    for i in range(len(plot_name)):
        # h_name = "hist" + str(i)
        h[i].GetXaxis().SetTitle(h_name[i])
        h[i].GetYaxis().SetTitle("Events")
        # h[i].GetYaxis().SetRangeUser(1e-1,h[i].GetMaximum()*200)
        h[i].Draw("hist")
        c1.SaveAs("Plots/" + input_f.replace(".root","_") + plot_name[i] + ".png")
        c1.SaveAs("Plots/" + input_f.replace(".root","_") + plot_name[i] + ".svg")
        c1.Clear()
        h[i].Reset("ICESM")


h_2d.GetXaxis().SetTitle("Parent id")
h_2d.GetYaxis().SetTitle("M_{Z'}")
# h_2d.GetXaxis().SetNdivisions(505, ROOT.kTRUE)
# h_2d.SetMarkerSize(0.2)
h_2d.Draw("COLZ")

c1.SaveAs("Plots/2d.png")
c1.SaveAs("Plots/2d.svg")
c1.Clear()