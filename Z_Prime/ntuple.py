import ROOT
import os
from math import fabs,pow,sqrt
from array import array
from operator import itemgetter

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

c1 = ROOT.TCanvas("c1", "canvas", 1000, 800)

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.05)
c1.SetTopMargin(0.05)
c1.SetBottomMargin(0.15)

# input variables initialization

# p4 = ROOT.TLorentzVector(1,1,0,5)

lep_id = [array('f',[0]) for i in range(3)]
lep_pt = [array('f',[0]) for i in range(3)]
lep_eta = [array('f',[0]) for i in range(3)]
lep_phi = [array('f',[0]) for i in range(3)]
lep_m = [array('f',[0]) for i in range(3)]
# lep_isoLoose = [array('i',[0]) for i in range(3)]
# lep_isoTight = [array('i',[0]) for i in range(3)]

# met_met = array('f',[0.])
# met_px = array('f',[0.])
# met_py = array('f',[0.])

weight = array('f',[0.])
xsec = array('f',[0.])

input_f_list = os.listdir(".")

for input_f in input_f_list:
    if "361" not in input_f:
    # if "data" not in input_f and "_16" not in input_f:
        continue

    print(input_f + " processing")

    f_in = ROOT.TFile(input_f)
    t = f_in.Get("tree_PFLOW")
    hInfo = f_in.Get("Hist/hInfo_PFlow")

    # t.SetBranchAddress("lep0_pt",lep_pt[0])
    # t.SetBranchAddress("lep1_pt",lep_pt[1])
    # t.SetBranchAddress("lep2_pt",lep_pt[2])
    # t.SetBranchAddress("lep0_eta",lep_eta[0])
    # t.SetBranchAddress("lep1_eta",lep_eta[1])
    # t.SetBranchAddress("lep2_eta",lep_eta[2])
    # t.SetBranchAddress("lep0_phi",lep_phi[0])
    # t.SetBranchAddress("lep1_phi",lep_phi[1])
    # t.SetBranchAddress("lep2_phi",lep_phi[2])
    # t.SetBranchAddress("lep0_m",lep_m[0])
    # t.SetBranchAddress("lep1_m",lep_m[1])
    # t.SetBranchAddress("lep2_m",lep_m[2])
    t.SetBranchAddress("lep0_charge",lep_id[0])
    t.SetBranchAddress("lep1_charge",lep_id[1])
    t.SetBranchAddress("lep2_charge",lep_id[2])

    # t.SetBranchAddress("lep0_isoLoose",lep_isoLoose[0])
    # t.SetBranchAddress("lep1_isoLoose",lep_isoLoose[1])
    # t.SetBranchAddress("lep2_isoLoose",lep_isoLoose[2])
    
    # t.SetBranchAddress("lep0_isoTight",lep_isoTight[0])
    # t.SetBranchAddress("lep1_isoTight",lep_isoTight[1])
    # t.SetBranchAddress("lep2_isoTight",lep_isoTight[2])
    
    # t.SetBranchAddress("met_tst",met_met)
    # t.SetBranchAddress("met_px_tst",met_px)
    # t.SetBranchAddress("met_py_tst",met_py)

    # t.SetBranchAddress("weight",weight)
    # t.SetBranchAddress("weight_xsection",xsec)

    h_name = [
                "M_{ll,1}",
                "M_{ll,2}",
                "P_{T,1}",
                "P_{T,2}",
                "P_{T,3}",
                # "M_{T,l1v}",
                # "M_{T,l2v}",
                # "M_{T,l3v}",
                ]

    h = [ROOT.TH1F() for i in range(len(h_name))]

    # for i in range(len(h_name)):
    #     # h_name = "hist" + str(i)
    #     h[i] = ROOT.TH1F(h_name[i],h_name[i],20,0,80)
    #     # h[i].SetLineColor(i+2)

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

    lep0_isoLoose = array('i',[0])
    lep1_isoLoose = array('i',[0])
    lep2_isoLoose = array('i',[0])

    lep0_isoTight = array('i',[0])
    lep1_isoTight = array('i',[0])
    lep2_isoTight = array('i',[0])

    met_met = array('f',[0.])
    # met_phi = array('f',[0.])
    met_px = array('f',[0.])
    met_py = array('f',[0.])

    mT = array('f',[0.])
    mll_1 = array('f',[0.])
    mll_2 = array('f',[0.])

    n_jets = array('i',[0])
    n_bjets = array('i',[0])

    weight = array('f',[0.])
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

    tt.Branch('met_met',met_met,'met_met/F')
    tt.Branch('met_px',met_px,'met_px/F')
    tt.Branch('met_py',met_py,'met_py/F')

    tt.Branch('mT',mT,'mT/F')
    tt.Branch('mll_1',mll_1,'mll_1/F')
    tt.Branch('mll_2',mll_2,'mll_2/F')

    tt.Branch('n_jets',n_jets,'n_jets/I')
    tt.Branch('n_bjets',n_bjets,'n_bjets/I')

    tt.Branch('weight',weight,'weight/F')
    tt.Branch('xsection',xsec,'xsec/F')

    # f_out.Close()
        
    xsec[0] = t.weight_xsection

    # print(t.GetEntries())
    
    lep = [ROOT.TLorentzVector(1,1,0,5) for i in range(3)]

    j = 0

    for i in range(t.GetEntries()):
        t.GetEntry(i)
        # if i > 50:
        #     break

        if t.lep0_FixedCutPflowLoose == 0 or t.lep1_FixedCutPflowLoose == 0 or  t.lep2_FixedCutPflowLoose == 0:
            continue

        # lep[0].SetPtEtaPhiM(lep_pt[0][0],lep_eta[0][0],lep_phi[0][0],lep_m[0][0])
        # lep[1].SetPtEtaPhiM(lep_pt[1][0],lep_eta[1][0],lep_phi[1][0],lep_m[1][0])
        # lep[2].SetPtEtaPhiM(lep_pt[2][0],lep_eta[2][0],lep_phi[2][0],lep_m[2][0])
        
        lep[0].SetPtEtaPhiM(t.lep0_pt,t.lep0_eta,t.lep0_phi,t.lep0_m)
        lep[1].SetPtEtaPhiM(t.lep1_pt,t.lep1_eta,t.lep1_phi,t.lep1_m)
        lep[2].SetPtEtaPhiM(t.lep2_pt,t.lep2_eta,t.lep2_phi,t.lep2_m)

        l_order = sorted(enumerate([lep[0].Pt(),lep[1].Pt(),lep[2].Pt()]),key=itemgetter(1))

        lep0 = lep[l_order[2][0]]
        lep1 = lep[l_order[1][0]]
        lep2 = lep[l_order[0][0]]

        if lep0.Pt() < 20:
            continue
            
        # print(lep_id[0][0],lep_id[1][0],lep_id[2][0])

        # print(lep[l_order[0][0]].Pt(),lep[l_order[1][0]].Pt(),lep[l_order[2][0]].Pt())

        if lep_id[0][0] != lep_id[1][0]:
            if lep_id[0][0] != lep_id[2][0]:
                if (lep[0]+lep[1]).M() > (lep[0]+lep[2]).M():
                    ll_1 = lep[0] + lep[1]
                    ll_2 = lep[0] + lep[2]
                else:
                    ll_2 = lep[0] + lep[1]
                    ll_1 = lep[0] + lep[2]
            else:
                if (lep[0]+lep[1]).M() > (lep[1]+lep[2]).M():
                    ll_1 = lep[0] + lep[1]
                    ll_2 = lep[1] + lep[2]
                else:
                    ll_2 = lep[0] + lep[1]
                    ll_1 = lep[1] + lep[2]
        else:
            if (lep[0]+lep[2]).M() > (lep[1]+lep[2]).M():
                ll_1 = lep[0] + lep[2]
                ll_2 = lep[1] + lep[2]
            else:
                ll_2 = lep[0] + lep[2]
                ll_1 = lep[1] + lep[2]
        
        if ll_2.M() < 4:
            continue

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
        lep0_id[0] = int(lep_id[l_order[2][0]][0])*13
        lep1_id[0] = int(lep_id[l_order[1][0]][0])*13
        lep2_id[0] = int(lep_id[l_order[0][0]][0])*13

        # print(lep1_pt[0])
        met_met[0] = t.met_tst
        met_px[0] = t.met_px_tst
        met_py[0] = t.met_py_tst

        lep0_isoLoose[0] = t.lep0_FixedCutPflowLoose
        lep1_isoLoose[0] = t.lep1_FixedCutPflowLoose
        lep2_isoLoose[0] = t.lep2_FixedCutPflowLoose
        
        lep0_isoTight[0] = t.lep0_FixedCutPflowTight
        lep1_isoTight[0] = t.lep1_FixedCutPflowTight
        lep2_isoTight[0] = t.lep2_FixedCutPflowTight

        mts = pow(lep[0].Et() + lep[1].Et() + lep[2].Et() + t.met_tst, 2) - pow( lep[0].Px() + lep[1].Px() + lep[2].Px() + t.met_px_tst , 2) - pow ( lep[0].Py() + lep[1].Py() + lep[2].Py() + t.met_py_tst , 2)
        mT[0] = sqrt(mts)

        mll_1[0] = ll_1.M()
        mll_2[0] = ll_2.M()

        n_jets[0] = t.n_jets
        n_bjets[0] = t.n_bjets

        if "data" in input_f:
            weight[0] = 1.0
        else:
            weight[0] = t.weight*hInfo.GetBinContent(1)*2.0/hInfo.GetEntries()/hInfo.GetBinContent(2)

        tt.Fill()
        j+=1

    tt.Write()
    f_out.Close()
    t.Delete()
    print "Event number:",j

# for i in range(len(h_name)):
#     # h_name = "hist" + str(i)
#     h[i].GetXaxis().SetTitle(h_name[i])
#     h[i].GetYaxis().SetTitle("Events")
#     h[i].Draw("hist")
#     c1.SaveAs("Plots/" + h_name[i] + ".png")
