import ROOT
import os
import time
from sympy import symbols, solve, evalf
import sympy
from math import fabs,pow,sqrt,isnan
from array import array
from operator import itemgetter

# def Angle(x1,y1,x2,y2):
#     angle = (x1*x2+y1*y2)/sqrt(x1*x1+y1*y1)/sqrt(x2*x2+y2*y2)
#     return(angle)

# def nv_pz(a,b,c):
#     x = symbols('x')
#     # root_list = solve( (a + sqrt( x**2 + b ) )**2 - x**2 - b - c, x )
#     # return root[0].evalf().as_real_imag()[0]
#     # print solve( (a + sympy.sqrt( x*x + b ) )*(a + sympy.sqrt( x*x + b ) ) - x*x - b - c, x )[0].evalf().as_real_imag()[0]
#     # print solve( (a + sympy.sqrt( x*x + b ) )**2 - x*x - b - c, x )
#     print a*a,b,c
#     return

#     root_list = solve( (a + sqrt( x**2 + b ) )**2 - x**2 - b - c, x )
    # for root in root_list:
    # if root[0].evalf().as_real_imag()[0] < root[1].evalf().as_real_imag()[0]: return root[0].evalf().as_real_imag()[0]
    # else: return root[1].evalf().as_real_imag()[0]

op = time.time()

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

lumi = 1.0
W_mass = 80.379
Z_mass = 91.1876
Zp_mass = [ 5, 9, 15, 19, 23, 27, 31, 35, 39, 45, 51, 54, 60, 66, 69, 75 ]
# Zp_mass = [ 5, 9, 15, 19, 23, 27, 31, 35, 39 ]
# Zp_mass = [ 39, 45, 51, 54, 60, 66, 69, 75 ]

WZ_list = [364253,364284,363358]
ZZ_list = [345666,345723,364285,364250,364283,345283,245706,363356,364254]
WW_list = [345718]
VVV_list = [364242,364243,364244,364245,364246,364247,364248,364249]
Zjets_list = [361107,361666,361667,361106,361664,361665]
# Zjets_list = [364100,364101,364102,364103,364104,364105,364106,364107,364108,364109,364110,364111,364112,364113,364198,364199,364200,364201,364202,364203]
t_list = [410644,410645,410658,410659]
Wt_list = [410646,410647]
tt_list = [410472]
ttV_list = [410155,410156,410157,410081,410218,410219]
Zgamma_list = [364500,364501,364502,364503,364504,364505,364506,364507,364508,364509]
# extra_list = [345705,345706,345714,364283,364284,364286,364287,364288,364289,364290,361108]

f_xsec = "Xsec/3muv.txt"
with open(f_xsec, "r") as f:
    xsec_list = f.readlines()

f = open("Plots/temp.csv","a")
f.seek(0)
f.truncate()

# input variables initialization

# p4 = ROOT.TLorentzVector(1,1,0,5)

lep_cid = [array('f',[0]) for i in range(3)]
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

sample_list = os.listdir(".")

# h_mll = ROOT.TH1F("m_ll","m_ll",500,0,100)
gaus = ROOT.TF1("gaus","gaus(0)+gaus(3)",0,100)
gaus.SetLineColor(ROOT.kRed)

bin = [0.8,2.1,3,3,5,5,5,5,7,7,7,7,12,12,12,12]

h_list = [[ROOT.TH1F()] for i in range(len(Zp_mass))]
for i in range(len(h_list)):
    h_list[i] = ROOT.TH1F(str(Zp_mass[i]),str(Zp_mass[i]),int(2500/bin[i]),0,100)

for sample in sample_list:
    if "muv" in sample:   continue
    if "root" not in sample:
    # if "361107" not in sample and "361666" not in sample and "361667" not in sample :
    # if "36428" not in sample:
        continue

    # if "Loose" not in sample and "LowPt" not in sample and "Med" not in sample:    continue
    # if "data" in sample or "muv" in sample: continue
    # if int(sample[:6]) not in Zjets_list and int(sample[:6]) not in Zgamma_list:   continue

    print(sample + " processing")

    f_in = ROOT.TFile(sample)
    t = f_in.Get("tree_PFLOW")
    hInfo = f_in.Get("Hist/hInfo_PFlow")

    t.SetBranchAddress("lep0_charge",lep_cid[0])
    t.SetBranchAddress("lep1_charge",lep_cid[1])
    t.SetBranchAddress("lep2_charge",lep_cid[2])
    
    t.SetBranchAddress("lep0_id",lep_id[0])
    t.SetBranchAddress("lep1_id",lep_id[1])
    t.SetBranchAddress("lep2_id",lep_id[2])

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

    for i in range(len(h_name)):
        # h_name = "hist" + str(i)
        h[i] = ROOT.TH1F(h_name[i],h_name[i],20,0,80)
        # h[i].SetLineColor(i+2)

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
    lep0_cid = array('i',[0])
    lep1_cid = array('i',[0])
    lep2_cid = array('i',[0])

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
    met_signif = array('f',[0.])

    mll_1 = array('f',[0.])
    mll_2 = array('f',[0.])
    mll_Z1 = array('f',[0.])
    mll_Z2 = array('f',[0.])
    mlll = array('f',[0.])

    HT = array('f',[0.])
    LT = array('f',[0.])
    VT = array('f',[0.])

    mT = array('f',[0.])
    mT_WZ = array('f',[0.])
    mT_Wvl = array('f',[0.])
    mT_vl = array('f',[0.])
    
    dR_1 = array('f',[0.])
    dR_2 = array('f',[0.])
    
    dPhi_1 = array('f',[0.])
    dPhi_2 = array('f',[0.])

    n_jets = array('i',[0])
    n_bjets = array('i',[0])

    weight = array('f',[0.])
    weight_g = array('f',[0.])
    xsec = array('f',[0.])
    # lumi = array('f',[0.])
    
    mass = array('f',[0.])
    bkg = array('i',[0])

    f_out = ROOT.TFile("Outputs/"+sample,"recreate")
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
    # tt.Branch('lep0_cid',lep0_cid,'lep0_cid/I')
    # tt.Branch('lep1_cid',lep1_cid,'lep1_cid/I')
    # tt.Branch('lep2_cid',lep2_cid,'lep2_cid/I')
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
    tt.Branch('met_signif',met_signif,'met_signif/F')

    tt.Branch('mll_1',mll_1,'mll_1/F')
    tt.Branch('mll_2',mll_2,'mll_2/F')
    tt.Branch('mll_Z1',mll_Z1,'mll_Z1/F')
    tt.Branch('mll_Z2',mll_Z2,'mll_Z2/F')
    tt.Branch('mlll',mlll,'mlll/F')

    tt.Branch('HT',HT,'HT/F')
    tt.Branch('LT',LT,'LT/F')
    tt.Branch('VT',VT,'VT/F')

    tt.Branch('mT',mT,'mT/F')
    tt.Branch('mT_WZ',mT_WZ,'mT_WZ/F')
    tt.Branch('mT_Wvl',mT_Wvl,'mT_Wvl/F')
    tt.Branch('mT_vl',mT_vl,'mT_vl/F')

    tt.Branch('dR_1',dR_1,'dR_1/F')
    tt.Branch('dR_2',dR_2,'dR_2/F')

    tt.Branch('dPhi_1',dPhi_1,'dPhi_1/F')
    tt.Branch('dPhi_2',dPhi_2,'dPhi_2/F')

    tt.Branch('n_jets',n_jets,'n_jets/I')
    tt.Branch('n_bjets',n_bjets,'n_bjets/I')

    tt.Branch('weight',weight,'weight/F')
    tt.Branch('weight_g',weight_g,'weight_g/F')
    tt.Branch('xsection',xsec,'xsec/F')
    # tt.Branch('lumi',lumi,'lumi/F')

    tt.Branch('mass',mass,'mass/F')
    tt.Branch('bkg',bkg,'bkg/I')

    xsec_flag = 0
    
    if "mc16a" in sample:   lumi = 36.0746
    elif "mc16d" in sample: lumi = 43.8137
    elif "mc16e" in sample: lumi = 58.4501

    # if "data" not in sample and "muv" not in sample:
    #     if int(sample[:6]) in Zjets_list:    lumi = 139
    # lumi = 36.0746
    
    if "muv" in sample:
        mass[0] = float(sample[sample.find("muvZp")+5:sample.find("muvZp")+8])
        id = Zp_mass.index(mass[0])
        if mass[0] > 10: continue
        # print id
        for xsec_value in xsec_list:
            # print float(sample[sample.find("mass")+5:-5])
            if float(sample[sample.find("muvZp")+5:sample.find("muvZp")+8]) == float(xsec_value.split()[0]):
                print "xsec:",xsec_value.split()[2],",coupling:",xsec_value.split()[1]
                xsec[0] = float(xsec_value.split()[2])*1e6
                weight_g[0] = xsec[0]/float(xsec_value.split()[1])/float(xsec_value.split()[1])
                # print weight_g[0]
                xsec_flag = 1

    # f_out.Close()
    # if mass[0] > 19: continue

    # print(mass[0])
    
    lep = [ROOT.TLorentzVector(1,1,0,5) for i in range(3)]
    met = ROOT.TLorentzVector(1,1,0,5)

    j = 0
    k = 0
    notrigger = 0

    for i in range(t.GetEntries()):
        t.GetEntry(i)
        # if i > 10:    break
        
        # if isnan(t.weight): continue
        
        # if fabs(lep_id[0][0]+lep_id[1][0]+lep_id[2][0]) != 39:
        #     continue
        # if fabs(lep_id[0][0]) != 13 or fabs(lep_id[1][0]) != 13 or fabs(lep_id[2][0]) != 13:    continue
        #     # print(lep_id[0][0],lep_id[1][0],lep_id[2][0])

        if fabs(lep_cid[0][0]+lep_cid[1][0]+lep_cid[2][0]) > 1: continue

        if t.lep0_FixedCutPflowLoose == 0 or t.lep1_FixedCutPflowLoose == 0 or  t.lep2_FixedCutPflowLoose == 0: continue
            
        if t.lep0_eta > 2.5 or t.lep1_eta > 2.5 or t.lep2_eta > 2.5:   continue

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

        if lep0.Pt() < 20:    continue

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

        if lep_cid[0][0]*lep_cid[1][0] == lep_cid[0][0]*lep_cid[2][0]:
            os = 0
        else:
            os = 2 if lep_cid[0][0]*lep_cid[1][0] > lep_cid[0][0]*lep_cid[2][0] else 1

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
            #   Not Z muon: ss_2

            mll_Z1[0] = ll_1.M()
            mll_Z2[0] = ll_2.M()

            mT_Wvl[0] = pow(lep[ss_2].Et() + t.met_tst, 2) - pow( lep[ss_2].Px() + t.met_px_tst , 2) - pow ( lep[ss_2].Py() + t.met_py_tst , 2)
            if mT_Wvl[0] > 0:   mT_Wvl[0] = sqrt(mT_Wvl[0])
            else:   mT_Wvl[0] = 0

            # dR_1[0] = sqrt(pow(lep[ss_1].Phi()-lep[os].Phi(),2)+pow(lep[ss_1].Eta()-lep[os].Eta(),2))
            # dR_2[0] = sqrt(pow(lep[ss_2].Phi()-lep[os].Phi(),2)+pow(lep[ss_2].Eta()-lep[os].Eta(),2))
        else:
            #   Not Z muon: ss_1
            mll_Z2[0] = ll_1.M()
            mll_Z1[0] = ll_2.M()

            mT_Wvl[0] = pow(lep[ss_1].Et() + t.met_tst, 2) - pow( lep[ss_1].Px() + t.met_px_tst , 2) - pow ( lep[ss_1].Py() + t.met_py_tst , 2)
            if mT_Wvl[0] > 0:   mT_Wvl[0] = sqrt(mT_Wvl[0])
            else:   mT_Wvl[0] = 0

            # dR_2[0] = sqrt(pow(lep[ss_1].Phi()-lep[os].Phi(),2)+pow(lep[ss_1].Eta()-lep[os].Eta(),2))
            # dR_1[0] = sqrt(pow(lep[ss_2].Phi()-lep[os].Phi(),2)+pow(lep[ss_2].Eta()-lep[os].Eta(),2))
        
        if ll_1.M() > ll_2.M():
        # if sqrt(pow(lep[ss_1].Phi()-lep[os].Phi(),2)+pow(lep[ss_1].Eta()-lep[os].Eta(),2)) > sqrt(pow(lep[ss_2].Phi()-lep[os].Phi(),2)+pow(lep[ss_2].Eta()-lep[os].Eta(),2)):
            
            #   W muon: ss_2

            mll_1[0] = ll_1.M()
            mll_2[0] = ll_2.M()

            mT_vl[0] = pow(lep[ss_2].Et() + t.met_tst, 2) - pow( lep[ss_2].Px() + t.met_px_tst , 2) - pow ( lep[ss_2].Py() + t.met_py_tst , 2)
            if mT_vl[0] > 0:   mT_vl[0] = sqrt(mT_vl[0])
            else:   mT_vl[0] = 0

            # pow((lep[os]+lep[ss_1]).Px() + t.met_px_tst, 2) + pow((lep[os]+lep[ss_1]).Py() + t.met_py_tst, 2)
            if lep[ss_2].Et < t.met_tst:    LT[0] = ( lep[os] + lep[ss_1] + lep[ss_2] ).Pt()
            else: LT[0] = lep[ss_2].Pt()

            dPhi_1[0] = (lep[os]+lep[ss_1]).Phi() - lep[ss_2].Phi()
            dPhi_2[0] = (lep[os]+lep[ss_2]).Phi() - lep[ss_1].Phi()

            dR_1[0] = sqrt( pow( (lep[os]+lep[ss_1]).Phi()-lep[ss_2].Phi(),2 )+pow( (lep[os]+lep[ss_1]).Eta()-lep[ss_2].Eta(),2) )
            dR_2[0] = sqrt( pow( (lep[os]+lep[ss_2]).Phi()-lep[ss_1].Phi(),2 )+pow( (lep[os]+lep[ss_2]).Eta()-lep[ss_1].Eta(),2) )

            # dPhi_1[0] = Angle( (lep[os]+lep[ss_1]).Px() + t.met_px_tst, (lep[os]+lep[ss_1]).Py() + t.met_py_tst, lep[ss_2].Px(), lep[ss_2].Py())
            # dPhi_2[0] = Angle( (lep[os]+lep[ss_1]).Px() , (lep[os]+lep[ss_1]).Py()  , t.met_px_tst, t.met_py_tst)

            # dR_1[0] = sqrt(pow(lep[ss_1].Phi()-lep[os].Phi(),2)+pow(lep[ss_1].Eta()-lep[os].Eta(),2))
            # dR_2[0] = sqrt(pow(lep[ss_2].Phi()-lep[os].Phi(),2)+pow(lep[ss_2].Eta()-lep[os].Eta(),2))
        else:
            #   W muon: ss_1

            mll_2[0] = ll_1.M()
            mll_1[0] = ll_2.M()

            mT_vl[0] = pow(lep[ss_1].Et() + t.met_tst, 2) - pow( lep[ss_1].Px() + t.met_px_tst , 2) - pow ( lep[ss_1].Py() + t.met_py_tst , 2)
            if mT_vl[0] > 0:   mT_vl[0] = sqrt(mT_vl[0])
            else:   mT_vl[0] = 0

            if lep[ss_1].Et < t.met_tst:    LT[0] = ( lep[os] + lep[ss_1] + lep[ss_2] ).Pt()
            else: LT[0] = lep[ss_1].Pt()

            dPhi_1[0] = (lep[os]+lep[ss_2]).Phi() - lep[ss_1].Phi()
            dPhi_2[0] = (lep[os]+lep[ss_1]).Phi() - lep[ss_2].Phi()

            dR_1[0] = sqrt( pow( (lep[os]+lep[ss_2]).Phi()-lep[ss_1].Phi(),2 ) + pow( (lep[os]+lep[ss_2]).Eta()-lep[ss_1].Eta(),2) )
            dR_2[0] = sqrt( pow( (lep[os]+lep[ss_1]).Phi()-lep[ss_2].Phi(),2 ) + pow( (lep[os]+lep[ss_1]).Eta()-lep[ss_2].Eta(),2) )

            # dPhi_1[0] = Angle( (lep[os]+lep[ss_2]).Px() + t.met_px_tst, (lep[os]+lep[ss_2]).Py() + t.met_py_tst, lep[ss_1].Px(), lep[ss_1].Py())
            # dPhi_2[0] = Angle( (lep[os]+lep[ss_2]).Px() , (lep[os]+lep[ss_2]).Py()  , t.met_px_tst, t.met_py_tst)

        # if  mll_1[0] != mll_Z1[0]:   print mll_1[0]

        # mll = 0.0
        # mll_0 = -1.0

        # for k in range(3):
        #     if k != os:
        #         # if lep_id[os][0] == lep_id[k][0]:    print("same",os,k,lep_id[os][0])
        #         mll = (lep[k] + lep[os]).M()
        #         if mll_0 < 0:
        #             mll_0 = mll
        
        # mll_1[0] = (lep[os]+lep[ss_1]).M() if (lep[os]+lep[ss_1]).M() >= (lep[os]+lep[ss_2]).M() else (lep[os]+lep[ss_2]).M()

        # mll_2[0] = (lep[os]+lep[ss_1]).M() if (lep[os]+lep[ss_1]).M() <= (lep[os]+lep[ss_2]).M() else (lep[os]+lep[ss_2]).M()
        
        # if ll_2.M() < 4:
        #     continue
            
        if mll_2[0] < 4:    continue

        # dR = 0
        # dRR = -1.0

        # for k in range(3):
        #     if k != os:
        #         dR = sqrt(pow(lep[k].Phi()-lep[os].Phi(),2)+pow(lep[k].Eta()-lep[os].Eta(),2))
        #         if dRR < 0:
        #             dRR = dR
        
        # dR_1[0] = dR if dR >= dRR else dRR
        # dR_2[0] = dR if dRR >= dR else dRR

        # print nv_pz( (lep[0]+lep[1]+lep[2]).E(), t.met_tst*t.met_tst, pow( (lep[0]+lep[1]+lep[2]).E(), 2) - pow( (lep[0]+lep[1]+lep[2]).M(), 2 ) + W_mass*W_mass )
        # print nv_pz( (lep[0]+lep[1]+lep[2]).E(), t.met_tst*t.met_tst, pow( (lep[0]+lep[1]+lep[2]).E(), 2) - pow( (lep[0]+lep[1]+lep[2]).M(), 2 ) + W_mass*W_mass )
        # nv_pz(1.0,2.0,3.0)
        # print (t.met_tst*t.met_tst - t.met_px_tst*t.met_px_tst - t.met_py_tst*t.met_py_tst)/t.met_tst

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
        lep0_id[0] = int(lep_id[l_order[2][0]][0])*int(lep_cid[l_order[2][0]][0])
        lep1_id[0] = int(lep_id[l_order[1][0]][0])*int(lep_cid[l_order[1][0]][0])
        lep2_id[0] = int(lep_id[l_order[0][0]][0])*int(lep_cid[l_order[0][0]][0])

        # print(lep1_pt[0])
        met_met[0] = t.met_tst
        met_px[0] = t.met_px_tst
        met_py[0] = t.met_py_tst
        met_signif[0] = t.met_signif

        lep0_isoLoose[0] = t.lep0_FixedCutPflowLoose
        lep1_isoLoose[0] = t.lep1_FixedCutPflowLoose
        lep2_isoLoose[0] = t.lep2_FixedCutPflowLoose
        
        lep0_isoTight[0] = t.lep0_FixedCutPflowTight
        lep1_isoTight[0] = t.lep1_FixedCutPflowTight
        lep2_isoTight[0] = t.lep2_FixedCutPflowTight

        n_jets[0] = t.n_jets
        n_bjets[0] = t.n_bjets

        # mts = pow(lep[0].Et() + lep[1].Et() + lep[2].Et() + t.met_tst, 2) - pow( lep[0].Px() + lep[1].Px() + lep[2].Px() + t.met_px_tst , 2) - pow ( lep[0].Py() + lep[1].Py() + lep[2].Py() + t.met_py_tst , 2)

        # dPhi_2[0] = Angle( (lep[0]+lep[1]+lep[2]).Px() , (lep[0]+lep[1]+lep[2]).Py() , t.met_px_tst, t.met_py_tst)
        
        mlll[0] = (lep[0] + lep[1] + lep[2]).M()

        HT[0] = lep[0].Pt() + lep[1].Pt() + lep[2].Pt() + t.met_tst
        VT[0] = sqrt( pow( (lep[0]+lep[1]+lep[2]).Px() + t.met_px_tst, 2) + pow( (lep[0]+lep[1]+lep[2]).Py() + t.met_py_tst, 2) )
        # print t.met_px_tst,(lep[0]+lep[1]+lep[2]).Px(),t.met_py_tst,(lep[0]+lep[1]+lep[2]).Py()
        # print (lep[0]+lep[1]+lep[2]).Px(),lep[0].Px()+lep[1].Px()+lep[2].Px()

        met.SetPxPyPzE(t.met_px_tst,t.met_py_tst,0,t.met_tst)
        # mT[0] = sqrt( pow( (lep[0] + lep[1] + lep[2] + met).Et() , 2) - VT[0]*VT[0] )

        mT[0] = sqrt( pow( (lep[0] + lep[1] + lep[2]).Et() + t.met_tst, 2) - VT[0]*VT[0] )

        lep[0].SetPtEtaPhiM(t.lep0_pt,0,t.lep0_phi,t.lep0_m)
        lep[1].SetPtEtaPhiM(t.lep1_pt,0,t.lep1_phi,t.lep1_m)
        lep[2].SetPtEtaPhiM(t.lep2_pt,0,t.lep2_phi,t.lep2_m)

        mT_WZ[0] = (lep[0] + lep[1] + lep[2] + met).M()

        # print (lep[0] + lep[1] + lep[2]).Et(), lep[0].Et() + lep[1].Et() + lep[2].Et(), (lep[0] + lep[1] + lep[2] + met).Et(), (lep[0] + lep[1] + lep[2] + met).Et()
        # print mT[0], (lep[0] + lep[1] + lep[2] + met).M()
        # mT[0] = sqrt(mts)
        # if fabs(t.met_tst*t.met_tst-t.met_px_tst*t.met_px_tst-t.met_py_tst*t.met_py_tst)>0.01: print(t.met_tst*t.met_tst-t.met_px_tst*t.met_px_tst-t.met_py_tst*t.met_py_tst)

        # if isnan(t.weight_trig): triSF = 1.0
        # else: triSF = t.weight_trig

        if "data" in sample:
            bkg[0] = 100
            mass[0] = Zp_mass[i%len(Zp_mass)]
            weight[0] = 1.0
        elif "muv" in sample:
            bkg[0] = 99
            weight[0] = t.weight/hInfo.GetBinContent(2)*lumi
            
            if fabs( ll_1.M() - mass[0] ) <  fabs( ll_2.M() - mass[0] ):    h_list[id].Fill(ll_1.M())
            else:    h_list[id].Fill(ll_2.M())

            # h[1].Fill(ll_2.M())
            # weight[0] = triSF * t.weight_pileup * t.weight_gen * t.weight_jvt * t.weight_jets * t.weight_exp/hInfo.GetBinContent(2) #   no trigger SF
        elif "_mc16" in sample:
            # weight[0] = t.weight*hInfo.GetBinContent(1)*2.0/hInfo.GetEntries()/hInfo.GetBinContent(2)*lumi
            # print t.weight/hInfo.GetBinContent(2), t.weight_pileup * t.weight_gen * t.weight_jvt * t.weight_jets * t.weight_exp/hInfo.GetBinContent(2)
            weight[0] = t.weight_pileup * t.weight_gen * t.weight_jvt * t.weight_jets * t.weight_exp*hInfo.GetBinContent(1)*2.0/hInfo.GetEntries()/hInfo.GetBinContent(2) * lumi   #   no trigger SF
            xsec[0] = t.weight_xsection
            
            mass[0] = Zp_mass[i%len(Zp_mass)]

            if int(sample[:6]) in WZ_list:    bkg[0] = 0
            elif int(sample[:6]) in ZZ_list:    bkg[0] = 1
            # elif int(sample[:6]) in WW_list:    k = 2
            elif int(sample[:6]) in Zjets_list:    bkg[0] = 2
            # elif int(sample[:6]) in t_list:     k = 5
            # elif int(sample[:6]) in Wt_list:     k = 6
            elif int(sample[:6]) in tt_list:    bkg[0] = 3
            elif int(sample[:6]) in ttV_list:    bkg[0] = 4
            elif int(sample[:6]) in Zgamma_list:    bkg[0] = 5
            # elif int(sample[:6]) in extra_list:    k = 5
            else:    bkg[0] = 6
        else:
            bkg[0] = 99
            weight[0] = t.weight/hInfo.GetBinContent(2)*lumi
            # weight[0] = t.weight_pileup * t.weight_gen * t.weight_jvt * t.weight_jets * t.weight_exp/hInfo.GetBinContent(2)*lumi #   no trigger SF
            xsec[0] = t.weight_xsection


        # if lep1_pt[0] < 10 or lep2_pt[0] < 6 or t.n_jets > 0 or mll_1[0] > 80:   continue
        # if lep1_pt[0] < 10 or lep2_pt[0] < 6:   continue
        # if lep2.Pt() < 10:    continue
        notrigger+=1
        try:
            if not t.passTrigger:   continue
            tt.Fill()
        except:
            tt.Fill()
        
        # h[0].Fill(ll_1.M())
        # h[1].Fill(ll_2.M())
        # h[2].Fill(lep1.Pt())
        # h[3].Fill(lep2.Pt())
        # h[4].Fill(lep3.Pt())

        j+=1
        
    # print xsec[0]
    try:
        print "Total number:",t.GetEntries(),"Event number:",j,"Trig Eff:",j*1.0/notrigger
        f.write(sample + "," + str(j*1.0/t.GetEntries()) + "\n")
    except ZeroDivisionError:
        print "No events"
    # try:
    #     print "Total number:",t.GetEntries(),"Event number:",j,"Trig Eff:",j*1.0/notrigger
    #     f.write(sample.replace(".root","") + "," + str(j*1.0/t.GetEntries()) + "," + str(h_mll.GetMean()) + "," + str(h_mll.GetRMS()) + "," + str(h_mll.GetEntries()) + "\n")
    # except ZeroDivisionError:
    #     print "No events"

    # if "muv" in sample and Zp_mass[len(h_list)] != mass[0]:
    if "muv" in sample and "16e" in sample:
        # h_list[id].GetXaxis().SetRangeUser( h_list[id].GetMean() - 3.0*h_list[id].GetRMS(),h_list[id].GetMean() + 3.0*h_list[id].GetRMS() )
        h_list[id].GetXaxis().SetRangeUser( h_list[id].GetMean() - bin[id],h_list[id].GetMean() + bin[id] )
        h_list[id].GetXaxis().SetTitle("Reco Zp Mass [GeV]")
        h_list[id].GetYaxis().SetTitle("Events")
        h_list[id].Draw("hist")

        gaus.SetParameters(h_list[id].GetMaximum()/2.0,h_list[id].GetMean()+0.5,h_list[id].GetRMS(),h_list[id].GetMaximum()/2.0,h_list[id].GetMean()-0.5,h_list[id].GetRMS())

        # r = h_list[id].Fit(gaus,"QS","C", h_list[id].GetMean() - 3.0*h_list[id].GetRMS(),h_list[id].GetMean() + 3.0*h_list[id].GetRMS() )
        r = h_list[id].Fit(gaus,"QS","C", h_list[id].GetMean() - bin[id],h_list[id].GetMean() + bin[id] )

        gaus.SetParameter(0,r.Value(0))
        gaus.SetParameter(1,r.Value(1))
        gaus.SetParameter(2,r.Value(2))
        gaus.SetParameter(3,r.Value(3))
        gaus.SetParameter(4,r.Value(4))
        gaus.SetParameter(5,r.Value(5))
        gaus.Draw("same")

        f.write( str(mass[0]) + "," + str(r.Value(2)*r.Value(0)/(r.Value(0)+r.Value(3)) + r.Value(5)*r.Value(3)/(r.Value(0)+r.Value(3))) + "," + str(r.Value(2)*r.Value(0)/(r.Value(0)+r.Value(3))) + "," + str(r.Value(5)*r.Value(3)/(r.Value(0)+r.Value(3))) + "\n")

        c1.SaveAs("Plots/" + str(int(mass[0])) + "Zp.png")
        c1.Clear()
        # print "==============================",h_list[id].GetMean(),h_list[id].GetRMS(),h_list[id].GetEntries()


    tt.Write()
    f_out.Close()
    t.Delete()

# for i in range(len(h_list)):
#     print "==============================",Zp_mass[i],h_list[i].GetMean(),h_list[i].GetRMS(),h_list[i].GetEntries()
# for i in range(len(h_name)):
#     # h_name = "hist" + str(i)
#     h[i].GetXaxis().SetTitle(h_name[i])
#     h[i].GetYaxis().SetTitle("Events")
#     h[i].Draw("hist")
#     c1.SaveAs("Plots/" + h_name[i] + ".png")

ed = time.time()
print "time consumption:",ed - op,"s"