import ROOT
import os
import time
from math import sqrt,pow,fabs,log
from array import array
import ctypes

def lep_type(type,origin):
    lep = 0
    if type==6 and (origin==10 or origin==12 or origin==13 or origin==14 or origin==15 or origin==22 or origin==43):    lep = 1

    ################### photon conv
    ## iso photon from a prompt photon sau Higgs
    if type==4 and origin==5: lep = 2
    if type==14 and origin==37: lep = 2
    
    ############## bkg electrons from ElProc from a prompt photon
    if type==4 and origin==7: lep = 2 
    
    ### bkg photon from UndrPhot; don't have many when IN::chFlipElLep[index]!=0 (Here there is a generator level photon (not gen electron : that later converts:
    if type==16 and origin==38: lep = 2
    ### allow to select photons from parton shower, i.e. brem photons not isolated

    #####################//  hadron decay
    ### bkg electrons from DalitzDec, LightMeson and StrangeMeson
    if type==4 and (origin==6 or origin==23 or origin==24): lep = 3
    
    ### bkg electrons / bkg photons from pi0, LightMeson and StrangeMeson
    if type==16 and (origin==42 or origin==23): lep = 3
    
    ### bkg muons from PionDecay, KaonDecay, LightMeson and StrangeMeson
    if type==8 and (origin==34 or origin==35 or origin==23 or origin==24): lep = 3
    
    ### hadrons
    if type==17: lep = 3

    #####################/ heavy flavor
    ### iso / nonIso muons from BottomMeson, BBbarMeson and BottomBaryon
    if (type==6 or type==7) and (origin==26 or origin==29 or origin==33): lep = 4
    
    ### non iso electrons from CharmedMeson, CharmedBaryon and CCbarMeson
    if type==3 and (origin==25 or origin==32 or origin==27): lep = 4
    
    ### bkg electrons from CCbarMeson
    if type==4 and origin==27: lep = 4 ### JPsi
    
    ### bkg photons from CharmedMeson and CCbarMeson
    if type==16 and (origin==25 or origin==27): lep = 4
    
    ### non iso muons from CharmedMeson, CharmedBaryon and CCbarMeson
    if type==7 and (origin==25 or origin==32 or origin==27): lep = 4
    
    ### iso or bkg muons from CCbarMeson
    if (type==6 or type==8) and origin==27: lep = 4 ### JPsi

    return lep

op = time.time()

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat("")

c1 = ROOT.TCanvas("c1", "canvas", 1000, 800)

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.05)
c1.SetTopMargin(0.05)
c1.SetBottomMargin(0.15)

#   Reader initialization

colorList = [
    ROOT.kGreen,
    ROOT.kBlue,
    ROOT.kTeal,
    ROOT.kYellow,
    # ROOT.kGray,
    ROOT.kCyan,

    ROOT.kViolet,
    ROOT.kRed,
    ROOT.kOrange,
    ROOT.kGreen + 2,
    ROOT.kAzure,
    ROOT.kPink + 4,
    
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
    ROOT.kGreen + 2,
    ROOT.kAzure + 2,
    ROOT.kViolet + 2,
    ROOT.kOrange + 2,
    ROOT.kCyan + 2,
    ROOT.kYellow + 2,
    ROOT.kTeal + 2,
    ROOT.kMagenta + 2,
    ROOT.kRed + 2,
    ROOT.kBlue + 2
]

f = open("Plots/results.txt","a")
f.seek(0)
f.truncate()

f1 = open("Plots/cutflow.txt","a")
f1.seek(0)
f1.truncate()

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

WZ_list = [364253,364284,363358,364289]
ZZ_list = [345666,345723,364285,364250,364283,345283,245706,363356,364254]
WW_list = [345718]
VVV_list = [364242,364253,364244,364245,364246,364247,364248,364249]
Zjets_list = [361107,361666,361667,361106,361664,361665]
# Zjets_list = [364100,364101,364102,364103,364104,364105,364106,364107,364108,364109,364110,364111,364112,364113,364198,364199,364200,364201,364202,364203]
t_list = [410644,410645,410658,410659]
Wt_list = [410646,410647]
tt_list = [410472]
ttV_list = [410155,410156,410157,410081,410218,410219]
Zgamma_list = [364500,364501,364502,364503,364504,364505,364506,364507,364508,364509]
extra_list = [345705,345706,345714,364283,364284,364286,364287,364288,364289,364290,361108]
# extra_list = [364505,364506,364507,364508,364509]

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
            "VT", #10
            "HT", #11

            "mT_met_l3", #12
            "mT_sum", #13
            "M_lll", #14
            "M_ll_SF", #15
            "M_ll_Z1", #16
            "M_ll_Z2", #17

            "dR_1", #18
            "dR_2", #19
            "Delta_Phi", #20
            "METS", #21
            "dPhi_1", #22
            "dPhi_2", #23

            "LT", #24
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
            # "m_{T,met&l1} [GeV]", #8
            # "m_{T,met&l2} [GeV]", #9
            "m_{T,met&l} [GeV]", #8
            "m_{T,met&l_{W}} [GeV]", #9
            "HVT [GeV]", #10
            "HT [GeV]", #11

            "m_{T,met&l3} [GeV]", #12
            "m_{T,sum} [GeV]", #13
            "M_{lll} [GeV]", #14
            "M_{SF} [GeV]", #15
            "M_{ll,Z1} [GeV]", #16
            "M_{ll,Z2} [GeV]", #17

            "dR_{ll,1}", #18
            "dR_{ll,2}", #19
            # "#Delta#Phi_{3l&MET}", #12
            # "weight", #14
            # "DNN", #20
            # "Num. of loose muon" , #20
            "Lepton type" , #20
            "METS", #21
            "dPhi_{1}", #22
            "dPhi_{2}", #23

            "LT [GeV]", #24
            ]

xrange = [
    (0,120), #0
    (0,80), #1
    (0,40), #2
    (0,80), #3
    # (0,100), #4
    # (0,100), #5
    # (0,300), #0
    # (0,200), #1
    # (0,100), #2
    # (0,300), #3
    (70,110), #4
    (0,100), #5

    (0,7), #6
    (0,3), #7
    (0,100), #8
    (0,80), #9
    (0,70), #10
    (40,200), #11

    (20,200), #12
    (40,200), #13
    # (0,300), #14
    # (0,200), #9
    # (0,100), #10
    # (0,800), #11
    (0,200), #14
    (80,100), #15
    (0,100), #16
    (0,100), #17

    (0,7), #18
    (0,7), #19
    # (-1.5,1.5), #20
    (-0.5,4.5), #20
    (0,10), #21
    (-7,7), #22
    (-7,7), #23

    (0,40), #24
]

nbin = [
    20, #0
    20, #1
    20, #2
    20, #3
    20, #4
    20, #5

    8, #6
    4, #7
    20, #8
    20, #9
    20, #10
    20, #11

    20, #12
    20, #13
    20, #14
    20, #15
    20, #16
    20, #17
    20, #18
    20, #19
    # 20, #20
    5, #20
    20, #21
    20, #22
    20, #23
    20, #24
]

fake_plot = [
            "PT_tmu", #25
            "PT_lmu", #26
            "Eta_tmu", #27
            "Eta_lmu", #28
            ]

fake_name = [
            "tight lepton P_{T}",
            "loose lepton P_{T}",
            "tight lepton #eta",
            "loose lepton #eta",
            ]

fake_range = [
    (0,80), #0
    (0,80), #0
    (-3,3), #0
    (-3,3), #0
]

fake_bin = [40,40,10,10]

# ff = [
#     [-2.75,0.170747280121,0.0293622333556],
#     [-2.0,0.151015728712,0.00765188317746],
#     [-1.0,0.17029350996,0.00653260666877],
#     [0.0,0.171640500426,0.00534381670877],
#     [1.0,0.163910999894,0.0063368268311],
#     [2.0,0.1567068398,0.00771223707125],
#     [2.75,0.170660853386,0.0286666452885]
# ]

ff_pt = [
    [ 0.0 , 0.167953014374 , 0.00146010622848 ],
[ 8.0 , 0.125795558095 , 0.00176774535794 ],
[ 10.0 , 0.108285866678 , 0.0021853193175 ],
[ 12.0 , 0.0916015580297 , 0.00254960404709 ],
[ 14.0 , 0.0812354385853 , 0.00303624500521 ],
[ 16.0 , 0.0703721940517 , 0.00347858876921 ],
[ 18.0 , 0.060889814049 , 0.0041483216919 ],
[ 20.0 , 0.0627787560225 , 0.00355704524554 ],
]

ff_eta = [
    [ 0.0 , 0.2237855196 , 0.00483670597896 ],
[ 0.5 , 0.235718697309 , 0.00439087301493 ],
[ 1.5 , 0.272843420506 , 0.00587551202625 ],
[ 2.5 , 0.452562958002 , 0.0300703346729 ],
]

# ff = [
#     [0.250838558715,0.265682128787,0.314494933403],[0.22564471669,0.214173152287,0.267416463312],[0.200097396436,0.21626652649,0.232843860964],[0.177938208341,0.192302272372,0.252871435441],[0.17831897644,0.184221231386,0.213165431174],[0.188683192334,0.162743050902,0.221331416932],[0.141642772647,0.143572476656,0.192757404758],[0.186530734919,0.163169955057,0.250436221873]
# ]   # non-prompt included

ff = [
[0.251007277971,0.265905253859,0.314663106691],[0.225954862342,0.214486465002,0.26773882349],[0.200489319344,0.216681901492,0.233332704445],[0.178553327308,0.192790596829,0.253755718556],[0.179482054951,0.185102008035,0.214601349068],[0.189963666949,0.1644746145,0.22367989416],[0.143491939627,0.145485844722,0.195348543708],[0.193398094034,0.17089523883,0.260909763354]
]   # prompt

do_fake = 0
est_fake = 0
app_fake = 1

print "do_fake:",do_fake,"est_fake:",est_fake,"app_fake:",app_fake

if do_fake + est_fake + app_fake > 1:
    print "FFFFFFake!"
    os.exit()

if do_fake:
    plot_name = plot_name + fake_plot
    h_name = h_name + fake_name
    xrange = xrange + fake_range
    nbin = nbin + fake_bin
else:   fake_plot = []

if est_fake:
    f_fake = ROOT.TFile("fake/fake.root","recreate")
    f_fake.Close()
    
    f_t = ROOT.TTree("nominal","nominal")
    f_t.SetDirectory(0)
    
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
    
    lep0_d0sig = array('f',[0])
    lep1_d0sig = array('f',[0])
    lep2_d0sig = array('f',[0])

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

    f_w = array('f',[0.])
    f_wg = array('f',[0.])
    f_e = array('f',[0.])
    
    mass = array('f',[0.])
    bkg = array('i',[0])

    f_t.Branch('lep0_pt',lep0_pt,'lep0_pt/F')
    f_t.Branch('lep1_pt',lep1_pt,'lep1_pt/F')
    f_t.Branch('lep2_pt',lep2_pt,'lep2_pt/F')
    f_t.Branch('lep0_eta',lep0_eta,'lep0_eta/F')
    f_t.Branch('lep1_eta',lep1_eta,'lep1_eta/F')
    f_t.Branch('lep2_eta',lep2_eta,'lep2_eta/F')
    f_t.Branch('lep0_phi',lep0_phi,'lep0_phi/F')
    f_t.Branch('lep1_phi',lep1_phi,'lep1_phi/F')
    f_t.Branch('lep2_phi',lep2_phi,'lep2_phi/F')
    f_t.Branch('lep0_E',lep0_E,'lep0_E/F')
    f_t.Branch('lep1_E',lep1_E,'lep1_E/F')
    f_t.Branch('lep2_E',lep2_E,'lep2_E/F')
    f_t.Branch('lep0_id',lep0_id,'lep0_id/I')
    f_t.Branch('lep1_id',lep1_id,'lep1_id/I')
    f_t.Branch('lep2_id',lep2_id,'lep2_id/I')

    f_t.Branch('lep0_d0sig',lep0_d0sig,'lep0_d0sig/F')
    f_t.Branch('lep1_d0sig',lep1_d0sig,'lep1_d0sig/F')
    f_t.Branch('lep2_d0sig',lep2_d0sig,'lep2_d0sig/F')
    
    f_t.Branch('lep0_isoLoose',lep0_isoLoose,'lep0_isoLoose/I')
    f_t.Branch('lep1_isoLoose',lep1_isoLoose,'lep1_isoLoose/I')
    f_t.Branch('lep2_isoLoose',lep2_isoLoose,'lep2_isoLoose/I')
    
    f_t.Branch('lep0_isoTight',lep0_isoTight,'lep0_isoTight/I')
    f_t.Branch('lep1_isoTight',lep1_isoTight,'lep1_isoTight/I')
    f_t.Branch('lep2_isoTight',lep2_isoTight,'lep2_isoTight/I')

    f_t.Branch('met_met',met_met,'met_met/F')
    f_t.Branch('met_px',met_px,'met_px/F')
    f_t.Branch('met_py',met_py,'met_py/F')
    f_t.Branch('met_signif',met_signif,'met_signif/F')

    f_t.Branch('mll_1',mll_1,'mll_1/F')
    f_t.Branch('mll_2',mll_2,'mll_2/F')
    f_t.Branch('mll_Z1',mll_Z1,'mll_Z1/F')
    f_t.Branch('mll_Z2',mll_Z2,'mll_Z2/F')
    f_t.Branch('mlll',mlll,'mlll/F')

    f_t.Branch('HT',HT,'HT/F')
    f_t.Branch('LT',LT,'LT/F')
    f_t.Branch('VT',VT,'VT/F')

    f_t.Branch('mT',mT,'mT/F')
    f_t.Branch('mT_WZ',mT_WZ,'mT_WZ/F')
    f_t.Branch('mT_Wvl',mT_Wvl,'mT_Wvl/F')
    f_t.Branch('mT_vl',mT_vl,'mT_vl/F')

    f_t.Branch('dR_1',dR_1,'dR_1/F')
    f_t.Branch('dR_2',dR_2,'dR_2/F')

    f_t.Branch('dPhi_1',dPhi_1,'dPhi_1/F')
    f_t.Branch('dPhi_2',dPhi_2,'dPhi_2/F')

    f_t.Branch('n_jets',n_jets,'n_jets/I')
    f_t.Branch('n_bjets',n_bjets,'n_bjets/I')

    f_t.Branch('weight',f_w,'weight/F')
    f_t.Branch('weight_g',f_wg,'weight_g/F')
    f_t.Branch('err',f_e,'err/F')

    f_t.Branch('mass',mass,'mass/F')
    f_t.Branch('bkg',bkg,'bkg/I')

process = [
    "WZ",
    "ZZ",
    # "WW",
    "Z+jets",
    # "top(s)",
    # "t",
    # "Wt",
    "tt",
    # "ttV",
    # "VVV",
    # "Offshell H",
    "Zgamma",
    # "Others",
    # "Total background", # N
    # "Loose",
    # "LowPt",
    # "Medium",
    "5 GeV", # N
    "9 GeV",
    "15 GeV",
    # "19 GeV",
    # "23 GeV",
    # "27 GeV",
    # "31 GeV",
    # "35 GeV",
    "39 GeV",
    # "45 GeV",
    # "51 GeV",
    # "54 GeV",
    # "60 GeV",
    # "66 GeV",
    "69 GeV",
    "75 GeV",
    # "200 GeV",

    # "200 GeV-g=0.1",
    # "300 GeV-g=0.1",
    # "400 GeV-g=0.1",
    # "500 GeV-g=0.1",
    # "200 GeV-g=10e-4",

    # "95 GeV",
    # "150 GeV",
    # "250 GeV",
    # "350 GeV",
    # "500 GeV",
    
    "Data" # N + N_sig
]

N_sig = 6
N = len(process) - N_sig - 1

amp = 1.0
# lumi = 139/36.1
lumi = 139

# print N

hs = [ROOT.THStack() for i in range(len(plot_name))]
# ratio = [ROOT.TH1F() for i in range(len(plot_name))]
h_bkg = [ROOT.TH1F() for i in range(len(plot_name))]
h = [[ROOT.TH1F() for j in range(len(plot_name))] for i in range(len(process))]
ratio = [[ROOT.TH1F() for j in range(len(plot_name))] for i in range(N_sig)]

# if do_fake:
h_fake = [ROOT.TH1F() for i in range(len(plot_name))]

xbins = array("f",[0,8,10,12,14,16,18,20,80])
ybins = array("f",[0.,0.5,1.5,2.5])

h_f2d = [ROOT.TH2D() for i in range(4)]
h_f2d[0] = ROOT.TH2D("d_t2d","d_t2d",len(xbins)-1,xbins,len(ybins)-1,ybins) #data tight
h_f2d[1] = ROOT.TH2D("d_l2d","d_l2d",len(xbins)-1,xbins,len(ybins)-1,ybins) #data loose
h_f2d[2] = ROOT.TH2D("m_t2d","m_t2d",len(xbins)-1,xbins,len(ybins)-1,ybins) #MC tight
h_f2d[3] = ROOT.TH2D("m_l2d","m_l2d",len(xbins)-1,xbins,len(ybins)-1,ybins) #MC loose

h_2d = ROOT.TH2D("2d","2d",40,0,80,20,0,100)
# print h 

for i in range(len(plot_name)):
    # nbin[i] = 2*nbin[i]
    hs[i] = ROOT.THStack(plot_name[i],plot_name[i])
    h_bkg[i] = ROOT.TH1F(plot_name[i],plot_name[i], nbin[i], xrange[i][0], xrange[i][1])
    
    # if app_fake:
    h_fake[i] = ROOT.TH1F(plot_name[i],plot_name[i], nbin[i], xrange[i][0], xrange[i][1])
    h_fake[i].SetFillColor(ROOT.kRed)
    
    for j in range(N+N_sig+1):
        h[j][i] = ROOT.TH1F(process[j] + plot_name[i], process[j] + plot_name[i], nbin[i], xrange[i][0], xrange[i][1])
        if j < N:
            h[j][i].SetFillColor(int(colorList[j]))
            # h[j][i].SetLineColor(int(colorList[j]))
            # h[j][i].SetLineWidth(2)

        elif j != N+N_sig:
            ratio[j-N][i] = ROOT.TH1F(process[j] + plot_name[i], process[j] + plot_name[i], nbin[i], xrange[i][0], xrange[i][1])
            ratio[j-N][i].SetLineColor(int(colorList[j]))
            # ratio[j-N][i].SetMarkerColor(int(colorList[N+j]))
            ratio[j-N][i].SetLineWidth(2)
            # ratio[j-N][i].SetMarkerSize(2)

            h[j][i].SetLineColor(int(colorList[j]))
            h[j][i].SetLineWidth(2)
        # h[j][i].SetLineColor(ROOT.kRed)

        else:
            h[j][i].SetLineColor(ROOT.kBlack)
            h[j][i].SetLineWidth(2)

# len(cut) = 5
cut = [
    # "$3\mu$ events",
    # "Isolation tight",
    "$3l$ events",
    "$\mu\mu\mu$",
    "$p_T$ requirement",
    "$\\rm {N_{jets}}=0$",
    # "$M_{ll,1}<80$ GeV",
    "80 GeV $<M_{ll,1}<100$ GeV",
]

Events = [array('f',(N+N_sig+2)*[0.]) for i in range(len(cut))]
Err = [array('f',(N+N_sig+2)*[0.]) for i in range(len(cut))]
Raw = [array('i',(N+N_sig+2)*[0]) for i in range(len(cut))]

MVA_flag = 0
# input_dir = "./NoCut/"
input_dir = "./fake/"
if MVA_flag: input_dir = "MVA/no_mlll/high_39/"

# if not MVA_flag:   sample_list  = os.listdir(".")
# else:   sample_list = os.listdir("MVA/")
sample_list = os.listdir(input_dir)

B = 0.0
B_e = 0.0
# S = 0.0
# S_e = 0.0
D = 0.0

# print(sample_list)
for sample in sample_list:
    if "root" not in sample:    continue
        
    # if "364254" in sample or "364285" in sample or "364250" in sample:
    # if "364253_mc16a" not in sample:    continue
    # if "mc" in sample:    continue
    # if "mc16a" not in sample and "data15" not in sample and "data16" not in sample:    continue

    # if "data" in sample:    continue
    # if "_mc16" not in sample or int(sample[:6]) not in Zjets_list:   continue
    
    # if "Loose" not in sample and "LowPt" not in sample and "Med" not in sample:    continue

    # if "_mc16a" in sample or "_mc16d" in sample or "_mc16e" in sample:    continue
    if "muvZp" not in sample and "mc16" in sample :
        # print int(sample[:6])
        if int(sample[:6]) in WZ_list:    k = 0
        elif int(sample[:6]) in ZZ_list:    k = 1
        # elif int(sample[:6]) in WW_list:    k = 2
        elif int(sample[:6]) in Zjets_list:    k = 2
        # elif int(sample[:6]) in t_list:     k = 5
        # elif int(sample[:6]) in Wt_list:     k = 6
        elif int(sample[:6]) in tt_list:    k = 3
        # elif int(sample[:6]) in ttV_list:    k = 4
        elif int(sample[:6]) in Zgamma_list:    k = 4
        elif int(sample[:6]) in extra_list:    k = 5
        else:   continue

    elif "Zp005" in sample:    k = N
    elif "Zp009" in sample:     k = N + 1
    elif "Zp035" in sample:    k = N + 2
    elif "Zp039" in sample:    k = N + 3
    # elif "Zp045" in sample:    k = N + 3
    # elif "Zp051" in sample:    k = N + 4
    # elif "Zp054" in sample:     k = N + 4
    # elif "Zp060" in sample:     k = N + 5
    # elif "Zp066" in sample:     k = N + 4
    elif "Zp069" in sample:     k = N + 4
    elif "Zp075" in sample:    k = N + 5
    # elif "350.root" in sample:    k = N + N_sig - 1
    elif "data" in sample:    k = N + N_sig
    elif "fake.root" in sample and app_fake:  k = N + N_sig + 1
    elif MVA_flag:    k = 0
    else:    continue
    
    # if k != N + N_sig + 1:  continue
    # if k != 3:  continue
    if est_fake:
        if k > 1 and k != N + N_sig:  continue
    if app_fake:
        if k > 1 and k < N:  continue

    # if k != N + N_sig + 1:  k = N
    # else: k = N + 1

    sample = input_dir + sample

    # if "Loose" in sample:   k = N
    # elif "LowPt" in sample:   k = N + 1
    # elif "Med" in sample:   k = N + 2
    # else: continue
    # if MVA_flag:
    #     if "bkg" in sample: k = 0

    events = array('f',len(cut)*[0.])
    err = array('f',len(cut)*[0.])
    raw = array('i',len(cut)*[0])

    print k,sample

    # sample = "NoCut/" + sample
    f_in = ROOT.TFile(sample)
    t = f_in.Get("nominal")

    # for i in range(t.GetEntries()):
    for entry in t:
        # print "Processing:{}%".format(round((i + 1) * 100 / t.GetEntries(),2)), "\r"
        
        # if i > 50:
        #     break

        # t.GetEntry(i)

        if MVA_flag:
            if "bkg" in sample:
                if t.bkg < 4:    k = t.bkg
                else:    k = 4
            if "sig" in sample:
                # if t.mass == 19:    k = N
                # elif "Zp019" in sample:     k = N + 1
                if t.mass == 60:    k = N + 1
                # if t.mass == 15:    k = N + 2
                # if t.mass == 39:    k = N + 3
                # if t.mass == 69:    k = N + 4
                # elif t.mass == 75:    k = N + 5
                else:    continue

        if k == N+N_sig:
            # weight = t.weight*lumi
            weight = 1.0
        # elif k >= N:
        #     weight = t.weight_g*lumi
        else:
            weight = t.weight
            # weight = t.weight*lumi
            # weight = 0.
        # print weight
        
        # if t.weight < 0.7:    continue

        # if do_fake:
        # if t.lep0_id != -t.lep1_id or abs(t.lep2_id) != 13:
        # if abs(t.lep0_id) + abs(t.lep1_id) + abs(t.lep2_id) != 39:
        if abs(t.lep0_id) != 11 or abs(t.lep1_id) != 11 or abs(t.lep2_id) != 13:
        # if abs(t.lep0_id) != 13 or abs(t.lep1_id) != 13 or abs(t.lep2_id) != 11:
            continue

        if "fake.root" not in sample:
            fake = [0,0,0]
            # if abs(t.lep0_id) == 13:
            if t.lep0_d0sig > 3 or t.lep0_isoLoose == 0:    fake[0] = 1
            if t.lep1_d0sig > 3 or t.lep1_isoLoose == 0:    fake[1] = 1
                # if t.lep0_d0sig > 3:    fake[0] = 1
                # if t.lep1_d0sig > 3:    fake[1] = 1
            # else:
            #     if t.lep0_d0sig > 5 or t.lep0_isoLoose == 0:    fake[0] = 1
            #     if t.lep1_d0sig > 5 or t.lep1_isoLoose == 0:    fake[1] = 1
                # if t.lep0_d0sig > 5:    fake[0] = 1
                # if t.lep1_d0sig > 5:    fake[1] = 1
            
            if t.lep2_d0sig > 3 or t.lep2_isoLoose == 0:    fake[2] = 1
            
            if fake[0] or fake[1]:   continue
            
            if est_fake:
                if fake[2] != 1:    continue
            if app_fake:
                if fake[2]:    continue
                
            # if sum(fake) != 0:  continue

        events[0] += weight
        err[0] += weight*weight
        raw[0] += 1

        Events[0][k] += weight
        Err[0][k] += weight*weight
        Raw[0][k] += 1

        # if fabs(lep[0].Eta()) > 2.5 or fabs(lep[1].Eta()) > 2.5 or fabs(lep[2].Eta()) > 2.5 or t.lep2_pt < 3:
        #     continue

        # if abs(t.lep0_id) == 11 and abs(t.lep0_eta) > 1.37 and  abs(t.lep0_eta) < 2.47: continue
        # if abs(t.lep1_id) == 11 and abs(t.lep1_eta) > 1.37 and  abs(t.lep1_eta) < 2.47: continue
        # if abs(t.lep2_id) == 11 and abs(t.lep2_eta) > 1.37 and  abs(t.lep2_eta) < 2.47: continue

        # if t.lep2_isoTight == 0 or t.lep1_isoTight == 0 or t.lep0_isoTight == 0:    continue
        # if t.lep0_isoTight != 0 or t.lep1_isoTight != 0 or t.lep2_isoTight != 0:    continue
        # if t.lep2_isoTight == 0:    continue
        events[1] += weight
        err[1] += weight*weight
        raw[1] += 1

        Events[1][k] += weight
        Err[1][k] += weight*weight
        Raw[1][k] += 1

        if t.lep1_pt < 20:    continue
        if t.lep2_pt < 6:    continue
        # if t.lep1_pt < 25:    continue

        events[2] += weight
        err[2] += weight*weight
        raw[2] += 1

        Events[2][k] += weight
        Err[2][k] += weight*weight
        Raw[2][k] += 1

        if t.n_jets > 0:    continue
        # if t.n_bjets > 0:    continue
        events[3] += weight
        err[3] += weight*weight
        raw[3] += 1

        Events[3][k] += weight
        Err[3][k] += weight*weight
        Raw[3][k] += 1

        # if t.mlll > 80:    continue
        # if t.LT > 8:
            # continue
        if t.mll_1 > 80:    continue
        # if t.mll_1 < 80 or t.mll_1 > 100:     continue
        # if t.mll_Z1 < 80 or t.mll_Z1 > 100:     continue
        # if t.mT_Wvl > 25:    continue
        # if t.met_met > 25:    continue
        if MVA_flag:
            # if t.SRBDT < -0.2: continue
            if 'bkg' in sample: h_2d.Fill(t.mll_2,t.lep1_pt)
            if t.DNN_60 < 0.6: continue

        if abs(t.lep2_id) == 13 and "mc" in sample:
            lepType = lep_type(t.lep2_truthType,t.lep2_truthOrigin)
        else:
            lepType = -1
        
        # lepType = 1

        if est_fake:
            for i in range(len(xbins)-1):
                if lep[2].Pt() >= xbins[i] and lep[2].Pt() < xbins[i+1]:
                    for l in range(len(ybins)-1):
                        if fabs(lep[2].Eta()) >= ybins[l] and lep[2].Eta() < ybins[l+1]: f_w[0] = weight*ff[i][l]

            # if lep[2].Pt() >= ff_pt[-1][0]:
            #     f_w[0] = weight*ff_pt[-1][1]
            #     # f_e[0] = weight*sqrt(ff[-1][1]*ff[-1][1]+ff[-1][2]*ff[-1][2])
            # else:
            #     for i in range(len(ff_pt)-1):
            #         if lep[2].Pt() >= ff_pt[i][0] and lep[2].Pt() < ff_pt[i+1][0]:
            #             f_w[0] = weight*ff_pt[i][1]
            #             f_e[0] = weight*sqrt(ff[i][1]*ff[i][1]+ff[i][2]*ff[i][2])

            if k != N + N_sig:  
                if lepType:    f_w[0] = -f_w[0]
                else: continue

            lep0_pt[0] = t.lep0_pt
            lep1_pt[0] = t.lep1_pt
            lep2_pt[0] = t.lep2_pt
            lep0_eta[0] = t.lep0_eta
            lep1_eta[0] = t.lep1_eta
            lep2_eta[0] = t.lep2_eta
            lep0_phi[0] = t.lep0_phi
            lep1_phi[0] = t.lep1_phi
            lep2_phi[0] = t.lep2_phi
            lep0_E[0] = t.lep0_E
            lep1_E[0] = t.lep1_E
            lep2_E[0] = t.lep2_E
            lep0_id[0] = t.lep0_id
            lep1_id[0] = t.lep1_id
            lep2_id[0] = t.lep2_id

            lep0_d0sig[0] = t.lep0_d0sig
            lep1_d0sig[0] = t.lep1_d0sig
            lep2_d0sig[0] = t.lep2_d0sig

            # print(lep1_pt[0])
            met_met[0] = t.met_met
            met_px[0] = t.met_px
            met_py[0] = t.met_py
            met_signif[0] = t.met_signif

            lep0_isoLoose[0] = t.lep0_isoLoose
            lep1_isoLoose[0] = t.lep1_isoLoose
            lep2_isoLoose[0] = t.lep2_isoLoose
            
            lep0_isoTight[0] = t.lep0_isoTight
            lep1_isoTight[0] = t.lep1_isoTight
            lep2_isoTight[0] = t.lep2_isoTight

            n_jets[0] = t.n_jets
            n_bjets[0] = t.n_bjets

            mll_1[0] = t.mll_1
            mll_2[0] = t.mll_2
            mll_Z1[0] = t.mll_Z1
            mll_Z2[0] = t.mll_Z2
            mlll[0] = t.mlll
            
            HT[0] = t.HT 
            LT[0] = t.LT 
            VT[0] = t.VT

            mT[0] = t.mT 
            mT_WZ[0] = t.mT_WZ 
            mT_Wvl[0] = t.mT_Wvl 
            mT_vl[0] = t.mT_vl 
            
            dR_1[0] = t.dR_1 
            dR_2[0] = t.dR_2 
            dPhi_1[0] = t.dPhi_1 
            dPhi_2[0] = t.dPhi_2

            bkg[0] = -1
                    
            f_t.Fill()

        if app_fake:
            if k == 0 or k == 1:
                if not lepType: continue

        events[4] += weight
        err[4] += weight*weight
        raw[4] += 1

        Events[4][k] += weight
        Err[4][k] += weight*weight
        Raw[4][k] += 1
        
        # if t.weight < 0.7:
        #     print t.weight

        # mT_vl1 = pow(lep[0].Et() + t.met_met, 2) - pow( lep[0].Px() + t.met_px , 2) - pow ( lep[0].Py() + t.met_py , 2)
        # mT_vl2 = pow(lep[1].Et() + t.met_met, 2) - pow( lep[1].Px() + t.met_px , 2) - pow ( lep[1].Py() + t.met_py , 2)
        # mT_vl3 = pow(lep[2].Et() + t.met_met, 2) - pow( lep[2].Px() + t.met_px , 2) - pow ( lep[2].Py() + t.met_py , 2)

        # mT = pow(sqrt(lep[0].Pt()*lep[0].Pt() + 0.105*0.105) + sqrt(lep[1].Pt()*lep[1].Pt() + 0.105*0.105) + sqrt(lep[2].Pt()*lep[2].Pt() + 0.105*0.105) + t.met_met, 2) - pow( lep[0].Px() + lep[1].Px() + lep[2].Px() + t.met_px , 2) - pow ( lep[0].Py() + lep[1].Py() + lep[2].Py() + t.met_py , 2)

        # mT = pow(lep[0].Et() + lep[1].Et() + lep[2].Et() + t.met_met, 2) - pow( lep[0].Px() + lep[1].Px() + lep[2].Px() + t.met_px , 2) - pow ( lep[0].Py() + lep[1].Py() + lep[2].Py() + t.met_py , 2)

        # if mT_vl1 < 0:  mT_vl1 = 0.0
        # if mT_vl2 < 0:  mT_vl2 = 0.0
        # if mT_vl3 < 0:  mT_vl3 = 0.0
        # if mT < 0:
        #     mT = 0.0
        
        # mT_vl1 = sqrt(mT_vl1)
        # mT_vl2 = sqrt(mT_vl2)
        # mT_vl3 = sqrt(mT_vl3)
        # mT_vl1 = 0.0
        # mT_vl2 = 0.0
        # mT_vl3 = 0.0
        # mT = sqrt(mT)

        lep[0].SetPtEtaPhiE(t.lep0_pt,t.lep0_eta,t.lep0_phi,t.lep0_E)
        lep[1].SetPtEtaPhiE(t.lep1_pt,t.lep1_eta,t.lep1_phi,t.lep1_E)
        lep[2].SetPtEtaPhiE(t.lep2_pt,t.lep2_eta,t.lep2_phi,t.lep2_E)
        # delta_phi = ((lep[0]+lep[1]+lep[2]).Px()*t.met_px + (lep[0]+lep[1]+lep[2]).Py()*t.met_py)/sqrt((lep[0]+lep[1]+lep[2]).Pt()*(lep[0]+lep[1]+lep[2]).Pt()*t.met_met*t.met_met)
        # ss_m = (lep[0]+lep[1]).M()
        variables = [
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
            # mT_vl1, #8
            # mT_vl2, #9
            t.mT_vl, #8
            t.mT_Wvl, #9
            t.VT, #10
            t.HT, #11

            # mT_vl3, #10
            t.mT_WZ, #12
            # t.mT_Wvl, #9
            t.mT, #13
            t.mlll, #14
            # t.mll_1 - t.mll_2, #16
            (lep[1]+lep[0]).M(), #15
            t.mll_Z1, #16
            t.mll_Z2, #17

            t.dR_1, #18
            t.dR_2, #19
            # 1.0,
            # 1.0,
            # 1.0,
            lepType,
            # sum(fake), #20
            # t.weight, #14
            # t.DNN_39, #14
            t.met_signif, #21
            t.dPhi_1, #21
            t.dPhi_2, #22

            t.LT, #24
        ]
        
        # if k != 2:   h_2d.Fill(t.mll_1,t.lep1_pt)
        # h_2d.Fill(t.mlll,t.lep0_pt)
        # print weight[0],xsec[0],lumi
        # if k > N-1:
        #     print t.weight*lumi
        
        # if lep[2].Eta() > 2.4: print lep[2].Eta()
        for j in range(len(plot_name)-len(fake_plot)):
            if k < N + N_sig:
                if k < N:
                    h[k][j].Fill(variables[j],weight)
                    h_bkg[j].Fill(variables[j],weight)
                    if j == 0:
                        B += weight
                        B_e += weight*weight
                    # h[k][j].Fill(variables[j],1)
                else:
                    h[k][j].Fill(variables[j],weight*amp)
                # print Events[k]
            elif k == N + N_sig:
                h[N+N_sig][j].Fill(variables[j],1.)
                if j == 0:  D += 1           
            else:
                h_fake[j].Fill(variables[j],weight)
                h_bkg[j].Fill(variables[j],weight)
                if j == 0:
                    B += weight
                    B_e += weight*weight

        # for j in range(3):
        if not do_fake:   continue
        if k != N+N_sig:
            if k < N:
                h[k][len(plot_name)-4+fake[2]].Fill(lep[2].Pt(),weight)
                h_bkg[len(plot_name)-4+fake[2]].Fill(lep[2].Pt(),weight)
                
                h[k][len(plot_name)-2+fake[2]].Fill(lep[2].Eta(),weight)
                h_bkg[len(plot_name)-2+fake[2]].Fill(lep[2].Eta(),weight)

                if (k == 0 or k == 1) and lepType == 1: h_f2d[fake[2]+2].Fill(lep[2].Pt(),fabs(lep[2].Eta()),weight)
            else:
                h[k][len(plot_name)-4+fake[2]].Fill(lep[2].Pt(),weight*amp)
                h[k][len(plot_name)-2+fake[2]].Fill(lep[2].Eta(),weight*amp)
        else:
            h[N+N_sig][len(plot_name)-4+fake[2]].Fill(lep[2].Pt(),1.)
            h[N+N_sig][len(plot_name)-2+fake[2]].Fill(lep[2].Eta(),1.)

            h_f2d[fake[2]].Fill(lep[2].Pt(),fabs(lep[2].Eta()),1.)

        # if k < N and k != N+N_sig:
        #     B += weight
        #     B_e += weight*weight
    # for j in range(len(cut)):
    #     Events[j][k] += events[j]
    #     Err[j][k] += err[j]
    #     Raw[j][k] += raw[j]

    #     print(k,t.mass,Events[j][k])

    f.write("\n"+sample)
    f.write("\n"+"Total:"+str(events[0])+"+-"+str(sqrt(err[0])))
    f.write("\n"+"Isolation tight:"+str(events[1])+"+-"+str(sqrt(err[1])))
    f.write("\n"+"Lep pt cut:"+str(events[2])+"+-"+str(sqrt(err[2])))
    f.write("\n"+"Njets cut:"+str(events[3])+"+-"+str(sqrt(err[3])))
    f.write("\n"+"Mll_1 cut:"+str(events[4])+"+-"+str(sqrt(err[4]))+"\n")

print "bkg:",B,"+-",sqrt(B_e),"Data:",D,"+-",sqrt(D)

if est_fake:
    f_fake = ROOT.TFile.Open("fake/fake.root","update")
    f_t.Write()
    f_fake.Write()
    f_fake.Close()

h_2d.GetXaxis().SetTitle("m_ll_1 [GeV]")
h_2d.GetYaxis().SetTitle("lep1_pt [GeV]")
h_2d.GetXaxis().SetNdivisions(505, ROOT.kTRUE)
h_2d.SetMarkerSize(0.2)
h_2d.Draw("COLZ")

c1.SaveAs("cor.png")
c1.Clear()

# print "Signal:",Events[N-1],"Background:",Events[N],"Data:",Events[N+1]
# flow = ["Process","Total","Isolation tight","l3 pt cut","Mll_1 cut"]
flow_bkg = ["Process"]
flow_sig = ["Mass"]
flow_bkg = flow_bkg + cut
flow_sig = flow_sig + cut

for i in range(N+N_sig+1):
    f.write("\n"+process[i])
    if i < N or i == N + N_sig:   flow_bkg[0] = flow_bkg[0] + "&    " + process[i]
    else:   flow_sig[0] = flow_sig[0] + "&  " + process[i]

    for j in range(len(cut)):
        f.write("\n"+cut[j]+":"+str(Events[j][i])+"+-"+str(sqrt(Err[j][i]))+",Raw number:"+str(Raw[j][i]))
        if i < N or i == N + N_sig:   flow_bkg[j+1] = flow_bkg[j+1] + "&    " + str(round(Events[j][i],2)) + "$\\pm$" + str(round(sqrt(Err[j][i]),2))
        else:   flow_sig[j+1] = flow_sig[j+1] + "&  " + str(round(Events[j][i],2)) + "$\\pm$" + str(round(sqrt(Err[j][i]),2))

    f.write("\n"+"Final:"+str(Events[4][i])+"+-"+str(sqrt(Err[4][i]))+",Raw number:"+str(Raw[4][i])+"\n")

f.write("\nbkg:"+str(B)+"+-"+str(sqrt(B_e))+",Data:"+str(D)+"+-"+str(sqrt(D))+"\n")

######################  ratio calculation
error = ctypes.c_double()
rf_max = -9999
for j in range(len(plot_name)):
    rf_max = -9999
    if h_fake[j].Integral(): hs[j].Add(h_fake[j])
    for i in range(N+N_sig+1):
        # print process[i],Events[i]
        # f.write(process[i]+" "+str(Events[i])+"\n")
        if i < N:
            hs[j].Add(h[i][j])
        elif i != N+N_sig:
            signif = 0
            for k in range(h[i][j].GetSize()):
                B = h_bkg[j].GetBinContent(k)
                D = h[N+N_sig][j].GetBinContent(k)

                B_e = h_bkg[j].GetBinError(k)
                D_e = h[N+N_sig][j].GetBinError(k)

                if D <= 0 or B <= 0 :    ratio[0][j].SetBinContent(k,0)
                else:   
                    ratio[0][j].SetBinContent(k,D/B)
                    # ratio[0][j].SetBinError( k, sqrt( pow(B_e/D,2) + pow(B/D/D*D_e,2) ) )
                    ratio[0][j].SetBinError( k, sqrt( pow(D_e/B,2) + pow(D/B/B*B_e,2) ) )


                ####    Significance
                # B = h_bkg[j].Integral( k, h_bkg[j].GetSize() )
                # B = h_bkg[j].IntegralAndError( j, h_bkg[j].GetSize(), error, "" )
                # B_e = error.value
                # S = h[i][j].Integral( k, h_bkg[j].GetSize() )/amp
        
                # # if B <= 0 or S <= 0:    ratio[i-N][j].SetBinContent(k,0)
                # # else:       
                # try:
                #     ratio[i-N][j].SetBinContent(k, sqrt( 2*( (S+B)*log( (S+B)*(B+B_e*B_e)/( B*B+(S+B)*B_e*B_e ) ) - B*B/B_e/B_e*log( 1+B_e*B_e*S/B/(B+B_e*B_e) ) ) ) )
                #     # ratio[i-N][j].SetBinContent(k, S/sqrt(B) )
                # except:
                #     ratio[i-N][j].SetBinContent(k,0)

                # if h_bkg[j].GetBinContent(k) > 0:  signif = sqrt( signif*signif + pow( h[i][j].GetBinContent(k)/amp/sqrt(h_bkg[j].GetBinContent(k)), 2) )

                # ratio[i-N][j].SetBinContent(k, signif)
                
                # if h_bkg[j].Integral(0,k) <= 0:    ratio[i-N][j].SetBinContent(k,0)
                # else:   ratio[i-N][j].SetBinContent(k,h[i][j].Integral(0,k)/amp/sqrt(h_bkg[j].Integral(0,k)))

                # if h_bkg[j].Integral(k,h[i][j].GetSize()) <= 0:    ratio[i-N][j].SetBinContent(k,0)
                # else:   ratio[i-N][j].SetBinContent(k,h[i][j].Integral(k,h[i][j].GetSize())/amp/sqrt(h_bkg[j].Integral(k,h[i][j].GetSize())))
                # print(bc_1,bc_2,e_1,e_2)

            f.write(process[i]+":"+str(ratio[i-N][j].GetXaxis().GetBinCenter(ratio[i-N][j].GetMaximumBin()))+"\n")
            if rf_max < ratio[i-N][j].GetMaximum(): rf_max = ratio[i-N][j].GetMaximum()

    f.write(plot_name[j]+"\n")
    # print rf_max
    # ratio[0][j].GetYaxis().SetRangeUser(1e-3,rf_max*1.1)
    ratio[0][j].GetYaxis().SetRangeUser(0.,2.)
    ratio[0][j].GetYaxis().SetNdivisions(505, ROOT.kTRUE)
    ratio[0][j].GetYaxis().SetTitle ("Data/MC")
    # ratio[0][j].GetYaxis().SetTitle ("Significance")
    ratio[0][j].GetXaxis().SetLabelSize(0.1)
    ratio[0][j].GetYaxis().SetLabelSize(0.1)
    ratio[0][j].GetYaxis().SetTitleSize(0.1)
    ratio[0][j].GetYaxis().SetTitleOffset(0.3)

for i in range(len(flow_bkg)):
    f1.write(flow_bkg[i]+"\\\\\n\\midrule\n")
for i in range(len(flow_sig)):
    f1.write(flow_sig[i]+"\\\\\n\\midrule\n")

# c1.SetLogy()
# c1.SetLogx()

c2 = ROOT.TCanvas("c2", "canvas", 1500, 800)

c2.SetLeftMargin(0.0)
c2.SetRightMargin(0.0)
c2.SetTopMargin(0.05)
c2.SetBottomMargin(0.0)
c2.Divide(3,2,0.001,0.001)

leg_cache = []

gl = ROOT.TGraph()
hl = ROOT.TH1F()

########### plot
for i in range(len(plot_name)):
    # pad = ROOT.TPad(plot_name[i],plot_name[i],0 ,0.3 ,1 ,1)
    # frame = ROOT.gPad.DrawFrame(xrange[i][0],1e-3,xrange[i][1],scale)

    pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,1)
    pad2 = ROOT.TPad("pad2","pad2",0,0.05,1,0.3)
    pad1.SetBottomMargin(0.1)
    pad1.SetLeftMargin(0.15)
    pad2.SetTopMargin (0.)
    pad2.SetBottomMargin(0.25)
    pad2.SetLeftMargin(0.15)

    c1.cd()
    pad1.Draw()
    pad1.cd()
    pad1.SetLogy()
    # ROOT.gPad.Update()

    # h_ref = ROOT.TH1F("ref","ref",nbin[i],xrange[i][0],xrange[i][1])
    h_bkg[i].GetXaxis().SetTitle(h_name[i])
    h_bkg[i].GetYaxis().SetTitle("Events")
    if h[N+N_sig][i].GetMaximum() <= 0: h_bkg[i].GetYaxis().SetRangeUser(1e-3,hs[i].GetMaximum()*200.0/2)
    # else:   h_bkg[i].GetYaxis().SetRangeUser(1e-3,h[N][i].GetMaximum()*200.0/2)
    else:   h_bkg[i].GetYaxis().SetRangeUser(1e-3,h[N+N_sig][i].GetMaximum()*200.0/2)
    # print(hs[i].GetMaximum()*3.0/2)
    # h_bkg[i].GetYaxis().SetRangeUser(1e-3,h[N][i].GetMaximum()*3.0/2)
    # h_bkg[i].GetYaxis().SetRangeUser(1,hs[i].GetMaximum()*200)
    h_bkg[i].Draw("hist")

    pt = ROOT.TPaveText(0.33,0.85,0.23,0.7,"NDC")
    pt.SetBorderSize(0)
    pt.SetFillColor(0)
    # pt.SetFillStyle(0)
    pt.SetTextFont(42)
    pt.SetTextAlign(12)
    pt.SetTextSize(0.04)

    text = ROOT.TText()
    text = pt.AddText("#scale[1.3]{#it{ATLAS}   Internal}")
    text = pt.AddText("Low mass Z")
    text = pt.AddText("13 TeV, 139 fb^{-1}")
    # text=pt.AddText(categories[type])

    pt.Draw("same")
    leg_cache.append(pt)

    leg = ROOT.TLegend(0.88,0.88,0.75,0.88-len(process)*0.03)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.03)

    # leg.AddEntry(h[0][i],process[0],"f")

    hs[i].Draw("same hist")

    if h_fake[i].Integral(): leg.AddEntry(h_fake[j],"fake","f")
    for j in range(N+N_sig):
        if j < N:
            # continue
            leg.AddEntry(h[j][i],process[j],"f")
            # h[j][i].Draw("same hist")
        # leg.AddEntry(h[j+1][i],process[j+1],"f")
        else:
            leg.AddEntry(h[j][i],process[j],"f")
            # print(h_name[i],j,h_bkg[i].GetMaximum(),h[j][i].Integral())
            
            # if h[j][i].Integral() > 0: h[j][i].Scale(h_bkg[i].Integral()/h[j][i].Integral())
            h[j][i].Draw("same hist")
    
            c1.cd()
            pad2.Draw()
            pad2.cd()
            if j - N == 0:
                ratio[0][i].Draw("L")
                gl.SetPoint(0,xrange[i][0],1)
                gl.SetPoint(1,xrange[i][1],1)
                gl.SetLineStyle(ROOT.kDashed)
                gl.Draw("L same")
            else: ratio[j-N][i].Draw("L same")

            c1.cd()
            pad1.Draw()
            pad1.cd()

        # h[j+1][i].Draw("same hist")

    h[N+N_sig][i].Draw("same E1")
    leg.AddEntry(h[N+N_sig][i],"Data","lep")

    leg.Draw("same")

    # c1.cd()
    # pad2.Draw()
    # pad2.cd()
    # ratio[0][i].GetYaxis().SetRangeUser(0.,1.2)
    # ratio[0][i].Draw("E1")
    
    c1.SaveAs("Plots/" + plot_name[i] + ".png")
    c1.SaveAs("Plots/" + plot_name[i] + ".svg")
    # c1.SaveAs("Plots/" + plot_name[i] + ".pdf")
    c1.Clear()

    leg_cache.append(leg)
    # leg_cache.append(hl)
    # hl.Reset("ICESM")

    pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,1)
    pad2 = ROOT.TPad("pad2","pad2",0,0.05,1,0.3)
    pad1.SetBottomMargin(0.1)
    pad1.SetLeftMargin(0.15)
    pad2.SetTopMargin (0.)
    pad2.SetBottomMargin(0.25)
    pad2.SetLeftMargin(0.15)

    c2.cd(i+1-6*int(i/6))
    pad1.Draw()
    pad1.cd()
    pad1.SetLogy()
    # if i+1-6*int(i/6) == 1:  h_ref.Draw("hist")
    # else: h_ref.Draw("same hist")
    h_bkg[i].Draw("hist")
    pt.Draw("same")
    hs[i].Draw("same hist")
    for j in range(N+N_sig):
        if j >= N:
            # continue
        #     leg.AddEntry(h[j][i],process[j],"f")
        #     # h[j][i].Draw("same hist")
        # # leg.AddEntry(h[j+1][i],process[j+1],"f")
        # else:
        #     leg.AddEntry(h[j][i],process[j],"f")
            # h[j][i] = h[j][i]*amp
            h[j][i].Draw("same hist")
    
            c2.cd(i+1-6*int(i/6))
            pad2.Draw()
            pad2.cd()
            if j - N == 0: 
                ratio[0][i].Draw("L")
                gl.SetPoint(0,xrange[i][0],1)
                gl.SetPoint(1,xrange[i][1],1)
                gl.SetLineStyle(ROOT.kDashed)
                gl.Draw("L same")
            else: ratio[j-N][i].Draw("L same")

            c2.cd(i+1-6*int(i/6))
            pad1.Draw()
            pad1.cd()

    h[N+N_sig][i].Draw("same E1")
    leg.Draw("same")

    if (i+1)%6 == 0:
        c2.SaveAs("Plots/Multi/" + str((i+1)/6) +"_dist.png")
        # c2.SaveAs("Plots/Multi/" + str((i+1)/6) +"_dist.svg")
        # c2.SaveAs("Plots/Multi/" + str((i+1)/6) +"_dist.pdf")
        c2.Clear()
        c2.Divide(3,2,0.01,0.01)
        
c2.SaveAs("Plots/Multi/" + str((i+1)/6+1) +"_dist.png")
# c2.SaveAs("Plots/Multi/" + str((i+1)/6+1) +"_dist.svg")
# c2.SaveAs("Plots/Multi/" + str((i+1)/6+1) +"_dist.pdf")
c2.Clear()
# c2.Divide(3,2,0.01,0.01)

c1.cd()

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.SetTopMargin(0.05)
c1.SetBottomMargin(0.15)

if do_fake:
    for i in range(4):
        h[N+N_sig][len(plot_name)-4+i].Add(h[0][len(plot_name)-4+i],-1)
        h[N+N_sig][len(plot_name)-4+i].Add(h[1][len(plot_name)-4+i],-1)
        f.write(fake_name[i]+"\n")
        for k in range(h[2][len(plot_name)-4+i].GetSize()):
            f.write(process[2]+","+str(h[2][len(plot_name)-4+i].GetBinCenter(k))+","+str(h[2][len(plot_name)-4+i].GetBinContent(k))+","+str(h[2][len(plot_name)-4+i].GetBinError(k))+",   Data,"+str(h[N+N_sig][len(plot_name)-4+i].GetBinCenter(k))+","+str(h[N+N_sig][len(plot_name)-4+i].GetBinContent(k))+","+str(h[N+N_sig][len(plot_name)-4+i].GetBinError(k))+"\n")

    h_f2d[0].Add(h_f2d[2],-1)
    h_f2d[1].Add(h_f2d[3],-1)
    # xaxis = h_f2d[0].GetXaxis
    for i in range(1,len(xbins)):
        for j in range(1,len(ybins)):
            f.write(str(h_f2d[0].GetXaxis().GetBinCenter(i))+","+str(h_f2d[0].GetYaxis().GetBinCenter(j))+","+str(h_f2d[0].GetBinContent(i,j))+","+str(h_f2d[1].GetBinContent(i,j))+","+str(h_f2d[0].GetBinContent(i,j)*1.0/h_f2d[1].GetBinContent(i,j))+"\n")
            
            h_f2d[2].SetBinContent(i,j,h_f2d[0].GetBinContent(i,j)*1.0/h_f2d[1].GetBinContent(i,j))
    
    h_f2d[2].GetXaxis().SetTitle("p_{T} [GeV]")
    h_f2d[2].GetYaxis().SetTitle("|#eta|")
    h_f2d[2].GetXaxis().SetNdivisions(505, ROOT.kTRUE)
    h_f2d[2].SetMarkerSize(0.2)
    h_f2d[2].Draw("COLZ")

    c1.SaveAs("ff.png")
    c1.Clear()

f.close()

ed = time.time()
print "time consumption:",ed - op,"s"