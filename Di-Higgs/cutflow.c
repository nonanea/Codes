#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TPad.h>
#include <THStack.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

const TString samplePath;
const TString treeName;
const TString plotStorePath = "Plots/";
const TString ext = ".root";
const int N = 6;//Num. of processes
int n_type;
float VectorAngle(TLorentzVector Init, TLorentzVector Fin);
float PlaneAngle(TLorentzVector l_0,TLorentzVector l_1,TLorentzVector l_2,TLorentzVector l_3);

vector<int> colorList{kYellow, kViolet, kOrange, kGreen, kBlue, kRed, kTeal, kMagenta, kCyan, kAzure,
                    kBlue + 1, kRed + 1, kGreen + 1, kTeal + 1, kPink + 1, kMagenta + 1, kViolet + 1, kAzure + 1, kCyan + 1, kYellow + 1};

void cutflow(TString samplePath = "samples/", TString treeName = "nominal", int n_type = 0)
{
    gErrorIgnoreLevel = kError;

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(.03);

    //Object definition
    int triggers,RunYear,quadlep_type,total_charge;
    ULong64_t totalEvents;
    Int_t nJets_OR_T, nJets_OR_T_MV2c10_60, nJets_OR_T_MV2c10_85, nJets_OR_DL1r_77, nJets_OR_DL1r_85, nJets_OR, nJets_OR_MV2c10_85;
    Float_t lep_ID_0,lep_ID_1,lep_ID_2,lep_ID_3,lep_ID_4;
    Float_t lep_Pt_0,lep_Pt_1,lep_Pt_2,lep_Pt_3,lep_Pt_4;
    Float_t lep_Eta_0,lep_Eta_1,lep_Eta_2,lep_Eta_3,lep_Eta_4;
    Float_t lep_Phi_0,lep_Phi_1,lep_Phi_2,lep_Phi_3,lep_Phi_4;
    Float_t lep_E_0,lep_E_1,lep_E_2,lep_E_3,lep_E_4;
    Float_t lep_topoEtcone20_0,lep_topoEtcone20_1,lep_topoEtcone20_2,lep_topoEtcone20_3,lep_topoEtcone20_4;
    Float_t lep_topoEtcone30_0,lep_topoEtcone30_1,lep_topoEtcone30_2,lep_topoEtcone30_3,lep_topoEtcone30_4;
    Float_t lep_topoEtcone40_0,lep_topoEtcone40_1,lep_topoEtcone40_2,lep_topoEtcone40_3,lep_topoEtcone40_4;
    Float_t lep_ptVarcone20_0,lep_ptVarcone20_1,lep_ptVarcone20_2,lep_ptVarcone20_3,lep_ptVarcone20_4;
    Float_t lep_ptVarcone30_0,lep_ptVarcone30_1,lep_ptVarcone30_2,lep_ptVarcone30_3,lep_ptVarcone30_4;
    Float_t met_met, met_phi, HT, HT_lep, HT_jets;
    Float_t lead_jetPt, lead_jetEta, lead_jetPhi, lead_jetE;
    Float_t sublead_jetPt, sublead_jetEta, sublead_jetPhi, sublead_jetE;
    Char_t lep_isTrigMatch_0, lep_isTrigMatchDLT_0, lep_isLooseLH_0, lep_isLoose_0;
    Char_t lep_isTrigMatch_1, lep_isTrigMatchDLT_1, lep_isLooseLH_1, lep_isLoose_1;
    Char_t lep_isTrigMatch_2, lep_isTrigMatchDLT_2, lep_isLooseLH_2, lep_isLoose_2;
    Char_t lep_isTrigMatch_3, lep_isTrigMatchDLT_3, lep_isLooseLH_3, lep_isLoose_3;
    Char_t lep_isTrigMatch_4, lep_isTrigMatchDLT_4, lep_isLooseLH_4, lep_isLoose_4;
    Int_t lep_isolationLoose_0, lep_isolationLoose_1, lep_isolationLoose_2, lep_isolationLoose_3;
    Char_t lep_isolationFCLoose_0, lep_isolationFCLoose_1, lep_isolationFCLoose_2, lep_isolationFCLoose_3;
    Int_t lep_isolationPflowLoose_0, lep_isolationPflowLoose_1, lep_isolationPflowLoose_2, lep_isolationPflowLoose_3;
    Int_t lep_plvWP_Loose_0, lep_plvWP_Loose_1, lep_plvWP_Loose_2, lep_plvWP_Loose_3;
    Float_t lep_sigd0PV_0, lep_sigd0PV_1, lep_sigd0PV_2, lep_sigd0PV_3;
    Float_t lep_Z0SinTheta_0, lep_Z0SinTheta_1, lep_Z0SinTheta_2, lep_Z0SinTheta_3;
    Double_t scale_nom, mcWeightOrg, mc_xSection, mc_rawXSection, pileupEventWeight_090, JVT_EventWeight,mc_kFactor, totalWeights;
    Float_t lepSFObjLoose, lepSFObjTight, custTrigSF_LooseID_FCLooseIso_SLT, custTrigSF_LooseID_FCLooseIso_DLT, custTrigSF_LooseID_FCLooseIso_SLTorDLT;
    Float_t  weight_mc, weight_jvt, weight_pileup, weight_bTagSF_DL1r_77;
    float totalEventsWeighted, total_weights = 0;
    vector<float> *mcEventWeights = 0, *m_jet_pt = 0, *m_jet_eta = 0, *m_jet_phi = 0, *m_jet_E = 0;
                     
    int i = 0, j = 0, k = 0;   
    int p = 0, q = 0;
    int type = 0, mark = 0;
    float S = 0.0, B = 0.0, B_error = 0.0, data = 0, cs_jet = 0.0, cs_lep_12 = 0.0, cs_lep_34 = 0.0, cs_pairs = 0.0, cs_Z_pair= 0.0;
    double Z_m = 91187.6, tau_m = 1776.86, object = 0.0;

    // int events[N] = {0};
    Float_t ossf[N+1] = {0}, trig_match[N+1] = {0}, loose[N+1] = {0}, iso[N+1] = {0}, jet_cut[N+1] = {0}, lep_cut[N+1] = {0}, jet_lep[N+1] = {0}, lep_sepa[N+1] = {0}, lep_pt[N+1] = {0}, pair[N+1] = {0}, Jpsi[N+1] = {0}, jet_num[N+1] = {0}, b_tag[N+1] = {0}, m_llll[N+1] = {0}, nentries_weighted[N+1] = {0}, events_weighted[N+1] = {0}, events[N+1] = {0};
    Float_t ossf_err[N+1] = {0}, trig_match_err[N+1] = {0}, loose_err[N+1] = {0}, iso_err[N+1] = {0}, jet_cut_err[N+1] = {0}, lep_cut_err[N+1] = {0}, jet_lep_err[N+1] = {0}, lep_sepa_err[N+1] = {0}, lep_pt_err[N+1] = {0}, pair_err[N+1] = {0}, Jpsi_err[N+1] = {0}, jet_num_err[N+1] = {0}, b_tag_err[N+1] = {0}, m_llll_err[N+1] = {0}, start_error[N+1] = {0}, end_error[N+1] = {0};

    const char *categories[5] = {"","4mu","2mu2e","2e2mu","4e"};
    const char *plotname[99] = {"nJets",//0
                            "nbJets",//1
                            "m_jj",//2
                            "pt_jj",//3
                            "m_4l",//4
                            "pt_4l",//5
                            "m_2l_leading",//6
                            "m_2l_subleading",//7
                            "cos_jet",//8
                            "cos_leading",//9
                            "cos_subleading",//10
                            "met_met",//11
                            "cos_pairs_plane",//12
                            "cos_Z_pair_Plane",//13
                            "d0",//14
                            "z0sintheta",//15
                            "object"//16
                            };
    const char *xaxis_name[99] = {
                                    // "M_01 [GeV]",
                                    // "M_02 [GeV]",
                                    // "M_03 [GeV]",
                                    // "M_12 [GeV]",
                                    // "M_13 [GeV]",
                                    // "M_23 [GeV]",
                                    // "lep_Phi_0",
                                    // "lep_Phi_1",
                                    // "lep_Phi_2",
                                    // "lep_Phi_3",
                                    // "lep_Eta_0",
                                    // "lep_Eta_1",
                                    // "lep_Eta_2",
                                    // "lep_Eta_3",
                                    // "lep_ID_{leading}",//16
                                    // "lep_ID_{subleading}",//16
                                    "nJets",//0
                                    "nbJets",//1
                                    "M_{jj} [GeV]",//2
                                    "P_{T,jj} [GeV]",//3
                                    "M_{4l} [GeV]",//4
                                    "P_{T,4l} [GeV]",//5
                                    "M_{leading_pair} [GeV]",//6
                                    "M_{subleading_pair} [GeV]",//7
                                    "cos#theta_{0}",//8
                                    "cos#theta_{1}",//9
                                    "cos#theta_{2}",//10
                                    "MET [GeV]",//11
                                    "cos#theta_{Pairs Plane}",//12
                                    "cos#theta_{Z&Pair Plane}",//13
                                    "Sigd_{0}",//14
                                    "z_{0}sin#theta [mm]",//15
                                    "lep_isoFC_0 [GeV]",//16
                                    "lep_isoFC_1 [GeV]",//17
                                    "lep_isoFC_2 [GeV]",//18
                                    "lep_isoFC_3 [GeV]",//19
                                    "lep_isoPflow_0 [GeV]",//20
                                    "lep_isoPflow_1 [GeV]",//21
                                    "lep_isoPflow_2 [GeV]",//22
                                    "lep_isoPflow_3 [GeV]",//23
                                    "lep_ptVarcone_0/lep_Pt_0",//24
                                    "lep_ptVarcone_1/lep_Pt_1",//25
                                    "lep_ptVarcone_2/lep_Pt_2",//26
                                    "lep_ptVarcone_3/lep_Pt_3",//27
                                    "lep_topoEtcone_0/lep_Pt_0",//28
                                    "lep_topoEtcone_1/lep_Pt_1",//29
                                    "lep_topoEtcone_2/lep_Pt_2",//30
                                    "lep_topoEtcone_3/lep_Pt_3",//31
                                    "#Delta#Phi_{l_{0}&j_{0}}",//32
                                    "#Delta#Phi_{l_{1}&j_{0}}",//33
                                    "#Delta#Phi_{l_{2}&j_{0}}",//34
                                    "#Delta#Phi_{l_{3}&j_{0}}",//35
                                    "cos#theta_{l_{0}&j_{0}}",//36
                                    "cos#theta_{l_{1}&j_{0}}",//37
                                    "cos#theta_{l_{2}&j_{0}}",//38
                                    "cos#theta_{l_{3}&j_{0}}",//39
                                    "object"//40
                                    };
    int nbin[99]={
                    8,//0
                    3,//1
                    15,//2
                    15,//3
                    15,//4
                    15,//5
                    15,//6
                    15,//7
                    15,//8
                    15,//9
                    15,//10
                    15,//11
                    15,//12
                    15,//13
                    15,//14
                    15,//15
                    2,//16
                    2,//17
                    2,//18
                    2,//19
                    2,//20
                    2,//21
                    2,//22
                    2,//23
                    20,//24
                    20,//25
                    20,//26
                    20,//27
                    20,//28
                    20,//29
                    20,//30
                    20,//31
                    10,//32
                    10,//33
                    10,//34
                    10,//35
                    10,//36
                    10,//37
                    10,//38
                    10,//39
                    20//40
                    };
    float xrange[99][2] = {
                            // {0,110},
                            {0,10},//0
                            {0,4},//1
                            {0,600},//2
                            {0,400},//3
                            // {100,140},//4
                            {80,140},//4
                            {0,400},//5
                            {0,110},//6
                            {0,70},//7
                            {-1.2,1.2},//8
                            {-1.2,1.2},//9
                            {-1.2,1.2},//10
                            {0,300},//11
                            {-1.2,1.2},//12
                            {-1.2,1.2},//13
                            {-7,7},//14
                            {-0.5,0.5},//15
                            {0,2},//16
                            {0,2},//17
                            {0,2},//18
                            {0,2},//19
                            {0,2},//20
                            {0,2},//21
                            {0,2},//22
                            {0,2},//23
                            {-2,2},//24
                            {-2,2},//25
                            {-2,2},//26
                            {-2,2},//27
                            {-2,2},//28
                            {-2,2},//29
                            {-2,2},//30
                            {-4,4},//31
                            {-7,7},//32
                            {-7,7},//33
                            {-7,7},//34
                            {-7,7},//35
                            {-1.2,1.2},//36
                            {-1.2,1.2},//37
                            {-1.2,1.2},//38
                            {-1.2,1.2},//39
                            {-7,7}//40
                            };
    
    int obj = 41;

    THStack *hs[obj];
    for(q = 0; q < obj; q++) hs[q] = new THStack(plotname[q],plotname[q]);

    TH1D *h_err[obj];
    TH1D *hList[obj][N+1];
    // TH1D *hmet_weightedList[N], *hmet_phi_weightedList[N], *hjj_delta_phiList[N], *hjj_delta_RList[N], *hleading_delta_phiList[N], *hsubleading_delta_phiList[N], *hjet_productList[N], *hleading_productList[N];

    for(q = 0; q < N+1; q++)
    {
        for(p = 0; p < obj; p++)
        {
            hList[p][q] = new TH1D(plotname[p],plotname[p],nbin[p],xrange[p][0],xrange[p][1]);
            h_err[p] = new TH1D("error","error",nbin[p],xrange[p][0],xrange[p][1]);
        }
    }

    //Output file definition
    float lep_ID_0_f[N+2] = {0},lep_ID_1_f[N+2] = {0},lep_ID_2_f[N+2] = {0},lep_ID_3_f[N+2] = {0},lep_ID_4_f[N+2] = {0};
    float lep_Pt_0_f[N+2] = {0},lep_Pt_1_f[N+2] = {0},lep_Pt_2_f[N+2] = {0},lep_Pt_3_f[N+2] = {0},lep_Pt_4_f[N+2] = {0};
    float lep_Eta_0_f[N+2] = {0},lep_Eta_1_f[N+2] = {0},lep_Eta_2_f[N+2] = {0},lep_Eta_3_f[N+2] = {0},lep_Eta_4_f[N+2] = {0};
    float lep_Phi_0_f[N+2] = {0},lep_Phi_1_f[N+2] = {0},lep_Phi_2_f[N+2] = {0},lep_Phi_3_f[N+2] = {0},lep_Phi_4_f[N+2] = {0};
    float lep_iso_0_f[N+2] = {0},lep_iso_1_f[N+2] = {0},lep_iso_2_f[N+2] = {0},lep_iso_3_f[N+2] = {0},lep_iso_4_f[N+2] = {0};
    float jet_E_0_f[N+2] = {0},jet_E_1_f[N+2] = {0};
    float jet_Pt_0_f[N+2] = {0},jet_Pt_1_f[N+2] = {0};
    float jet_Eta_0_f[N+2] = {0},jet_Eta_1_f[N+2] = {0};
    float jet_Phi_0_f[N+2] = {0},jet_Phi_1_f[N+2] = {0};
    float HT_f[N+2] = {0}, HT_lep_f[N+2] = {0}, HT_jets_f[N+2] = {0}, met_met_f[N+2] = {0}, met_phi_f[N+2] = {0}, m_12_f[N+2] = {0}, m_34_f[N+2] = {0}, p_4l_f[N+2] = {0}, p_jj_f[N+2] = {0}, m_4l_f[N+2] = {0}, m_jj_f[N+2] = {0}, cs_jet_f[N+2] = {0}, cs_lep_12_f[N+2] = {0}, cs_lep_34_f[N+2] = {0}, cs_pairs_f[N+2] = {0}, cs_Z_pair_f[N+2] = {0}, Dphi_met_jets_f[N+2] = {0};
    float njets_f[N+2] = {0},nbjets_f[N+2] = {0};
    float mcWeight_f[N+2] = {0};

    //signal
    // float lep_ID_0_s,lep_ID_1_s,lep_ID_2_s,lep_ID_3_s,lep_ID_4_s;
    // float lep_Pt_0_s,lep_Pt_1_s,lep_Pt_2_s,lep_Pt_3_s,lep_Pt_4_s;
    // float lep_Eta_0_s,lep_Eta_1_s,lep_Eta_2_s,lep_Eta_3_s,lep_Eta_4_s;
    // float lep_Phi_0_s,lep_Phi_1_s,lep_Phi_2_s,lep_Phi_3_s,lep_Phi_4_s;
    // float jet_E_0_s,jet_E_1_s;
    // float jet_Pt_0_s,jet_Pt_1_s;
    // float jet_Eta_0_s,jet_Eta_1_s;
    // float jet_Phi_0_s,jet_Phi_1_s;
    // float HT_s, HT_lep_s, HT_jets_s, met_met_s, met_phi_s, m_12_s, m_34_s, p_4l_s, p_jj_s, m_4l_s, m_jj_s, cs_jet_s, cs_lep_12_s, cs_lep_34_s, cs_pairs_s, cs_Z_pair_s, Dphi_met_jets_s;
    // float njets_s,nbjets_s;
    // float mcWeight_s;

    // //background
    // float lep_ID_0_b,lep_ID_1_b,lep_ID_2_b,lep_ID_3_b,lep_ID_4_b;
    // float lep_Pt_0_b,lep_Pt_1_b,lep_Pt_2_b,lep_Pt_3_b,lep_Pt_4_b;
    // float lep_Eta_0_b,lep_Eta_1_b,lep_Eta_2_b,lep_Eta_3_b,lep_Eta_4_b;
    // float lep_Phi_0_b,lep_Phi_1_b,lep_Phi_2_b,lep_Phi_3_b,lep_Phi_4_b;
    // float jet_E_0_b,jet_E_1_b;
    // float jet_Pt_0_b,jet_Pt_1_b;
    // float jet_Eta_0_b,jet_Eta_1_b;
    // float jet_Phi_0_b,jet_Phi_1_b;
    // float HT_b, HT_lep_b, HT_jets_b, met_met_b, met_phi_b, m_12_b, m_34_b, p_4l_b, p_jj_b, m_4l_b, m_jj_b, cs_jet_b, cs_lep_12_b, cs_lep_34_b, cs_pairs_b, cs_Z_pair_b, Dphi_met_jets_b;
    // float njets_b,nbjets_b;
    // float mcWeight_b;
    int bkg;
    //float lep_E_0,lep_E_1,lep_E_2,lep_E_3,lep_E_4;

    TString rootFile;
    TString plotPath = samplePath;
    void *sampleDir = gSystem->OpenDirectory(samplePath);

    // Check Plot directory exist or not, if not, create it
    if (!gSystem->OpenDirectory(plotPath + plotStorePath + categories[n_type] + "_/"))
        gSystem->MakeDirectory(plotPath + plotStorePath + categories[n_type] + "_/");

    plotPath = plotPath + plotStorePath + categories[n_type] + "_/" + categories[type] + "_";

    cout<<plotPath<<endl;

    TFile *pre_f=new TFile( plotPath + "Preselection.root","recreate");
    
    TTree *tt[N+2];

    tt[0] = new TTree("tt","tt");
    tt[1] = new TTree("VV","VV");
    tt[2] = new TTree("ttV","ttV");
    tt[3] = new TTree("single_H","single_H");
    tt[4] = new TTree("Zjets","Zjets");
    // tt[5] = new TTree("VVV","VVV");
    tt[N-1] = new TTree("signal","signal");
    tt[N] = new TTree("data","data");
    tt[N+1] = new TTree("background","background");

    for(i = 0; i < N + 2; i++)
    {
        tt[i]->Branch("lep_ID_0",&lep_ID_0_f[i]);
        tt[i]->Branch("lep_Pt_0",&lep_Pt_0_f[i]);
        tt[i]->Branch("lep_Eta_0",&lep_Eta_0_f[i]);
        tt[i]->Branch("lep_Phi_0",&lep_Phi_0_f[i]);
        tt[i]->Branch("lep_iso_0",&lep_iso_0_f[i]);

        tt[i]->Branch("lep_ID_1",&lep_ID_1_f[i]);
        tt[i]->Branch("lep_Pt_1",&lep_Pt_1_f[i]);
        tt[i]->Branch("lep_Eta_1",&lep_Eta_1_f[i]);
        tt[i]->Branch("lep_Phi_1",&lep_Phi_1_f[i]);
        tt[i]->Branch("lep_iso_1",&lep_iso_1_f[i]);
        
        tt[i]->Branch("lep_ID_2",&lep_ID_2_f[i]);
        tt[i]->Branch("lep_Pt_2",&lep_Pt_2_f[i]);
        tt[i]->Branch("lep_Eta_2",&lep_Eta_2_f[i]);
        tt[i]->Branch("lep_Phi_2",&lep_Phi_2_f[i]);
        tt[i]->Branch("lep_iso_1",&lep_iso_1_f[i]);
        
        tt[i]->Branch("lep_ID_3",&lep_ID_3_f[i]); 
        tt[i]->Branch("lep_Pt_3",&lep_Pt_3_f[i]);
        tt[i]->Branch("lep_Eta_3",&lep_Eta_3_f[i]);
        tt[i]->Branch("lep_iso_3",&lep_iso_3_f[i]);
        
        tt[i]->Branch("jet_E_0",&jet_E_0_f[i]);
        tt[i]->Branch("jet_Pt_0",&jet_Pt_0_f[i]);
        tt[i]->Branch("jet_Eta_0",&jet_Eta_0_f[i]);
        tt[i]->Branch("jet_Phi_0",&jet_Phi_0_f[i]);

        tt[i]->Branch("jet_E_1",&jet_E_1_f[i]);
        tt[i]->Branch("jet_Pt_1",&jet_Pt_1_f[i]);
        tt[i]->Branch("jet_Eta_1",&jet_Eta_1_f[i]);
        tt[i]->Branch("jet_Phi_1",&jet_Phi_1_f[i]);
        
        tt[i]->Branch("m_12",&m_12_f[i]);
        tt[i]->Branch("m_34",&m_34_f[i]);
        tt[i]->Branch("m_4l",&m_4l_f[i]);
        tt[i]->Branch("m_jj",&m_jj_f[i]);
        tt[i]->Branch("p_4l",&p_4l_f[i]);
        tt[i]->Branch("p_jj",&p_jj_f[i]);
        tt[i]->Branch("cs_jet",&cs_jet_f[i]);
        tt[i]->Branch("cs_lep_12",&cs_lep_12_f[i]);
        tt[i]->Branch("cs_lep_34",&cs_lep_34_f[i]);
        tt[i]->Branch("cs_pairs",&cs_pairs_f[i]);
        tt[i]->Branch("cs_Z_pair",&cs_Z_pair_f[i]);
        tt[i]->Branch("met_met",&met_met_f[i]);
        tt[i]->Branch("HT",&HT_f[i]);
        tt[i]->Branch("HT_lep",&HT_lep_f[i]);
        tt[i]->Branch("Dphi_met_jets",&Dphi_met_jets_f[i]);

        tt[i]->Branch("njets",&njets_f[i]);
        tt[i]->Branch("nbjets",&nbjets_f[i]);

        tt[i]->Branch("mcWeight",&mcWeight_f[i]);
    }

    tt[N+1]->Branch("bkg",&bkg,"bkg/I");

    const char *entry;
    ofstream myfile, yieldfile;
    myfile.open (plotPath + "Results.txt");

    // const char *process_list[N+1] = {"tt","VV","ttV","Higgs","Z+jets","VVV","Signal","data"};
    const char *process_list[N+1] = {"tt","VV","ttV","Higgs","Z+jets","Signal","data"};

    while((entry = (char *)gSystem->GetDirEntry(sampleDir)))
    {
        rootFile = entry;
        //cout<<rootFile<<endl;

        if(rootFile.EndsWith(ext))
        {
            //cout<<p<<endl;
            cout<<samplePath+rootFile<<endl;
            myfile<<samplePath+rootFile<<endl;

            if(rootFile.Contains("_tt_"))   q = 0;
            // else if(rootFile.Contains("_ggZZ_"))   q = 1;
            else if(rootFile.Contains("_VV_"))   q = 1;
            else if(rootFile.Contains("_ttV_"))   q = 2;
            else if(rootFile.Contains("_Single_H_"))   q = 3;
            else if(rootFile.Contains("_Z+jets_"))   q = 4;
            // else if(rootFile.Contains("_VVV_"))   q = 5;
            else if(rootFile.Contains("_Signal_"))   q = N - 1;
            else if(rootFile.Contains("data"))   q = N;
            else
            {
                cout<<"Unknown file!"<<endl;
                continue;
            }

            // if(!(rootFile.Contains("34228"))) continue;
            // if( q != 0 && q != N && q!= 4)    continue;
            if(q != N-1)  continue;

            float cache[12] = {
                                trig_match[q],//0
                                // loose[q],
                                iso[q],//1
                                lep_sepa[q],//2
                                lep_pt[q],//3
                                ossf[q],//4
                                Jpsi[q],//5
                                jet_num[q],//6
                                b_tag[q],//7
                                m_llll[q],//8
                                nentries_weighted[q],//9
                                events_weighted[q],//10
                                events[q]//11
                                };

            float err_cache[12] = {
                                    trig_match_err[q],
                                    // loose_err[q],
                                    iso_err[q],
                                    lep_sepa_err[q],
                                    lep_pt_err[q],
                                    ossf_err[q],
                                    Jpsi_err[q],
                                    jet_num_err[q],
                                    b_tag_err[q],
                                    m_llll_err[q],
                                    start_error[q],
                                    end_error[q]
                                };

            int nentries, check = 0;
            double ab_eff = 0;
            int n_lep = 0;
            float lep_id[7] = {0};
            Float_t lumi = 0.0, mcWeight = 0.0, check_weighted = 0.0;
            TLorentzVector l,l_2,l_4,jet[2],lepton[7];

            TFile *f = new TFile(samplePath+rootFile,"read");
            TList *tl = (TList*) f->GetListOfKeys();

            // for(i = 0; i < l->GetEntries(); i++)
            // s.push_back((*l->At(i)).GetName());
            //TTree *t = (TTree*) f->Get("nominal");
            //TTree *t1 = (TTree*) f->Get("sumWeights");
            // TTree *t = (TTree*) f->Get("quadlep");
            TTree *t = (TTree*) f->Get(treeName);

            //t1->SetBranchAddress("totalEventsWeighted",&totalEventsWeighted);
            //t1->SetBranchAddress("totalEvents",&totalEvents);

            t->SetBranchAddress("RunYear",&RunYear);
            t->SetBranchAddress("quadlep_type",&quadlep_type);
            t->SetBranchAddress("total_charge",&total_charge);
            
            t->SetBranchAddress("lep_ID_0",&lep_ID_0);
            t->SetBranchAddress("lep_Pt_0",&lep_Pt_0);
            t->SetBranchAddress("lep_Eta_0",&lep_Eta_0);
            t->SetBranchAddress("lep_Phi_0",&lep_Phi_0);
            t->SetBranchAddress("lep_sigd0PV_0",&lep_sigd0PV_0);
            t->SetBranchAddress("lep_Z0SinTheta_0",&lep_Z0SinTheta_0);
            t->SetBranchAddress("lep_topoEtcone20_0",&lep_topoEtcone20_0);
            t->SetBranchAddress("lep_topoEtcone30_0",&lep_topoEtcone30_0);
            t->SetBranchAddress("lep_topoEtcone40_0",&lep_topoEtcone40_0);
            t->SetBranchAddress("lep_ptVarcone20_0",&lep_ptVarcone20_0);
            t->SetBranchAddress("lep_ptVarcone30_0",&lep_ptVarcone30_0);
            //t->SetBranchAddress("lep_E_0",&lep_E_0);

            t->SetBranchAddress("lep_ID_1",&lep_ID_1);
            t->SetBranchAddress("lep_Pt_1",&lep_Pt_1);
            t->SetBranchAddress("lep_Eta_1",&lep_Eta_1);
            t->SetBranchAddress("lep_Phi_1",&lep_Phi_1);
            t->SetBranchAddress("lep_sigd0PV_1",&lep_sigd0PV_1);
            t->SetBranchAddress("lep_Z0SinTheta_1",&lep_Z0SinTheta_1);
            t->SetBranchAddress("lep_topoEtcone20_1",&lep_topoEtcone20_1);
            t->SetBranchAddress("lep_topoEtcone30_1",&lep_topoEtcone30_1);
            t->SetBranchAddress("lep_topoEtcone40_1",&lep_topoEtcone40_1);
            t->SetBranchAddress("lep_ptVarcone20_1",&lep_ptVarcone20_1);
            t->SetBranchAddress("lep_ptVarcone30_1",&lep_ptVarcone30_1);
            //t->SetBranchAddress("lep_E_1",&lep_E_1);

            t->SetBranchAddress("lep_ID_2",&lep_ID_2);
            t->SetBranchAddress("lep_Pt_2",&lep_Pt_2);
            t->SetBranchAddress("lep_Eta_2",&lep_Eta_2);
            t->SetBranchAddress("lep_Phi_2",&lep_Phi_2);
            t->SetBranchAddress("lep_sigd0PV_2",&lep_sigd0PV_2);
            t->SetBranchAddress("lep_Z0SinTheta_2",&lep_Z0SinTheta_2);
            t->SetBranchAddress("lep_topoEtcone20_2",&lep_topoEtcone20_2);
            t->SetBranchAddress("lep_topoEtcone30_2",&lep_topoEtcone30_2);
            t->SetBranchAddress("lep_topoEtcone40_2",&lep_topoEtcone40_2);
            t->SetBranchAddress("lep_ptVarcone20_2",&lep_ptVarcone20_2);
            t->SetBranchAddress("lep_ptVarcone30_2",&lep_ptVarcone30_2);
            //t->SetBranchAddress("lep_E_2",&lep_E_2);

            t->SetBranchAddress("lep_ID_3",&lep_ID_3);
            t->SetBranchAddress("lep_Pt_3",&lep_Pt_3);
            t->SetBranchAddress("lep_Eta_3",&lep_Eta_3);
            t->SetBranchAddress("lep_Phi_3",&lep_Phi_3);
            t->SetBranchAddress("lep_sigd0PV_3",&lep_sigd0PV_3);
            t->SetBranchAddress("lep_Z0SinTheta_3",&lep_Z0SinTheta_3);
            t->SetBranchAddress("lep_topoEtcone20_3",&lep_topoEtcone20_3);
            t->SetBranchAddress("lep_topoEtcone30_3",&lep_topoEtcone30_3);
            t->SetBranchAddress("lep_topoEtcone40_3",&lep_topoEtcone40_3);
            t->SetBranchAddress("lep_ptVarcone20_3",&lep_ptVarcone20_3);
            t->SetBranchAddress("lep_ptVarcone30_3",&lep_ptVarcone30_3);
            //t->SetBranchAddress("lep_E_3",&lep_E_3);
            
            t->SetBranchAddress("lep_ID_4",&lep_ID_4);
            t->SetBranchAddress("lep_Pt_4",&lep_Pt_4);
            t->SetBranchAddress("lep_Eta_4",&lep_Eta_4);
            t->SetBranchAddress("lep_Phi_4",&lep_Phi_4);
            //t->SetBranchAddress("lep_E_4",&lep_E_4);

            t->SetBranchAddress("met_met",&met_met);
            t->SetBranchAddress("met_phi",&met_phi);
            t->SetBranchAddress("HT",&HT);
            t->SetBranchAddress("HT_lep",&HT_lep);
            t->SetBranchAddress("HT_jets",&HT_jets);

            /*t->SetBranchAddress("tau_pt_0",&tau_pt_0);
            t->SetBranchAddress("tau_eta_0",&tau_eta_0);
            t->SetBranchAddress("tau_phi_0",&tau_phi_0);
            t->SetBranchAddress("tau_E_0",&tau_E_0);
            t->SetBranchAddress("tau_charge_0",&tau_charge_0);

            t->SetBranchAddress("tau_pt_1",&tau_pt_1);
            t->SetBranchAddress("tau_eta_1",&tau_eta_1);
            t->SetBranchAddress("tau_phi_1",&tau_phi_1);
            t->SetBranchAddress("tau_E_1",&tau_E_1);
            t->SetBranchAddress("tau_charge_1",&tau_charge_1);
            
            //t->SetBranchAddress("m_truth_jet_pt",&m_truth_jet_pt);
            t->SetBranchAddress("m_jet_pt",&m_jet_pt);
            t->SetBranchAddress("m_jet_eta",&m_jet_eta);
            t->SetBranchAddress("m_jet_phi",&m_jet_phi);
            t->SetBranchAddress("m_jet_E",&m_jet_E);*/

            t->SetBranchAddress("lead_jetPt",&lead_jetPt);
            t->SetBranchAddress("lead_jetEta",&lead_jetEta);
            t->SetBranchAddress("lead_jetPhi",&lead_jetPhi);
            t->SetBranchAddress("lead_jetE",&lead_jetE);

            t->SetBranchAddress("sublead_jetPt",&sublead_jetPt);
            t->SetBranchAddress("sublead_jetEta",&sublead_jetEta);
            t->SetBranchAddress("sublead_jetPhi",&sublead_jetPhi);
            t->SetBranchAddress("sublead_jetE",&sublead_jetE);

            t->SetBranchAddress("lep_isTrigMatch_0",&lep_isTrigMatch_0);
            t->SetBranchAddress("lep_isTrigMatchDLT_0",&lep_isTrigMatchDLT_0);
            t->SetBranchAddress("lep_isLooseLH_0",&lep_isLooseLH_0);
            t->SetBranchAddress("lep_isLoose_0",&lep_isLoose_0);
            t->SetBranchAddress("lep_isolationLoose_0",&lep_isolationLoose_0);
            t->SetBranchAddress("lep_isolationFCLoose_0",&lep_isolationFCLoose_0);
            t->SetBranchAddress("lep_isolationPflowLoose_0",&lep_isolationPflowLoose_0);
            t->SetBranchAddress("lep_plvWP_Loose_0",&lep_plvWP_Loose_0);
            
            t->SetBranchAddress("lep_isTrigMatch_1",&lep_isTrigMatch_1);
            t->SetBranchAddress("lep_isTrigMatchDLT_1",&lep_isTrigMatchDLT_1);
            t->SetBranchAddress("lep_isLooseLH_1",&lep_isLooseLH_1);
            t->SetBranchAddress("lep_isLoose_1",&lep_isLoose_1);
            t->SetBranchAddress("lep_isolationLoose_1",&lep_isolationLoose_1);
            t->SetBranchAddress("lep_isolationFCLoose_1",&lep_isolationFCLoose_1);
            t->SetBranchAddress("lep_isolationPflowLoose_1",&lep_isolationPflowLoose_1);
            t->SetBranchAddress("lep_plvWP_Loose_1",&lep_plvWP_Loose_1);
            
            t->SetBranchAddress("lep_isTrigMatch_2",&lep_isTrigMatch_2);
            t->SetBranchAddress("lep_isTrigMatchDLT_2",&lep_isTrigMatchDLT_2);
            t->SetBranchAddress("lep_isLooseLH_2",&lep_isLooseLH_2);
            t->SetBranchAddress("lep_isLoose_2",&lep_isLoose_2);
            t->SetBranchAddress("lep_isolationLoose_2",&lep_isolationLoose_2);
            t->SetBranchAddress("lep_isolationFCLoose_2",&lep_isolationFCLoose_2);
            t->SetBranchAddress("lep_isolationPflowLoose_2",&lep_isolationPflowLoose_2);
            t->SetBranchAddress("lep_plvWP_Loose_2",&lep_plvWP_Loose_2);
            
            t->SetBranchAddress("lep_isTrigMatch_3",&lep_isTrigMatch_3);
            t->SetBranchAddress("lep_isTrigMatchDLT_3",&lep_isTrigMatchDLT_3);
            t->SetBranchAddress("lep_isLooseLH_3",&lep_isLooseLH_3);
            t->SetBranchAddress("lep_isLoose_3",&lep_isLoose_3);
            t->SetBranchAddress("lep_isolationLoose_3",&lep_isolationLoose_3);
            t->SetBranchAddress("lep_isolationFCLoose_3",&lep_isolationFCLoose_3);
            t->SetBranchAddress("lep_isolationPflowLoose_3",&lep_isolationPflowLoose_3);
            t->SetBranchAddress("lep_plvWP_Loose_3",&lep_plvWP_Loose_3);
            
            t->SetBranchAddress("lep_isTrigMatch_4",&lep_isTrigMatch_4);
            t->SetBranchAddress("lep_isTrigMatchDLT_4",&lep_isTrigMatchDLT_4);
            t->SetBranchAddress("lep_isLooseLH_4",&lep_isLooseLH_4);
            t->SetBranchAddress("lep_isLoose_4",&lep_isLoose_4);
            //t->SetBranchAddress("lep_isolationLoose_4",&lep_isolationLoose_4);

            //t->SetBranchAddress("nJets_OR_T",&nJets_OR_T);
            //t->SetBranchAddress("nJets_OR_T_MV2c10_60",&nJets_OR_T_MV2c10_60);
            //t->SetBranchAddress("nJets_OR_T_MV2c10_85",&nJets_OR_T_MV2c10_85);
            t->SetBranchAddress("nJets_OR",&nJets_OR);
            //t->SetBranchAddress("nJets_OR_MV2c10_85",&nJets_OR_MV2c10_85);
            t->SetBranchAddress("nJets_OR_DL1r_77",&nJets_OR_DL1r_77);
            t->SetBranchAddress("nJets_OR_DL1r_85",&nJets_OR_DL1r_85);

            t->SetBranchAddress("scale_nom",&scale_nom);
            t->SetBranchAddress("totalWeights",&totalWeights);
            //t->SetBranchAddress("pileupEventWeight_090",&pileupEventWeight_090);
            //t->SetBranchAddress("JVT_EventWeight",&JVT_EventWeight);
            //t->SetBranchAddress("mcWeightOrg",&mcWeightOrg);
            t->SetBranchAddress("weight_pileup",&weight_pileup);
            t->SetBranchAddress("weight_jvt",&weight_jvt);
            t->SetBranchAddress("weight_mc",&weight_mc);
            t->SetBranchAddress("weight_bTagSF_DL1r_77",&weight_bTagSF_DL1r_77);
            //t->SetBranchAddress("mcEventWeights",&mcEventWeights);
            t->SetBranchAddress("mc_xSection",&mc_xSection);
            t->SetBranchAddress("mc_rawXSection",&mc_rawXSection);
            t->SetBranchAddress("mc_kFactor",&mc_kFactor);
            t->SetBranchAddress("lepSFObjLoose",&lepSFObjLoose);
            t->SetBranchAddress("lepSFObjTight",&lepSFObjTight);
            t->SetBranchAddress("custTrigSF_LooseID_FCLooseIso_SLT",&custTrigSF_LooseID_FCLooseIso_SLT);
            t->SetBranchAddress("custTrigSF_LooseID_FCLooseIso_DLT",&custTrigSF_LooseID_FCLooseIso_DLT);
            t->SetBranchAddress("custTrigSF_LooseID_FCLooseIso_SLTorDLT",&custTrigSF_LooseID_FCLooseIso_SLTorDLT);

            nentries = t->GetEntries();

            for(i = 0; i < nentries; i++)
            {
                vector<TLorentzVector> lep;
                vector<float> lep_ID;
                vector<int> jet_id;
                int ID[4] = {0};
                
                t->GetEntry(i);
    
                double delta_m = 0, m_12 = 0, m_34 = 0, delta_R = 0, sum_weights = 1;

                //if(quadlep_type < 1)   continue;
                //if(total_charge != 0)   continue;
                //if(i > 100)    break;
                //MCweight: "(36074.6*(RunYear==2015||RunYear==2016)+43813.7*(RunYear==2017)+ 58450.1*(RunYear==2018))*pileupEventWeight_090*JVT_EventWeight*mc_xSection*mc_kFactor*mcWeightOrg/totalweightedevents"
                //cout<<RunYear<<endl;

                if(RunYear == 2015 || RunYear == 2016)  lumi = 36074.6;
                else if(RunYear == 2017) lumi = 43813.7;
                else if(RunYear == 2018) lumi = 58450.1;
                else
                {
                    cout<<"No lumi found!"<<endl;
                    break;
                }

                // weight_jvt = 1.0;
                // weight_bTagSF_DL1r_77 = 1.0;
                //if(rootFile.Contains("450578")) mcWeight = lumi*pileupEventWeight_090*JVT_EventWeight*mcWeightOrg*4.8e-6*scale_nom/mc_xSection;
                // if(rootFile.Contains("450578")) mcWeight = lumi*weight_pileup*weight_jvt*weight_mc*9.68e-6*mc_kFactor/totalWeights;
                if(q == 0 || q == 1 || q == 2 || q == 4)    mcWeight = lumi*weight_pileup*weight_jvt*weight_mc*weight_bTagSF_DL1r_77*lepSFObjLoose*custTrigSF_LooseID_FCLooseIso_SLTorDLT*mc_xSection/totalWeights;
                
                else if(q == N) mcWeight = 1.0;
                
                else if(rootFile.Contains("342284")) mcWeight = lumi*weight_pileup*weight_jvt*weight_mc*weight_bTagSF_DL1r_77*lepSFObjLoose*custTrigSF_LooseID_FCLooseIso_SLTorDLT*1.37*mc_kFactor/totalWeights;
                
                else if(rootFile.Contains("342285")) mcWeight = lumi*weight_pileup*weight_jvt*weight_mc*weight_bTagSF_DL1r_77*lepSFObjLoose*custTrigSF_LooseID_FCLooseIso_SLTorDLT*0.88*mc_kFactor/totalWeights;

                else if(rootFile.Contains("346345")) mcWeight = lumi*weight_pileup*weight_jvt*weight_mc*weight_bTagSF_DL1r_77*lepSFObjLoose*custTrigSF_LooseID_FCLooseIso_SLTorDLT*0.056448*mc_kFactor/totalWeights;

                else if(rootFile.Contains("341471")) mcWeight = lumi*weight_pileup*weight_jvt*weight_mc*weight_bTagSF_DL1r_77*lepSFObjLoose*custTrigSF_LooseID_FCLooseIso_SLTorDLT*0.012989*mc_kFactor/totalWeights;
                
                else if(rootFile.Contains("341488")) mcWeight = lumi*weight_pileup*weight_jvt*weight_mc*weight_bTagSF_DL1r_77*lepSFObjLoose*custTrigSF_LooseID_FCLooseIso_SLTorDLT*0.001010*mc_kFactor/totalWeights;

                else mcWeight = lumi*weight_pileup*weight_jvt*weight_mc*weight_bTagSF_DL1r_77*lepSFObjLoose*custTrigSF_LooseID_FCLooseIso_SLTorDLT*9.68e-6/totalWeights; // Signal

                // else if(rootFile.Contains("410472"))    mcWeight = 729.77/87.7076*lumi*weight_pileup*weight_jvt*weight_mc*weight_bTagSF_DL1r_77*lepSFObjLoose*custTrigSF_LooseID_FCLooseIso_SLTorDLT*mc_xSection*mc_kFactor/totalWeights;


                //mcWeight = lumi*pileupEventWeight_090*JVT_EventWeight*mcWeightOrg*scale_nom;
                //else    mcWeight = lumi*weight_pileup*weight_jvt*weight_mc*scale_nom;
                
                // else    mcWeight = lumi*weight_pileup*weight_jvt*weight_mc*mc_xSection*mc_kFactor/totalWeights;
                
                //cout<<weight_mc*weight_pileup*weight_jvt*scale_nom/mc_xSection<<endl;
                //cout<<mc_kFactor<<endl;
                nentries_weighted[q] = nentries_weighted[q] + mcWeight;
                start_error[q] = start_error[q] +mcWeight*mcWeight;

                jet[0].SetPtEtaPhiE(lead_jetPt,lead_jetEta,lead_jetPhi,lead_jetE);
                jet[1].SetPtEtaPhiE(sublead_jetPt,sublead_jetEta,sublead_jetPhi,sublead_jetE);

                lepton[0].SetPtEtaPhiM(lep_Pt_0,lep_Eta_0,lep_Phi_0,0);
                lepton[1].SetPtEtaPhiM(lep_Pt_1,lep_Eta_1,lep_Phi_1,0); 
                lepton[2].SetPtEtaPhiM(lep_Pt_2,lep_Eta_2,lep_Phi_2,0); 
                lepton[3].SetPtEtaPhiM(lep_Pt_3,lep_Eta_3,lep_Phi_3,0);
                //lepton[4].SetPtEtaPhiM(lep_Pt_4,lep_Eta_4,lep_Phi_4,0);

                lep_id[0] = 1.0*lep_ID_0;
                lep_id[1] = 1.0*lep_ID_1;
                lep_id[2] = 1.0*lep_ID_2;
                lep_id[3] = 1.0*lep_ID_3;
                //lep_id[4] = 1.0*lep_ID_4;

                delta_m = 999999;
                m_12 = 0.0;
                m_34 = 0.0;

                for(j = 0; j < 4; j++)
                {
                    lep.push_back(lepton[j]);
                    lep_ID.push_back(lep_id[j]);
                }
                n_lep = lep.size();
                //if(n_lep < 4)   continue;

                //Cutflow start
                j = 0;
                // if(lep_isTrigMatchDLT_0)    j++;
                // if(lep_isTrigMatchDLT_1)    j++;
                // if(lep_isTrigMatchDLT_2)    j++;
                // if(lep_isTrigMatchDLT_3)    j++;
                if(!(lep_isTrigMatch_0||lep_isTrigMatch_1||lep_isTrigMatch_2||lep_isTrigMatch_3||lep_isTrigMatchDLT_0||lep_isTrigMatchDLT_1||lep_isTrigMatchDLT_2||lep_isTrigMatchDLT_3))   continue;
                // if(!(lep_isTrigMatch_0&&lep_isTrigMatch_1&&lep_isTrigMatch_2&&lep_isTrigMatch_3) && j < 2)   continue;
                trig_match[q]=trig_match[q]+mcWeight;
                trig_match_err[q]=trig_match_err[q]+mcWeight*mcWeight;

                // if(abs(lep_ID[ID[0]]) == 13 && abs(lep_ID[ID[2]]) == 13)    type = 1;
                // else if(abs(lep_ID[ID[0]]) == 13 && abs(lep_ID[ID[2]]) == 11)   type = 2;
                // else if(abs(lep_ID[ID[0]]) == 11 && abs(lep_ID[ID[2]]) == 13)   type = 3;
                // else if(abs(lep_ID[ID[0]]) == 11 && abs(lep_ID[ID[2]]) == 11)   type = 4;
                // else    type = 0;

                // if(n_type == 0) 
                type = 0;

                // if(type != n_type)
                // {
                //     type = n_type;
                //     continue;
                // }

                //if(lep_isolationFCLoose_0)    cout<<lep_isolationFCLoose_0<<endl;
                //if(lep_isolationPflowLoose_0 == 0)    cout<<lep_isolationPflowLoose_0<<endl;                
                //if((abs(lep_ID_0) == 11 && !lep_isolationFCLoose_0) 
                //    || (abs(lep_ID_1) == 11 && !lep_isolationFCLoose_1) 
                //    || (abs(lep_ID_2) == 11 && !lep_isolationFCLoose_2) 
                //    || (abs(lep_ID_3) == 11 && !lep_isolationFCLoose_3)) continue;

                //if((abs(lep_ID_0) == 13 && lep_isolationPflowLoose_0 == 0)
                //    || (abs(lep_ID_1) == 13 && lep_isolationPflowLoose_1 == 0)
                //    || (abs(lep_ID_2) == 13 && lep_isolationPflowLoose_2 == 0)
                //    || (abs(lep_ID_3) == 13 && lep_isolationPflowLoose_3 == 0)) continue;

                j = 0;
                // if(lep_isolationPflowLoose_0)   j++;
                // if(lep_isolationPflowLoose_1)   j++;
                // if(lep_isolationPflowLoose_2)   j++;
                // if(lep_isolationPflowLoose_3)   j++;
                // if(j < 4)   continue; 

                // if(lep_isolationPflowLoose_3 == 0) continue;
                // if(lep_isolationPflowLoose_1 == 0 || lep_isolationPflowLoose_2 == 0 || (lep_isolationPflowLoose_2 == 0 && lep_isolationPflowLoose_3 == 0)) continue;

                // if(!lep_isolationFCLoose_2 && !lep_isolationFCLoose_3) continue;
                // if(lep_plvWP_Loose_2 == 0 && lep_plvWP_Loose_3 == 0) continue;
                if(lep_isolationPflowLoose_2 == 0 && lep_isolationPflowLoose_3 == 0) continue;
                
                iso[q]=iso[q]+mcWeight;
                iso_err[q]=iso_err[q]+mcWeight*mcWeight;
                
                //if(fabs(lep_sigd0PV_0) > 5) continue;
                //if(fabs(lep_sigd0PV_1) > 5) continue;
                //if(fabs(lep_sigd0PV_2) > 5) continue;
                //if(fabs(lep_sigd0PV_3) > 5) continue;
                //if(fabs(lep_Z0SinTheta_0) > 0.5) continue;
                //if(fabs(lep_Z0SinTheta_1) > 0.5) continue;
                //if(fabs(lep_Z0SinTheta_2) > 0.5) continue;
                //if(fabs(lep_Z0SinTheta_3) > 0.5) continue;

                /*for(j = 0; j < 2; j++)
                {
                    k = 0;
                    if((jet[j].Pt() < 20000 || fabs(jet[j].Eta()) > 2.5) && (jet[j].Pt() < 30000 || fabs(jet[j].Eta()) < 2.5 || fabs(jet[j].Eta()) > 4.5))
                    {
                        k = 1;
                        break;
                    }
                }
                //if(k == 1)  continue;
                //cout<<"outside:"<<i<<","<<j<<",jetpt:"<<jet[j-1].Pt()<<endl;
                jet_cut[q]=jet_cut[q]+mcWeight;
                jet_cut_err[q]=jet_cut_err[q]+mcWeight*mcWeight;*/
                
                delta_R = 10.0;

                check++;
                check_weighted = check_weighted + mcWeight;
    
                //lepton separation
                for(j = 1; j < n_lep; j++)
                {
                    for(k = 0; k < j; k++)
                    {
                        if(delta_R < 0.1) break;

                        delta_R = sqrt((lep[j].Eta()-lep[k].Eta())*(lep[j].Eta()-lep[k].Eta())+(lep[j].Phi()-lep[k].Phi())*(lep[j].Phi()-lep[k].Phi()));
                    }
                }
                if(delta_R < 0.1) continue;
                lep_sepa[q]=lep_sepa[q]+mcWeight;
                lep_sepa_err[q]=lep_sepa_err[q]+mcWeight*mcWeight;

                k = 0;
                //letpon kinematic
                for(j = 0; j < n_lep; j++)
                {
                    if(lep[j].Pt() > 20000) k++;
                }
                if(k < 1)  continue;
                else k = 0;

                for(j = 0; j < n_lep; j++)
                {
                    if(lep[j].Pt() > 15000) k++;
                }
                if(k < 2)  continue;
                else k = 0;

                for(j = 0; j < n_lep; j++)
                {
                    if(lep[j].Pt() > 10000) k++;
                }
                if(k < 3)  continue;
                else k = 0;

                // for(j = 0; j < n_lep; j++)
                // {
                //     if(lep[j].Pt() > 5000) k++;
                // }
                // if(k < 4)  continue;
                // else k = 0;

                lep_pt[q]=lep_pt[q]+mcWeight;
                lep_pt_err[q]=lep_pt_err[q]+mcWeight*mcWeight;

                //lepton pair selection
                for(j = 1; j < n_lep; j++)
                {
                    for(k = 0; k < j; k++)
                    {
                        if(lep_ID[j] == -lep_ID[k] && lep_ID[j] != 0 )
                        {
                            l_2 = lep[j]+lep[k];
                            if(delta_m > fabs(Z_m - l_2.M()))
                            {
                                m_12 = l_2.M();
                                ID[0] = j;
                                ID[1] = k;
                                
                                delta_m = fabs(Z_m - l_2.M());
                            }
                        }
                    }
                }
                // if(lep_ID[ID[0]] != -lep_ID[ID[1]]) continue;

                // for(j = 1; j < n_lep; j++)
                // {
                //     if(j == ID[0] || j == ID[1])  continue;
                //     for(k = 0; k < j; k++)
                //     {
                //         if(k == ID[0] || k == ID[1])    continue;
                //         if(lep_ID[j] == -lep_ID[k] && lep_ID[j] != 0 )
                //         {
                //             ID[2] = j;
                //             ID[3] = k;
                //             l_2 = lep[j]+lep[k];
                //             m_34 = l_2.M();
                //         }
                //     }
                // }

                if(ID[1] == 3 - ID[0])
                {
                    ID[2] = abs(2 - ID[0]);
                    ID[3] = 3 - ID[2];
                }
                else
                {
                    ID[2] = 3 - ID[0];
                    ID[3] = 3 - ID[1];
                }
                m_34 = (lep[ID[2]]+lep[ID[3]]).M();

                if(lep_ID[ID[2]] != -lep_ID[ID[3]]) continue;
                ossf[q]=ossf[q]+mcWeight;
                ossf_err[q]=ossf_err[q]+mcWeight*mcWeight;
                
                //cout<<lep_ID[ID[0]]<<lep_ID[ID[1]]<<lep_ID[ID[2]]<<lep_ID[ID[3]]<<endl;

                if(m_34 < 5000 ) continue;

                if(lep_ID[ID[0]] == -lep_ID[ID[2]])
                {
                    if((lep[ID[0]]+lep[ID[2]]).M() < 5000 || (lep[ID[1]]+lep[ID[3]]).M() < 5000 )  continue;
                    // if((lep[ID[0]]+lep[ID[2]]).M() < (lep[ID[1]]+lep[ID[3]]).M())   object = (lep[ID[0]]+lep[ID[2]]).M();
                    // else    object = (lep[ID[1]]+lep[ID[3]]).M();
                }
                if(lep_ID[ID[0]] == -lep_ID[ID[3]])
                {
                    if((lep[ID[0]]+lep[ID[3]]).M() < 5000 || (lep[ID[1]]+lep[ID[2]]).M() < 5000 )  continue;
                    // if((lep[ID[0]]+lep[ID[3]]).M() < (lep[ID[1]]+lep[ID[2]]).M())   object = (lep[ID[0]]+lep[ID[3]]).M();
                    // else    object = (lep[ID[1]]+lep[ID[2]]).M();
                }
                Jpsi[q]=Jpsi[q]+mcWeight;
                Jpsi_err[q]=Jpsi_err[q]+mcWeight*mcWeight;

                // if(nJets_OR > 0) continue;
                if(nJets_OR < 2) continue;
                jet_num[q] = jet_num[q] + mcWeight;
                jet_num_err[q] = jet_num_err[q] +mcWeight*mcWeight;

                // if(nJets_OR_DL1r_77 > 0)  continue;
                if(nJets_OR_DL1r_77 < 1)  continue;
                b_tag[q]=b_tag[q]+mcWeight;
                b_tag_err[q]=b_tag_err[q]+mcWeight*mcWeight;

                // if(m_12 > 75000 && m_12 < 100000) continue;
                // if(m_12 > 50000) continue;
                // if((jet[0]+jet[1]).M() < 280000) continue;
                // if(met_met > 80000)    continue;

                l_4=lep[ID[0]]+lep[ID[1]]+lep[ID[2]]+lep[ID[3]];

                if(l_4.M() < 107000 || l_4.M() > 133000) continue;
                // if(l_4.M() > 107000 && l_4.M() < 133000) continue;
                // if(l_4.M() < 300000) continue;
                // if(l_4.M() < 95000 || l_4.M() > 120000) continue;
                // if(l_4.M() > 120000) continue;

                m_llll[q]=m_llll[q]+mcWeight;
                m_llll_err[q]=m_llll_err[q]+mcWeight*mcWeight;

                TVector3 l_boost = l_4.BoostVector();
                TVector3 jet_boost = (jet[0]+jet[1]).BoostVector();
                TLorentzVector z_axis;
                z_axis.SetPxPyPzE(0,0,1,2);

                TLorentzVector l_0 = lep[ID[0]], l_1 = lep[ID[1]], Z_0 = lep[ID[0]]+lep[ID[1]], l_2 = lep[ID[2]], l_3 = lep[ID[3]], Z_1 = lep[ID[2]]+lep[ID[3]];
                
                l_0.Boost(-l_boost);
                l_1.Boost(-l_boost);
                Z_0.Boost(-l_boost);

                l_2.Boost(-l_boost);
                l_3.Boost(-l_boost);
                Z_1.Boost(-l_boost);
                
                //TLorentzVector jet_0 = jet[0].Boost(-jet_boost);
                //TLorentzVector jet_1 = jet[1].Boost(-jet_boost);
                
                cs_jet = VectorAngle(jet[0]+jet[1],jet[0]);

                if(lep_ID[ID[0]] < 0)  cs_lep_12 = VectorAngle(Z_0,l_0);
                else    cs_lep_12 = VectorAngle(Z_0,l_1);
                
                if(lep_ID[ID[2]] < 0)  cs_lep_34 = VectorAngle(Z_1,l_2);
                else    cs_lep_34 = VectorAngle(Z_1,l_3);

                //cout<<DecayAngle((lep[ID[0]]+lep[ID[1]]),lep[ID[0]])<<","<<DecayAngle((lep[ID[0]]+lep[ID[1]]),lep[ID[1]])<<endl;

                cs_Z_pair = PlaneAngle(z_axis,Z_0,l_0,l_1);

                cs_pairs = PlaneAngle(l_0,l_1,l_2,l_3);

                if(q != 4)
                {
                    lep_ID_0_f[q] = lep_ID[0];
                    lep_Pt_0_f[q] = lep[0].Pt();
                    lep_Eta_0_f[q] = lep[0].Eta();
                    lep_Phi_0_f[q] = lep[0].Phi();
                    lep_iso_0_f[q] = lep_isolationFCLoose_0*1.0;

                    lep_ID_1_f[q] = lep_ID[1];
                    lep_Pt_1_f[q] = lep[1].Pt();
                    lep_Eta_1_f[q] = lep[1].Eta();
                    lep_Phi_1_f[q] = lep[1].Phi();
                    lep_iso_1_f[q] = lep_isolationFCLoose_1*1.0;

                    lep_ID_2_f[q] = lep_ID[2];
                    lep_Pt_2_f[q] = lep[2].Pt();
                    lep_Eta_2_f[q] = lep[2].Eta();
                    lep_Phi_2_f[q] = lep[2].Phi();
                    lep_iso_2_f[q] = lep_isolationFCLoose_2*1.0;

                    lep_ID_3_f[q] = lep_ID[3];
                    lep_Pt_3_f[q] = lep[3].Pt();
                    lep_Eta_3_f[q] = lep[3].Eta();
                    lep_Phi_3_f[q] = lep[3].Phi();
                    lep_iso_3_f[q] = lep_isolationFCLoose_3*1.0;
                    
                    jet_E_0_f[q] = jet[0].E();
                    jet_Pt_0_f[q] = jet[0].Pt();
                    jet_Eta_0_f[q] = jet[0].Eta();
                    jet_Phi_0_f[q] = jet[0].Phi();
                    jet_E_1_f[q] = jet[1].E();
                    jet_Pt_1_f[q] = jet[1].Pt();
                    jet_Eta_1_f[q] = jet[1].Eta();
                    jet_Phi_1_f[q] = jet[1].Phi();

                    m_12_f[q] = m_12;
                    m_34_f[q] = m_34;
                    m_4l_f[q] = l_4.M()/1000;
                    m_jj_f[q] = (jet[0]+jet[1]).M()/1000;
                    p_4l_f[q] = l_4.Pt()/1000;
                    p_jj_f[q] = (jet[0]+jet[1]).Pt()/1000;
                    cs_jet_f[q] = cs_jet;
                    cs_lep_12_f[q] = cs_lep_12;
                    cs_lep_34_f[q] = cs_lep_34;
                    cs_pairs_f[q] = cs_pairs;
                    cs_Z_pair_f[q] = cs_Z_pair;
                    met_met_f[q] = met_met;
                    HT_f[q] = HT;
                    HT_lep_f[q] = HT_lep;
                    Dphi_met_jets_f[q] = met_phi - jet[0].Phi();

                    njets_f[q] = nJets_OR*1.0;
                    nbjets_f[q] = nJets_OR_DL1r_77*1.0;

                    mcWeight_f[q] = mcWeight;

                    tt[q]->Fill();
                }

                else if(q != N - 1 && q != N)
                {
                    lep_ID_0_f[N+1] = lep_ID[0];
                    lep_Pt_0_f[N+1] = lep[0].Pt();
                    lep_Eta_0_f[N+1] = lep[0].Eta();
                    lep_Phi_0_f[N+1] = lep[0].Phi();
                    lep_iso_0_f[N+1] = lep_isolationFCLoose_0*1.0;

                    lep_ID_1_f[N+1] = lep_ID[1];
                    lep_Pt_1_f[N+1] = lep[1].Pt();
                    lep_Eta_1_f[N+1] = lep[1].Eta();
                    lep_Phi_1_f[N+1] = lep[1].Phi();
                    lep_iso_1_f[N+1] = lep_isolationFCLoose_1*1.0;

                    lep_ID_2_f[N+1] = lep_ID[2];
                    lep_Pt_2_f[N+1] = lep[2].Pt();
                    lep_Eta_2_f[N+1] = lep[2].Eta();
                    lep_Phi_2_f[N+1] = lep[2].Phi();
                    lep_iso_2_f[N+1] = lep_isolationFCLoose_2*1.0;

                    lep_ID_3_f[N+1] = lep_ID[3];
                    lep_Pt_3_f[N+1] = lep[3].Pt();
                    lep_Eta_3_f[N+1] = lep[3].Eta();
                    lep_Phi_3_f[N+1] = lep[3].Phi();
                    lep_iso_3_f[N+1] = lep_isolationFCLoose_3*1.0;
                    
                    jet_E_0_f[N+1] = jet[0].E();
                    jet_Pt_0_f[N+1] = jet[0].Pt();
                    jet_Eta_0_f[N+1] = jet[0].Eta();
                    jet_Phi_0_f[N+1] = jet[0].Phi();
                    jet_E_1_f[N+1] = jet[1].E();
                    jet_Pt_1_f[N+1] = jet[1].Pt();
                    jet_Eta_1_f[N+1] = jet[1].Eta();
                    jet_Phi_1_f[N+1] = jet[1].Phi();

                    m_12_f[N+1] = m_12;
                    m_34_f[N+1] = m_34;
                    m_4l_f[N+1] = l_4.M()/1000;
                    m_jj_f[N+1] = (jet[0]+jet[1]).M()/1000;
                    p_4l_f[N+1] = l_4.Pt()/1000;
                    p_jj_f[N+1] = (jet[0]+jet[1]).Pt()/1000;
                    cs_jet_f[N+1] = cs_jet;
                    cs_lep_12_f[N+1] = cs_lep_12;
                    cs_lep_34_f[N+1] = cs_lep_34;
                    cs_pairs_f[N+1] = cs_pairs;
                    cs_Z_pair_f[N+1] = cs_Z_pair;
                    met_met_f[N+1] = met_met;
                    HT_f[N+1] = HT;
                    HT_lep_f[N+1] = HT_lep;
                    Dphi_met_jets_f[N+1] = met_phi - jet[0].Phi();

                    njets_f[N+1] = nJets_OR*1.0;
                    nbjets_f[N+1] = nJets_OR_DL1r_77*1.0;

                    mcWeight_f[N+1] = mcWeight;

                    bkg = q;

                    tt[N+1]->Fill();
                }

                // if(q == N - 1)
                // {
                    // lep_ID_0_s = lep_ID[ID[0]];
                    // lep_Pt_0_s = lep[ID[0]].Pt();
                    // lep_Eta_0_s = lep[ID[0]].Eta();
                    // lep_Phi_0_s = lep[ID[0]].Phi();
                    // lep_ID_1_s = lep_ID[ID[1]];
                    // lep_Pt_1_s = lep[ID[1]].Pt();
                    // lep_Eta_1_s = lep[ID[1]].Eta();
                    // lep_Phi_1_s = lep[ID[1]].Phi();
                    // lep_ID_2_s = lep_ID[ID[2]];
                    // lep_Pt_2_s = lep[ID[2]].Pt();
                    // lep_Eta_2_s = lep[ID[2]].Eta();
                    // lep_Phi_2_s = lep[ID[2]].Phi();
                    // lep_ID_3_s = lep_ID[ID[3]];
                    // lep_Pt_3_s = lep[ID[3]].Pt();
                    // lep_Eta_3_s = lep[ID[3]].Eta();
                    // lep_Phi_3_s = lep[ID[3]].Phi();

                //     lep_ID_0_s = lep_ID[0];
                //     lep_Pt_0_s = lep[0].Pt();
                //     lep_Eta_0_s = lep[0].Eta();
                //     lep_Phi_0_s = lep[0].Phi();
                //     lep_ID_1_s = lep_ID[1];
                //     lep_Pt_1_s = lep[1].Pt();
                //     lep_Eta_1_s = lep[1].Eta();
                //     lep_Phi_1_s = lep[1].Phi();
                //     lep_ID_2_s = lep_ID[2];
                //     lep_Pt_2_s = lep[2].Pt();
                //     lep_Eta_2_s = lep[2].Eta();
                //     lep_Phi_2_s = lep[2].Phi();
                //     lep_ID_3_s = lep_ID[3];
                //     lep_Pt_3_s = lep[3].Pt();
                //     lep_Eta_3_s = lep[3].Eta();
                //     lep_Phi_3_s = lep[3].Phi();
                    
                //     jet_E_0_s = jet[0].E();
                //     jet_Pt_0_s = jet[0].Pt();
                //     jet_Eta_0_s = jet[0].Eta();
                //     jet_Phi_0_s = jet[0].Phi();
                //     jet_E_1_s = jet[1].E();
                //     jet_Pt_1_s = jet[1].Pt();
                //     jet_Eta_1_s = jet[1].Eta();
                //     jet_Phi_1_s = jet[1].Phi();

                //     m_12_s = m_12;
                //     m_34_s = m_34;
                //     m_4l_s = l_4.M()/1000;
                //     m_jj_s = (jet[0]+jet[1]).M()/1000;
                //     p_4l_s = l_4.Pt()/1000;
                //     p_jj_s = (jet[0]+jet[1]).Pt()/1000;
                //     cs_jet_s = cs_jet;
                //     cs_lep_12_s = cs_lep_12;
                //     cs_lep_34_s = cs_lep_34;
                //     cs_pairs_s = cs_pairs;
                //     cs_Z_pair_s = cs_Z_pair;
                //     met_met_s = met_met;
                //     HT_s = HT;
                //     HT_lep_s = HT_lep;
                //     Dphi_met_jets_s = met_phi - jet[0].Phi();

                //     njets_s = nJets_OR*1.0;
                //     nbjets_s = nJets_OR_DL1r_77*1.0;

                //     mcWeight_s = mcWeight;

                //     t_sig->Fill();
                // }

                // else
                // {
                //     lep_ID_0_b = lep_ID[0];
                //     lep_Pt_0_b = lep[0].Pt();
                //     lep_Eta_0_b = lep[0].Eta();
                //     lep_Phi_0_b = lep[0].Phi();
                //     lep_ID_1_b = lep_ID[1];
                //     lep_Pt_1_b = lep[1].Pt();
                //     lep_Eta_1_b = lep[1].Eta();
                //     lep_Phi_1_b = lep[1].Phi();
                //     lep_ID_2_b = lep_ID[2];
                //     lep_Pt_2_b = lep[2].Pt();
                //     lep_Eta_2_b = lep[2].Eta();
                //     lep_Phi_2_b = lep[2].Phi();
                //     lep_ID_3_b = lep_ID[3];
                //     lep_Pt_3_b = lep[3].Pt();
                //     lep_Eta_3_b = lep[3].Eta();
                //     lep_Phi_3_b = lep[3].Phi();
                    
                //     jet_E_0_b = jet[0].E();
                //     jet_Pt_0_b = jet[0].Pt();
                //     jet_Eta_0_b = jet[0].Eta();
                //     jet_Phi_0_b = jet[0].Phi();
                //     jet_E_1_b = jet[1].E();
                //     jet_Pt_1_b = jet[1].Pt();
                //     jet_Eta_1_b = jet[1].Eta();
                //     jet_Phi_1_b = jet[1].Phi();

                //     m_12_b = m_12;
                //     m_34_b = m_34;
                //     m_4l_b = l_4.M()/1000;
                //     m_jj_b = (jet[0]+jet[1]).M()/1000;
                //     p_4l_b = l_4.Pt()/1000;
                //     p_jj_b = (jet[0]+jet[1]).Pt()/1000;
                //     cs_jet_b = cs_jet;
                //     cs_lep_12_b = cs_lep_12;
                //     cs_lep_34_b = cs_lep_34;
                //     cs_pairs_b = cs_pairs;
                //     cs_Z_pair_b = cs_Z_pair;
                //     met_met_b = met_met;
                //     HT_b = HT;
                //     HT_lep_b = HT_lep;
                //     Dphi_met_jets_b = met_phi - jet[0].Phi();

                //     njets_b = nJets_OR*1.0;
                //     nbjets_b = nJets_OR_DL1r_77*1.0;

                //     mcWeight_b = 25.62/4.4*mcWeight;

                //     bkg = q;

                //     t_bkg->Fill();
                // }
                
                // m_12 = 0;
                // for(j = 1; j < n_lep; j++)
                // {
                //     for(k = 0; k < j; k++)
                //     {
                //         l_2 = lep[j]+lep[k];
                //         if(l_2.M() > m_12 && lep_ID[j]*lep_ID[k] < 0)
                //         {
                //             m_12 = l_2.M();
                //             ID[0] = j;
                //             ID[1] = k;
                //         }
                //     }
                // }

                double variables[99] = {
                                        // m_12/1000,
                                        // (lep[2]+lep[3]).M()/1000,
                                        // lep_Phi_0,
                                        // lep_Phi_1,
                                        // lep_Phi_2,
                                        // lep_Phi_3,
                                        // lep_Eta_0,
                                        // lep_Eta_1,
                                        // lep_Eta_2,
                                        // lep_Eta_3,
                                        // lep_ID[ID[0]],//16
                                        // lep_ID[ID[2]],//16
                                        nJets_OR*1.0,//0
                                        nJets_OR_DL1r_77*1.0,//1
                                        (jet[0]+jet[1]).M()/1000.0,//2
                                        (jet[0]+jet[1]).Pt()/1000.0,//3
                                        l_4.M()/1000.0,//4
                                        l_4.Pt()/1000.0,//5
                                        m_12/1000.0,//6
                                        m_34/1000.0,//7
                                        cs_jet,//8
                                        cs_lep_12,//9
                                        cs_lep_34,//10
                                        met_met/1000.0,//11
                                        cs_pairs,//12
                                        cs_Z_pair,//13
                                        lep_sigd0PV_1,//14
                                        lep_Z0SinTheta_1,//15
                                        lep_isolationFCLoose_0*1.0,//16
                                        lep_isolationFCLoose_1*1.0,//17
                                        lep_isolationFCLoose_2*1.0,//18
                                        lep_isolationFCLoose_3*1.0,//19
                                        int(lep_isolationPflowLoose_0)*1.0,//20
                                        int(lep_isolationPflowLoose_1)*1.0,//21
                                        int(lep_isolationPflowLoose_2)*1.0,//22
                                        int(lep_isolationPflowLoose_3)*1.0,//23
                                        lep_ptVarcone30_0/lep_Pt_0,//24
                                        lep_ptVarcone30_1/lep_Pt_1,//25
                                        lep_ptVarcone30_2/lep_Pt_2,//26
                                        lep_ptVarcone30_3/lep_Pt_3,//27
                                        lep_topoEtcone30_0/lep_Pt_0,//28
                                        lep_topoEtcone30_1/lep_Pt_1,//29
                                        lep_topoEtcone30_2/lep_Pt_2,//30
                                        lep_topoEtcone30_3/lep_Pt_3,//31
                                        lepton[0].Phi()-jet[0].Phi(),//32
                                        lepton[1].Phi()-jet[0].Phi(),//33
                                        lepton[2].Phi()-jet[0].Phi(),//34
                                        lepton[3].Phi()-jet[0].Phi(),//35
                                        // lepton[0].Phi()-met_phi,//32
                                        // lepton[1].Phi()-met_phi,//33
                                        // lepton[2].Phi()-met_phi,//34
                                        // lepton[3].Phi()-met_phi,//35
                                        VectorAngle(lepton[0],jet[1]),//32
                                        VectorAngle(lepton[1],jet[1]),//33
                                        VectorAngle(lepton[2],jet[1]),//34
                                        VectorAngle(lepton[3],jet[1]),//35
                                        // mcWeight//32
                                        met_phi - jet[0].Phi()//24
                                        };

                for(p = 0; p < obj; p++)
                {
                    //h_err[p]->SetLineColor(kRed);
                    if(q == N-1)
                    {
                        hList[p][q]->Fill(variables[p],100*mcWeight);
                        hList[p][q]->SetLineColor(colorList[q]);
                        hList[p][q]->SetLineStyle(kDashed);
                    }
                    else
                    {
                        hList[p][q]->Fill(variables[p],mcWeight);
                        hList[p][q]->SetFillColor(colorList[q]);
                    }
                }

                events[q] = events[q] + 1.0;
                events_weighted[q] = events_weighted[q] + mcWeight;
                end_error[q] = end_error[q] + mcWeight*mcWeight;
                //ab_eff = ab_eff + pileupEventWeight_090*JVT_EventWeight*mcWeightOrg*scale_nom/mc_xSection/mc_kFactor;
                ab_eff = ab_eff + weight_pileup*weight_jvt*weight_mc*scale_nom/mc_xSection/mc_kFactor;
            }

            t->Delete("");

            cout<<q<<endl;
            cout<<"nentries:"<<nentries<<",events:"<<events[q] - cache[11]<<endl;
            cout<<"efficiency:"<<(events[q] - cache[11])*1.0/nentries<<endl;
            cout<<"trigger match:"<<trig_match[q] - cache[0]<<"+-"<<sqrt(trig_match_err[q] - err_cache[0])<<endl;
            cout<<"lep isolation:"<<iso[q] - cache[1]<<"+-"<<sqrt(iso_err[q] - err_cache[1])<<endl;
            cout<<"lep sepa:"<<lep_sepa[q] - cache[2]<<"+-"<<sqrt(lep_sepa_err[q] - err_cache[2])<<endl;
            cout<<"lep kine:"<<lep_pt[q] - cache[3]<<"+-"<<sqrt(lep_pt_err[q] - err_cache[3])<<endl;
            cout<<"pairs ossf:"<<ossf[q] - cache[4]<<"+-"<<sqrt(ossf_err[q] - err_cache[4])<<endl;
            cout<<"J/psi veto:"<<Jpsi[q] - cache[5]<<"+-"<<sqrt(Jpsi_err[q] - err_cache[5])<<endl;
            cout<<"jet number:"<<jet_num[q] - cache[6]<<"+-"<<sqrt(jet_num_err[q] - err_cache[6])<<endl;
            cout<<"b-tagging:"<<b_tag[q] - cache[7]<<"+-"<<sqrt(b_tag_err[q] - err_cache[7])<<endl;
            cout<<"4l mass:"<<m_llll[q] - cache[8]<<"+-"<<sqrt(m_llll_err[q] - err_cache[8])<<endl;
            cout<<"nentries_weighted:"<<nentries_weighted[q] - cache[9]<<"+-"<<sqrt(start_error[q] - err_cache[9])<<", events_weighted:"<<events_weighted[q] - cache[10]<<"+-"<<sqrt(end_error[q] - err_cache[10])<<endl;
            cout<<"efficiency_weighted:"<<events_weighted[q]*1.0/nentries_weighted[q]<<endl;
            cout<<"absoulute efficiency:"<<ab_eff<<endl;
            cout<<"Xsection:"<<mc_xSection<<endl;
            cout<<"raw_Xsection:"<<mc_rawXSection<<endl;
            //cout<<"total_weights:"<<mc_kFactor*mc_xSection/scale_nom<<",%f"<<totalWeights<<endl;
            Printf("total_weights:%.2f,%.2f\n",mc_kFactor*mc_xSection/scale_nom,totalWeights);
            //cout<<"check:"<<check<<",check_weighted:"<<check_weighted<<endl;
            
            myfile<<q<<endl;
            myfile<<"nentries:"<<nentries<<",events:"<<events[q] - cache[11]<<endl;
            myfile<<"efficiency:"<<(events[q] - cache[11])*1.0/nentries<<endl;
            myfile<<"trigger match:"<<trig_match[q] - cache[0]<<"+-"<<sqrt(trig_match_err[q] - err_cache[0])<<endl;
            myfile<<"lep isolation:"<<iso[q] - cache[1]<<"+-"<<sqrt(iso_err[q] - err_cache[1])<<endl;
            myfile<<"lep sepa:"<<lep_sepa[q] - cache[2]<<"+-"<<sqrt(lep_sepa_err[q] - err_cache[2])<<endl;
            myfile<<"lep kine:"<<lep_pt[q] - cache[3]<<"+-"<<sqrt(lep_pt_err[q] - err_cache[3])<<endl;
            myfile<<"pairs ossf:"<<ossf[q] - cache[4]<<"+-"<<sqrt(ossf_err[q] - err_cache[4])<<endl;
            myfile<<"J/psi veto:"<<Jpsi[q] - cache[5]<<"+-"<<sqrt(Jpsi_err[q] - err_cache[5])<<endl;
            myfile<<"jet number:"<<jet_num[q] - cache[6]<<"+-"<<sqrt(jet_num_err[q] - err_cache[6])<<endl;
            myfile<<"b-tagging:"<<b_tag[q] - cache[7]<<"+-"<<sqrt(b_tag_err[q] - err_cache[7])<<endl;
            myfile<<"4l mass:"<<m_llll[q] - cache[8]<<"+-"<<sqrt(m_llll_err[q] - err_cache[8])<<endl;
            myfile<<"nentries_weighted:"<<nentries_weighted[q] - cache[9]<<"+-"<<sqrt(start_error[q] - err_cache[9])<<", events_weighted:"<<events_weighted[q] - cache[10]<<"+-"<<sqrt(end_error[q] - err_cache[10])<<endl;
            // myfile<<"trigger match:"<<trig_match[q] - cache[0]<<"+-"<<sqrt(trig_match_err[q] - err_cache[0])<<endl;
            // myfile<<"lep isolation:"<<iso[q] - cache[2]<<"+-"<<sqrt(iso_err[q] - err_cache[2])<<endl;
            // myfile<<"lep sepa:"<<lep_sepa[q] - cache[3]<<"+-"<<sqrt(lep_sepa_err[q] - err_cache[3])<<endl;
            // myfile<<"lep kine:"<<lep_pt[q] - cache[4]<<"+-"<<sqrt(lep_pt_err[q] - err_cache[4])<<endl;
            // myfile<<"J/psi veto:"<<Jpsi[q] - cache[5]<<"+-"<<sqrt(Jpsi_err[q] - err_cache[5])<<endl;
            // myfile<<"jet number:"<<jet_num[q] - cache[6]<<"+-"<<sqrt(jet_num_err[q] - err_cache[6])<<endl;
            // myfile<<"b-tagging:"<<b_tag[q] - cache[7]<<"+-"<<sqrt(b_tag_err[q] - err_cache[7])<<endl;
            // myfile<<"4l mass:"<<m_llll[q] - cache[8]<<"+-"<<sqrt(m_llll_err[q] - err_cache[8])<<endl;
            // myfile<<"nentries_weighted:"<<nentries_weighted[q] - cache[9]<<"+-"<<sqrt(start_error[q] - err_cache[9])<<", events_weighted:"<<events_weighted[q] - cache[10]<<"+-"<<sqrt(end_error[q] - err_cache[10])<<endl;
            myfile<<"efficiency_weighted:"<<events_weighted[q]*1.0/nentries_weighted[q]<<endl;
            myfile<<"absoulute efficiency:"<<ab_eff<<endl;
            myfile<<"Xsection:"<<mc_xSection<<endl;
            myfile<<"raw_Xsection:"<<mc_rawXSection<<endl;
            myfile<<"total_weights:"<<mc_kFactor*mc_xSection/scale_nom<<","<<totalWeights<<endl;
            myfile<<endl;
            // Printf("total_weights:%.2f,%.2f\n",mc_kFactor*mc_xSection/scale_nom,totalWeights);
            //cout<<"check:"<<check<<",check_weighted:"<<check_weighted<<endl;
        }
    }

    pre_f->Write();
    pre_f->Close();

    TPaveText *pt;
    // pt = new TPaveText(.05,.1,.9,.8);
    // pt->SetBorderSize(0);
    // pt->SetFillColor(0);
    // pt->SetFillStyle(0);
    // pt->SetTextAlign(12);
    // pt->SetTextSize(0.03);
    // pt->AddText("tt VV ttV Higgs Z+jets VVV Signal");
    
    yieldfile.open (plotPath + "yields.txt");

    string flow[13];
    // flow = "Yields&Normalized Events&Trigger Matching&Lepton Isolation&Leptons Separation&Leptons Kinematics&Pair OSSF&J/$\\psi$ Veto&$\\ge2$ Jets&$\\ge1$ b Jets&107 GeV$<M_{4l}<$133 GeV&Final Entries Number&Relative Efficiency\\\\";
    
    flow[0] = "Yields";
    flow[1] = "Normalized Events";
    flow[2] = "Trigger Matching";
    flow[3] = "Lepton Isolation";
    flow[4] = "Leptons Separation";
    flow[5] = "Leptons Kinematics";
    flow[6] = "Pair OSSF";
    flow[7] = "J/$\\psi$ Veto";
    flow[8] = "$\\ge2$ Jets";
    flow[9] = "$\\ge1$ b Jets";
    flow[10] = "107 GeV$<M_{4l}<$133 GeV";
    flow[11] = "Final Entries Number";
    flow[12] = "Relative Efficiency";
    // yieldfile<<flow<<"\n"<<"\\midrule"<<endl;

    const char *flowList[13] = {"","Normalized Events","Trigger Matching","Leptons Isolation","Leptons Separation","Leptons Kinematics","Pair OSSF","J/ #psi veto","#geq 2 jets","#geq 1 b jets","107 GeV <M_{ 4l} <133 GeV","Final Entries Number","Relative Efficiency"};

    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 1600, 400);
    //TText *t = new TText();
    c1->SetLeftMargin(0.2);
    c1->SetRightMargin(0.2);
    c1->SetTopMargin(0.2);

    auto *lv = new TLine;
    // lv->SetLineColor(kBlue);
    lv->DrawLine(0,0.84,0.15+0.1*(N+2),0.84);
    lv->Draw();

    // lv->DrawLine(0.13,0.9-0.03*12,0.13,0.9);
    // lv->Draw();

    pt = new TPaveText(0.12,0.9-0.06*12,0.03,.9);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextFont(22);
    pt->SetTextSize(0.045);
    pt->SetTextAlign(22);

    for(p = 0; p < 13; p++) pt->AddText(flowList[p]);
    pt->Draw();

    for(p = 0; p < N+1 ; p++)
    {
        lv->DrawLine(0.1*p+0.16,0.9-0.06*12,0.1*p+0.16,0.9);
        lv->Draw();

        pt = new TPaveText(0.1*p+0.20,0.9-0.06*12,0.26+0.1*p,.9);
        pt->SetBorderSize(0);
        pt->SetFillColor(0);
        pt->SetFillStyle(0);
        pt->SetTextFont(22);
        pt->SetTextSize(0.045);
        pt->SetTextAlign(32);

        // entry = process_list[p];
        pt->AddText(process_list[p]);
        flow[0] = flow[0] + "&" + process_list[p];

        // pt->SetTextAlign(32);
        if(p != N)
        {
            flow[1] = flow[1] + "&" + (Form("%.3f$\\pm$%.3f",nentries_weighted[p],sqrt(start_error[p])));
            flow[2] = flow[2] + "&" + (Form("%.3f$\\pm$%.3f",trig_match[p],sqrt(trig_match_err[p])));
            flow[3] = flow[3] + "&" + (Form("%.3f$\\pm$%.3f",iso[p],sqrt(iso_err[p])));
            flow[4] = flow[4] + "&" + (Form("%.3f$\\pm$%.3f",lep_sepa[p],sqrt(lep_sepa_err[p])));
            flow[5] = flow[5] + "&" + (Form("%.3f$\\pm$%.3f",lep_pt[p],sqrt(lep_pt_err[p])));
            flow[6] = flow[6] + "&" + (Form("%.3f$\\pm$%.3f",ossf[p],sqrt(ossf_err[p])));
            flow[7] = flow[7] + "&" + (Form("%.3f$\\pm$%.3f",Jpsi[p],sqrt(Jpsi_err[p])));
            flow[8] = flow[8] + "&" + (Form("%.3f$\\pm$%.3f",jet_num[p],sqrt(jet_num_err[p])));
            flow[9] = flow[9] + "&" + (Form("%.3f$\\pm$%.3f",b_tag[p],sqrt(b_tag_err[p])));
            flow[10] = flow[10] + "&" + (Form("%.3f$\\pm$%.3f",m_llll[p],sqrt(m_llll_err[p])));
            flow[11] = flow[11] + "&" + (Form("%.0f",events[p]));
            flow[12] = flow[12] + "&" + (Form("%.4f",events_weighted[p]*1.0/nentries_weighted[p]));

            pt->AddText(Form("%.3f #pm%.3f",nentries_weighted[p],sqrt(start_error[p])));
            pt->AddText(Form("%.3f #pm%.3f",trig_match[p],sqrt(trig_match_err[p])));
            pt->AddText(Form("%.3f #pm%.3f",iso[p],sqrt(iso_err[p])));
            pt->AddText(Form("%.3f #pm%.3f",lep_sepa[p],sqrt(lep_sepa_err[p])));
            pt->AddText(Form("%.3f #pm%.3f",lep_pt[p],sqrt(lep_pt_err[p])));
            pt->AddText(Form("%.3f #pm%.3f",ossf[p],sqrt(ossf_err[p])));
            pt->AddText(Form("%.3f #pm%.3f",Jpsi[p],sqrt(Jpsi_err[p])));
            pt->AddText(Form("%.3f #pm%.3f",jet_num[p],sqrt(jet_num_err[p])));
            pt->AddText(Form("%.3f #pm%.3f",b_tag[p],sqrt(b_tag_err[p])));
            pt->AddText(Form("%.3f #pm%.3f",m_llll[p],sqrt(m_llll_err[p])));
            pt->AddText(Form("%.0f",events[p]));
            pt->AddText(Form("%.4f",events_weighted[p]*1.0/nentries_weighted[p]));
        }

        else
        {
            flow[1] = flow[1] + "&" + (Form("%.0f$\\pm$%.0f",nentries_weighted[p],sqrt(start_error[p])));
            flow[2] = flow[2] + "&" + (Form("%.0f$\\pm$%.0f",trig_match[p],sqrt(trig_match_err[p])));
            flow[3] = flow[3] + "&" + (Form("%.0f$\\pm$%.0f",iso[p],sqrt(iso_err[p])));
            flow[4] = flow[4] + "&" + (Form("%.0f$\\pm$%.0f",lep_sepa[p],sqrt(lep_sepa_err[p])));
            flow[5] = flow[5] + "&" + (Form("%.0f$\\pm$%.0f",lep_pt[p],sqrt(lep_pt_err[p])));
            flow[6] = flow[6] + "&" + (Form("%.0f$\\pm$%.0f",ossf[p],sqrt(ossf_err[p])));
            flow[7] = flow[7] + "&" + (Form("%.0f$\\pm$%.0f",Jpsi[p],sqrt(Jpsi_err[p])));
            flow[8] = flow[8] + "&" + (Form("%.0f$\\pm$%.0f",jet_num[p],sqrt(jet_num_err[p])));
            flow[9] = flow[9] + "&" + (Form("%.0f$\\pm$%.0f",b_tag[p],sqrt(b_tag_err[p])));
            flow[10] = flow[10] + "&" + (Form("%.0f$\\pm$%.0f",m_llll[p],sqrt(m_llll_err[p])));
            flow[11] = flow[11] + "&" + (Form("%.0f",events[p]));
            flow[12] = flow[12] + "&" + (Form("%.4f",events_weighted[p]*1.0/nentries_weighted[p]));

            pt->AddText(Form("%.0f",nentries_weighted[p]));
            pt->AddText(Form("%.0f",trig_match[p]));
            pt->AddText(Form("%.0f",iso[p]));
            pt->AddText(Form("%.0f",lep_sepa[p]));
            pt->AddText(Form("%.0f",lep_pt[p]));
            pt->AddText(Form("%.0f",ossf[p]));
            pt->AddText(Form("%.0f",Jpsi[p]));
            pt->AddText(Form("%.0f",jet_num[p]));
            pt->AddText(Form("%.0f",b_tag[p]));
            pt->AddText(Form("%.0f",m_llll[p]));
            pt->AddText(Form("%.0f",events[p]));
            pt->AddText(Form("%.4f",events_weighted[p]*1.0/nentries_weighted[p]));
        }
        // pt->AddText(Form("%.3f#pm%.3f %.3f#pm%.3f %.3f#pm%.3f %.3f#pm%.3f %.3f#pm%.3f %.3f#pm%.3f %.3f#pm%.3f %.3f#pm%.3f %.3f#pm%.3f %d",nentries_weighted[p],start_error[p],trig_match[p],trig_match_err[p],iso[p],iso_err[p],lep_sepa[p],lep_sepa_err[p],lep_pt[p],lep_pt_err[p],Jpsi[p],Jpsi_err[p],jet_num[p],jet_num_err[p],b_tag[p],b_tag_err[p],m_llll[p],m_llll_err[p],events[p]));

        pt->Draw();

        // flow = process_list[p] + "," + Form("%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d,%.5f",nentries_weighted[p],trig_match[p],iso[p],lep_sepa[p],lep_pt[p],Jpsi[p],jet_num[p],b_tag[p],m_llll[p],events[p],events_weighted[p]*1.0/nentries_weighted[p]);

        // flow = process_list[p];
        // flow = flow + "&" + (Form("%.3f$\\pm$%.3f&%.3f$\\pm$%.3f&%.3f$\\pm$%.3f&%.3f$\\pm$%.3f&%.3f$\\pm$%.3f&%.3f$\\pm$%.3f&%.3f$\\pm$%.3f&%.3f$\\pm$%.3f&%.3f$\\pm$%.3f&%.0f\\\\",nentries_weighted[p],start_error[p],trig_match[p],trig_match_err[p],lep_sepa[p],lep_sepa_err[p],lep_pt[p],lep_pt_err[p],iso[p],iso_err[p],Jpsi[p],Jpsi_err[p],jet_num[p],jet_num_err[p],b_tag[p],b_tag_err[p],m_llll[p],m_llll_err[p],events[p]));

        cout<<p<<endl;
        cout<<"Final events:"<<events[p]<<endl;
        cout<<"trigger match:"<<trig_match[p]<<"+-"<<sqrt(trig_match_err[p])<<endl;
        cout<<"lep isolation:"<<iso[p]<<"+-"<<sqrt(iso_err[p])<<endl;
        cout<<"lep sepa:"<<lep_sepa[p]<<"+-"<<sqrt(lep_sepa_err[p])<<endl;
        cout<<"lep kine:"<<lep_pt[p]<<"+-"<<sqrt(lep_pt_err[p])<<endl;
        cout<<"pair ossf:"<<ossf[p]<<"+-"<<sqrt(ossf_err[p])<<endl;
        cout<<"J/psi veto:"<<Jpsi[p]<<"+-"<<sqrt(Jpsi_err[p])<<endl;
        cout<<"jet number:"<<jet_num[p]<<"+-"<<sqrt(jet_num_err[p])<<endl;
        cout<<"b-tagging:"<<b_tag[p]<<"+-"<<sqrt(b_tag_err[p])<<endl;
        cout<<"4l mass:"<<m_llll[p]<<"+-"<<sqrt(m_llll_err[p])<<endl;
        cout<<"nentries_weighted:"<<nentries_weighted[p]<<"+-"<<sqrt(start_error[p])<<",events_weighted:"<<events_weighted[p]<<"+-"<<sqrt(end_error[p])<<endl;
        cout<<"efficiency_weighted:"<<events_weighted[p]*1.0/nentries_weighted[p]<<endl;
        cout<<endl;

        myfile<<p<<endl;
        myfile<<"Final events:"<<events[p]<<endl;
        myfile<<"trigger match:"<<trig_match[p]<<"+-"<<sqrt(trig_match_err[p])<<endl;
        myfile<<"lep isolation:"<<iso[p]<<"+-"<<sqrt(iso_err[p])<<endl;
        myfile<<"lep sepa:"<<lep_sepa[p]<<"+-"<<sqrt(lep_sepa_err[p])<<endl;
        myfile<<"lep kine:"<<lep_pt[p]<<"+-"<<sqrt(lep_pt_err[p])<<endl;
        myfile<<"pair ossf:"<<ossf[p]<<"+-"<<sqrt(ossf_err[p])<<endl;
        myfile<<"J/psi veto:"<<Jpsi[p]<<"+-"<<sqrt(Jpsi_err[p])<<endl;
        myfile<<"jet number:"<<jet_num[p]<<"+-"<<sqrt(jet_num_err[p])<<endl;
        myfile<<"b-tagging:"<<b_tag[p]<<"+-"<<sqrt(b_tag_err[p])<<endl;
        myfile<<"4l mass:"<<m_llll[p]<<"+-"<<sqrt(m_llll_err[p])<<endl;
        myfile<<"nentries_weighted:"<<nentries_weighted[p]<<"+-"<<sqrt(start_error[p])<<",events_weighted:"<<events_weighted[p]<<"+-"<<sqrt(end_error[p])<<endl;
        myfile<<"efficiency_weighted:"<<events_weighted[p]*1.0/nentries_weighted[p]<<endl;
        myfile<<endl;

        if(p < N - 1)   
        {
            B = B + events_weighted[p];
            B_error = B_error + end_error[p]*end_error[p];
            for(q = 0; q < obj; q++)
            {
                hs[q]->Add(hList[q][p]);
                h_err[q]->Add(hList[q][p]);
            }
        }
    
        else if(p == N-1) S = events_weighted[p];
        else    data = events_weighted[p];
    }

    for(p = 0; p < 13; p++)  yieldfile<<flow[p]<<"\\\\"<<"\n"<<"\\midrule"<<endl;

    c1->SaveAs(plotPath+"yields.svg");
    c1->SaveAs(plotPath+"yields.png");
    c1->SaveAs(plotPath+"yields.pdf");
    // c1->SaveAs(plotPath+"yields.txt");
    c1->Clear();

    cout<<"S:"<<S<<",B:"<<B<<",B error:"<<sqrt(B_error)<<",S/sqrt(B):"<<S/sqrt(B)<<",data:"<<data<<endl;
    myfile<<"S:"<<S<<",B:"<<B<<",B error:"<<B_error<<",S/sqrt(B):"<<S/sqrt(B)<<",data:"<<data<<endl;
    //cout<<",poi"<<sqrt(2*((S+B)log(1+S/B)-S))<<endl;

    TCanvas *c2 = new TCanvas("c2", "Graph Draw Options", 3600, 1600);
    //TText *t = new TText();
    c2->SetLeftMargin(0.25);
    c2->SetRightMargin(0.25);
    c2->SetTopMargin(0.25);
    c2->Divide(4,2,0.01,0.01);

    float scale;
    TH1* frame;
    TText *text;
    TLegend *leg;

    TRatioPlot *rp;

    //cout<<type<<endl;
    //cout<<categories[type]<<endl;
    gStyle->SetOptStat("");
    for(q = 0; q < obj; q++)
    {
        c2->cd(q+1-8*int(q/8));
        // cout<<q+1-8*int(q/8)<<endl;

        scale = (hList[q][N]->GetMaximum())*3.0/2;
        if(scale == 0)  scale = (hs[q]->GetMaximum())*3.0/2;    // No data
        frame = gPad->DrawFrame(xrange[q][0],1e-5,xrange[q][1],scale);
        gPad->Update();
        frame->GetXaxis()->SetNdivisions(505, kTRUE);
        frame->GetYaxis()->SetTitle("Events");
        frame->GetXaxis()->SetTitle(xaxis_name[q]);
        // c1->SetLogy();

        // frame = gPad->DrawFrame(xrange[q][0],1e-5,xrange[q][1],scale);
        // gPad->Update();

        // hList[q][N]->GetXaxis()->SetLimits(xrange[q][0],xrange[q][1]);
        // hList[q][N]->GetYaxis()->SetRangeUser(1e-5,scale);
        // hList[q][N]->GetXaxis()->SetNdivisions(505, kTRUE);
        // // hList[q][N]->GetYaxis()->SetTitle("Events");
        // hList[q][N]->GetXaxis()->SetTitle(xaxis_name[q]);

        // rp = new TRatioPlot(hList[q][N]);
        // rp->GetLowerRefYaxis()->SetTitle("Ratio");
        // rp->GetUpperRefYaxis()->SetTitle("Events");
        // rp->GetLowerRefGraph()->SetMinimum(1.5);
        // rp->GetLowerRefGraph()->SetMaximum(0.5);

        // rp->Draw();

        // h_err[q]->Draw("same hist, E1");

        pt=new TPaveText(0.33,0.92,0.23,0.70,"NDC");
        pt->SetBorderSize(0);
        pt->SetFillColor(0);
        pt->SetFillStyle(0);
        pt->SetTextAlign(12);
        pt->SetTextSize(0.03);

        text=pt->AddText("#scale[1.3]{#it{ATLAS} internal}");
        text=pt->AddText("HH#rightarrow b#bar{b}ZZ#rightarrow b#bar{b}4l");
        text=pt->AddText("13 TeV, 139 fb^{-1}");
        text=pt->AddText(categories[type]);

        pt->Draw("same");
        
        leg = new TLegend(0.92,0.92,0.7,0.68);
        leg->SetFillColor(0);
        leg->SetLineColor(0);
        leg->SetTextSize(0.03);
        leg->AddEntry(hList[q][0],"t#bar{t}","f");
        leg->AddEntry(hList[q][1],"VV","f");
        leg->AddEntry(hList[q][2],"ttV","f");
        leg->AddEntry(hList[q][3],"Single H","f");
        leg->AddEntry(hList[q][4],"Z+jets","f");
        // leg->AddEntry(hList[q][5],"VVV","f");
        leg->AddEntry(hList[q][N-1],"Signal*100","f");
        leg->AddEntry(hList[q][N],"Run-2 Data","lep");
        leg->Draw("same");
        
        hs[q]->Draw("same hist");
        hList[q][N-1]->Draw("same hist");
        hList[q][N]->Draw("same E1");

        if( (q+1)%8 == 0 )
        {
            c2->SaveAs(plotPath + (q+1)/8 +"_dist.c");
            c2->SaveAs(plotPath + (q+1)/8 +"_dist.png");
            c2->SaveAs(plotPath + (q+1)/8 +"_dist.svg");
            c2->SaveAs(plotPath + (q+1)/8 +"_dist.pdf");
            c2->Clear();
            c2->Divide(4,2,0.01,0.01);
        }
        // c2->SaveAs(plotPath+plotname[q]+"_dist.png");
        // c2->SaveAs(plotPath+plotname[q]+"_dist.svg");
        // c2->Clear();
    }
    c2->SaveAs(plotPath + ((q+1)/8+1) + "_dist.c");
    c2->SaveAs(plotPath + ((q+1)/8+1) + "_dist.png");
    c2->SaveAs(plotPath + ((q+1)/8+1) + "_dist.eps");
    c2->SaveAs(plotPath + ((q+1)/8+1) + "_dist.pdf");
    c2->Clear();
}

float VectorAngle(TLorentzVector Init, TLorentzVector Fin)
{
    /*float b,r,p_x,p_y,p_xx,cs;

    r = Init.E()/Init.M();
    b = sqrt(1-1/r/r);
    p_x = (Fin.Px()*Init.Px()+Fin.Py()*Init.Py()+Fin.Pz()*Init.Pz())/sqrt(Init.Px()*Init.Px()+Init.Py()*Init.Py()+Init.Pz()*Init.Pz());
    p_y = sqrt(Fin.Px()*Fin.Px()+Fin.Py()*Fin.Py()+Fin.Pz()*Fin.Pz()-p_x*p_x);
    p_xx = r*(p_x - b*Fin.E());
    cs = p_xx/sqrt(p_xx*p_xx+p_y*p_y);*/
    //cout<<p_x<<","<<p_xx<<","<<p_y<<","<<cs<<endl;
    float cs = 0.0;

    cs = (Fin.Px()*Init.Px()+Fin.Py()*Init.Py()+Fin.Pz()*Init.Pz())/sqrt(Init.Px()*Init.Px()+Init.Py()*Init.Py()+Init.Pz()*Init.Pz())/sqrt(Fin.Px()*Fin.Px()+Fin.Py()*Fin.Py()+Fin.Pz()*Fin.Pz());

    return cs;
}

float PlaneAngle(TLorentzVector l_0,TLorentzVector l_1,TLorentzVector l_2,TLorentzVector l_3)
{
    float cs = 0.0;

    TVector3 l_01,l_23;
    l_01.SetXYZ(l_0.Py()*l_1.Pz()-l_1.Py()*l_0.Pz(), l_1.Px()*l_0.Pz()-l_0.Px()*l_1.Pz(), l_1.Py()*l_0.Px()-l_0.Py()*l_1.Px());
    l_23.SetXYZ(l_2.Py()*l_3.Pz()-l_3.Py()*l_2.Pz(), l_3.Px()*l_2.Pz()-l_2.Px()*l_3.Pz(), l_3.Py()*l_2.Px()-l_2.Py()*l_3.Px());

    cs = l_01*l_23/sqrt((l_01*l_01)*(l_23*l_23));

    return cs;
}

int main(int argc, char **argv)
{
    if (argc == 1 || argc > 2)
        cout << "Usage: " << argv[0] << " dataPath ";
    else
    {
        TString dataPath(argv[1]);
        if (!dataPath.EndsWith("/"))
            dataPath = dataPath + "/";
        // string dataPath (argv[1]);
        cout << dataPath << endl;
        cutflow(dataPath,n_type);
    }
    return 0;
}
