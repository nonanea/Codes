#include<math.h>

const int N = 7;

// void input(TString treename)
void input()
{
  gErrorIgnoreLevel = kError;

  // This loads the library
  TMVA::Tools::Instance();

  // Create the Reader object
  // : This will evaluate our data according to TMVA machine learning result
  TMVA::Reader *reader_s = new TMVA::Reader( "!Color:!Silent" );
  // TMVA::Reader *reader_b = new TMVA::Reader( "!Color:!Silent" );
  // TMVA::Reader *reader_d = new TMVA::Reader( "!Color:!Silent" );

  // Create a set of variables and declare them to the reader

  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
  Float_t lep_Pt_0_s, lep_Pt_1_s, lep_Pt_2_s, lep_Pt_3_s, lep_Etcone30_0_s, lep_Etcone30_1_s, lep_Etcone30_2_s, lep_Etcone30_3_s, lep_Eta_0_s, lep_Eta_1_s, lep_Eta_2_s, lep_Eta_3_s, jet_E_0_s, jet_E_1_s, jet_Pt_0_s, jet_Pt_1_s, jet_Phi_0_s, jet_Phi_1_s, m_12_s, m_34_s, p_4l_s, p_jj_s, m_4l_s, m_jj_s, cs_pairs_s, cs_Z_pair_s, cs_lep_12_s, cs_lep_34_s, met_met_s, HT_s, HT_lep_s, Dphi_met_jets_s, njets_s, nbjets_s;

  Float_t lep_Pt_0_b, lep_Pt_1_b, lep_Pt_2_b, lep_Pt_3_b, lep_Phi_0_b, lep_Phi_1_b, lep_Phi_2_b, lep_Phi_3_b, jet_E_0_b, jet_E_1_b, jet_Pt_0_b, jet_Pt_1_b, jet_Phi_0_b, jet_Phi_1_b, m_12_b, m_34_b, p_4l_b, p_jj_b, m_4l_b, m_jj_b, cs_pairs_b, cs_Z_pair_b, cs_lep_12_b, cs_lep_34_b, met_met_b, HT_b, HT_lep_b, Dphi_met_jets_b, njets_b, nbjets_b;
  Float_t lep_Pt_0_d, lep_Pt_1_d, lep_Pt_2_d, lep_Pt_3_d, lep_Phi_0_d, lep_Phi_1_d, lep_Phi_2_d, lep_Phi_3_d, jet_E_0_d, jet_E_1_d, jet_Pt_0_d, jet_Pt_1_d, jet_Phi_0_d, jet_Phi_1_d, m_12_d, m_34_d, p_4l_d, p_jj_d, m_4l_d, m_jj_d, cs_pairs_d, cs_Z_pair_d, cs_lep_12_d, cs_lep_34_d, met_met_d, HT_d, HT_lep_d, Dphi_met_jets_d, njets_d, nbjets_d;
  Float_t weight_s, weight_b;

  Float_t bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1down,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1down,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1down,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1down,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1down,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1down,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1down,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1down,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1down,bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1down,bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1down,lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_DOWN_AT,lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_DOWN_AT,lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_DOWN_AT,lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_DOWN_AT,lepSF_SF_El_Reco_DOWN_AT,lepSF_SF_El_ID_LooseAndBLayerLH_DOWN_AT,lepSF_SF_El_ID_TightLH_DOWN_AT,lepSF_SF_El_Iso_FCLoose_DOWN_AT,lepSF_SF_El_PLVLoose_DOWN_AT,lepSF_SF_El_PLVTight_DOWN_AT,lepSF_SF_Mu_TTVA_STAT_DOWN_AT,lepSF_SF_Mu_TTVA_SYST_DOWN_AT,lepSF_SF_Mu_ID_Loose_STAT_DOWN_AT,lepSF_SF_Mu_ID_Medium_STAT_DOWN_AT,lepSF_SF_Mu_ID_Loose_SYST_DOWN_AT,lepSF_SF_Mu_ID_Medium_SYST_DOWN_AT,lepSF_SF_Mu_ID_Loose_STAT_LOWPT_DOWN_AT,lepSF_SF_Mu_ID_Medium_STAT_LOWPT_DOWN_AT,lepSF_SF_Mu_ID_Loose_SYST_LOWPT_DOWN_AT,lepSF_SF_Mu_ID_Medium_SYST_LOWPT_DOWN_AT,lepSF_SF_Mu_Iso_FCLoose_SYST_DOWN_AT,lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_down_AT,lepSF_SF_Mu_PLVLoose_DOWN_AT,lepSF_SF_Mu_PLVTight_DOWN_AT,custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down,custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1down,custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1down,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1down,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1down;

  Float_t bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1up,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1up,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1up,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1up,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1up,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1up,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1up,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1up,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1up,bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1up,bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1up,lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_UP_AT,lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_UP_AT,lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_UP_AT,lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_UP_AT,lepSF_SF_El_Reco_UP_AT,lepSF_SF_El_ID_LooseAndBLayerLH_UP_AT,lepSF_SF_El_ID_TightLH_UP_AT,lepSF_SF_El_Iso_FCLoose_UP_AT,lepSF_SF_El_PLVLoose_UP_AT,lepSF_SF_El_PLVTight_UP_AT,lepSF_SF_Mu_TTVA_STAT_UP_AT,lepSF_SF_Mu_TTVA_SYST_UP_AT,lepSF_SF_Mu_ID_Loose_STAT_UP_AT,lepSF_SF_Mu_ID_Medium_STAT_UP_AT,lepSF_SF_Mu_ID_Loose_SYST_UP_AT,lepSF_SF_Mu_ID_Medium_SYST_UP_AT,lepSF_SF_Mu_ID_Loose_STAT_LOWPT_UP_AT,lepSF_SF_Mu_ID_Medium_STAT_LOWPT_UP_AT,lepSF_SF_Mu_ID_Loose_SYST_LOWPT_UP_AT,lepSF_SF_Mu_ID_Medium_SYST_LOWPT_UP_AT,lepSF_SF_Mu_Iso_FCLoose_SYST_UP_AT,lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_UP_AT,lepSF_SF_Mu_PLVLoose_UP_AT,lepSF_SF_Mu_PLVTight_UP_AT,
  custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up,custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1up,custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1up,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1up,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1up;

  Float_t sum_1, sum_2;

  reader_s->AddVariable( "lep_Pt_0:= lep_Pt_0", &lep_Pt_0_s);
  reader_s->AddVariable( "lep_Pt_1:= lep_Pt_1", &lep_Pt_1_s);
  reader_s->AddVariable( "lep_Pt_2:= lep_Pt_2", &lep_Pt_2_s);
  reader_s->AddVariable( "lep_Pt_3:= lep_Pt_3", &lep_Pt_3_s);
  reader_s->AddVariable( "lep_Etcone30_0:= lep_Etcone30_0", &lep_Etcone30_0_s);
  reader_s->AddVariable( "lep_Etcone30_1:= lep_Etcone30_1", &lep_Etcone30_1_s);
  reader_s->AddVariable( "lep_Etcone30_2:= lep_Etcone30_2", &lep_Etcone30_2_s);
  reader_s->AddVariable( "lep_Etcone30_3:= lep_Etcone30_3", &lep_Etcone30_3_s);
  // reader_s->AddVariable( "lep_Eta_0:= lep_Eta_0", &lep_Eta_0_s);
  reader_s->AddVariable( "lep_Eta_1:= lep_Eta_1", &lep_Eta_1_s);
  reader_s->AddVariable( "lep_Eta_2:= lep_Eta_2", &lep_Eta_2_s);
  reader_s->AddVariable( "lep_Eta_3:= lep_Eta_3", &lep_Eta_3_s);
  // reader_s->AddVariable( "jet_E_0:= jet_E_0", &jet_E_0_s);
  // reader_s->AddVariable( "jet_E_1:= jet_E_1", &jet_E_1_s);
  reader_s->AddVariable( "jet_Pt_0:= jet_Pt_0", &jet_Pt_0_s);
  // reader_s->AddVariable( "jet_Pt_1:= jet_Pt_1", &jet_Pt_1_s);
  // reader_s->AddVariable( "jet_Phi_0:= jet_Phi_0", &jet_Phi_0_s);
  // reader_s->AddVariable( "jet_Phi_1:= jet_Phi_1", &jet_Phi_1_s);
  reader_s->AddVariable( "m_12:= m_12", &m_12_s);
  reader_s->AddVariable( "m_34:= m_34", &m_34_s);
  reader_s->AddVariable( "m_4l:= m_4l", &m_4l_s);
  reader_s->AddVariable( "m_jj:= m_jj", &m_jj_s);
  // reader_s->AddVariable( "p_4l:= p_4l", &p_4l_s);
  reader_s->AddVariable( "p_jj:= p_jj", &p_jj_s);
  reader_s->AddVariable( "HT:= HT", &HT_s);
  // reader_s->AddVariable( "HT_lep:= HT_lep", &HT_lep_s);
  // reader_s->AddVariable( "cs_lep_12:= cs_lep_12", &cs_lep_12_s);
  // reader_s->AddVariable( "cs_lep_34:= cs_lep_34", &cs_lep_34_s);
  // reader_s->AddVariable( "cs_pairs:= cs_pairs", &cs_pairs_s);
  reader_s->AddVariable( "met_met:= met_met", &met_met_s);
  reader_s->AddVariable( "Dphi_met_jets:= Dphi_met_jets", &Dphi_met_jets_s);
  //reader_s->AddVariable( "cs_Z_pair:= cs_Z_pair", &cs_Z_pair_s);
  // reader_s->AddVariable( "njets:= njets", &njets_s);
  reader_s->AddVariable( "nbjets:= nbjets", &nbjets_s);

  // We are going to evaluate only with MLP method
  // reader_s->BookMVA( "BDT method", "dataset/weights/TMVAClassification_BDT.weights.xml" );
  // reader_b->BookMVA( "BDT method", "dataset/weights/TMVAClassification_BDT.weights.xml" );
  // reader_d->BookMVA( "BDT method", "dataset/weights/TMVAClassification_BDT.weights.xml" );

  reader_s->BookMVA("BDT method", "dataset/weights/nominal_BDT.weights.xml" );
  
  // output
  float bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1down_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1down_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1down_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1down_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1down_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1down_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1down_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1down_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1down_f,bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1down_f,bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1down_f,lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_DOWN_AT_f,lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_DOWN_AT_f,lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_DOWN_AT_f,lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_DOWN_AT_f,lepSF_SF_El_Reco_DOWN_AT_f,lepSF_SF_El_ID_LooseAndBLayerLH_DOWN_AT_f,lepSF_SF_El_ID_TightLH_DOWN_AT_f,lepSF_SF_El_Iso_FCLoose_DOWN_AT_f,lepSF_SF_El_PLVLoose_DOWN_AT_f,lepSF_SF_El_PLVTight_DOWN_AT_f,lepSF_SF_Mu_TTVA_STAT_DOWN_AT_f,lepSF_SF_Mu_TTVA_SYST_DOWN_AT_f,lepSF_SF_Mu_ID_Loose_STAT_DOWN_AT_f,lepSF_SF_Mu_ID_Medium_STAT_DOWN_AT_f,lepSF_SF_Mu_ID_Loose_SYST_DOWN_AT_f,lepSF_SF_Mu_ID_Medium_SYST_DOWN_AT_f,lepSF_SF_Mu_ID_Loose_STAT_LOWPT_DOWN_AT_f,lepSF_SF_Mu_ID_Medium_STAT_LOWPT_DOWN_AT_f,lepSF_SF_Mu_ID_Loose_SYST_LOWPT_DOWN_AT_f,lepSF_SF_Mu_ID_Medium_SYST_LOWPT_DOWN_AT_f,lepSF_SF_Mu_Iso_FCLoose_SYST_DOWN_AT_f,lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_down_AT_f,lepSF_SF_Mu_PLVLoose_DOWN_AT_f,lepSF_SF_Mu_PLVTight_DOWN_AT_f,custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down_f,custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1down_f,custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1down_f,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down_f,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1down_f,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1down_f;
  
  float bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1up_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1up_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1up_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1up_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1up_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1up_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1up_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1up_f,bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1up_f,bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1up_f,bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1up_f,lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_UP_AT_f,lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_UP_AT_f,lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_UP_AT_f,lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_UP_AT_f,lepSF_SF_El_Reco_UP_AT_f,lepSF_SF_El_ID_LooseAndBLayerLH_UP_AT_f,lepSF_SF_El_ID_TightLH_UP_AT_f,lepSF_SF_El_Iso_FCLoose_UP_AT_f,lepSF_SF_El_PLVLoose_UP_AT_f,lepSF_SF_El_PLVTight_UP_AT_f,lepSF_SF_Mu_TTVA_STAT_UP_AT_f,lepSF_SF_Mu_TTVA_SYST_UP_AT_f,lepSF_SF_Mu_ID_Loose_STAT_UP_AT_f,lepSF_SF_Mu_ID_Medium_STAT_UP_AT_f,lepSF_SF_Mu_ID_Loose_SYST_UP_AT_f,lepSF_SF_Mu_ID_Medium_SYST_UP_AT_f,lepSF_SF_Mu_ID_Loose_STAT_LOWPT_UP_AT_f,lepSF_SF_Mu_ID_Medium_STAT_LOWPT_UP_AT_f,lepSF_SF_Mu_ID_Loose_SYST_LOWPT_UP_AT_f,lepSF_SF_Mu_ID_Medium_SYST_LOWPT_UP_AT_f,lepSF_SF_Mu_Iso_FCLoose_SYST_UP_AT_f,lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_UP_AT_f,lepSF_SF_Mu_PLVLoose_UP_AT_f,lepSF_SF_Mu_PLVTight_UP_AT_f,custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up_f,custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1up_f,custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1up_f,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up_f,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1up_f,custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1up_f;

  Int_t i, j, k, l, bkg, dsid, dsid_s, nbin = 8;
  Double_t range[2] = {-1.0,1.0};
  Double_t bin_edges[17] = {-0.6,-0.5,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.3,0.4};
  Double_t ZH[4][2] = {{0.35,1.97959887981},{0.55,1.28266811371},{0.75,1.24872112274},{0.95,0.531870007515}};
  Double_t ZH_sum = ZH[0][1] + ZH[1][1] + ZH[2][1] + ZH[3][1];

  // f_s = "Ntuples/signal/" + process_list[i] + "_Preselection.root";
  // f_b = "Ntuples/background/" + process_list[i] + "_Preselection.root";
  // histname = file.ReplaceAll("_Preselection.root","") + "BDT";

  // Prepare input tree for application
  // TFile *input = TFile::Open( filename,"read" );
  const char *process_list[N+1] = {"tt","VV","ttV","Higgs","Zjets","ggF","VBF","data"};
  const char *region[N] = {"tt","ttV","VVHiggs","Zjets","SR","VR"};
  // const char *region[N] = {"SR"};

  const char* entry;
  void *file_dir;
  TString file_path, output_path, ntuple, histname, syst_name;
  TFile *f, *output;
  TTree *tree, *tBDT;
  TH1F *histBDT;
  // for(i = 0; i < N+1; i++)  histBDT[i] = new TH1F("BDT",process_list[i],nbin,range[0],range[1]);

  Float_t bdt_eval_s = 0.0;
  // double sum;
  k = 4;
  // for(k = 0; k == 0; k++)
  // for(k = 0; k < N; k++)
  {
    for(i = 0; i < N+1; i++)
    // i = 3;
    {
      file_path = "Ntuples/";
      // file_path = "test/";
      file_path = file_path + region[k] + "/" + process_list[i] + "/";
      
      output_path = "Results/";
      output = new TFile(output_path + region[k] + "/" + process_list[i]+".root","recreate");
      output->Close();

      file_dir = gSystem->OpenDirectory(file_path);
      while((entry = (char *)gSystem->GetDirEntry(file_dir)))
      {
        ntuple = entry;

        // if(ntuple.EndsWith("root"))
        if(ntuple.Contains("nominal"))
        // ||ntuple.Contains("Rho")||ntuple.Contains("EffectiveNP_3")||ntuple.Contains("EffectiveNP_6"))
        {
          syst_name = ntuple;
          syst_name = string(syst_name).substr(0,string(syst_name).length()-18);

          // output = TFile::Open(output_path + region[k] + "/" + process_list[i]+".root","update");
          tBDT = new TTree(syst_name, syst_name);
          tBDT->SetDirectory(0);

          tBDT->Branch(TString(region[k]) + "BDT",&bdt_eval_s);
          tBDT->Branch("weight",&weight_b);
          // tBDT->Branch("dsid",&dsid_s);

          if(ntuple.Contains("nominal"))
          {
            // down
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1down_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1down_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1down_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1down_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1down_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1down_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1down_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1down_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1down_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1down",&bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1down_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1down",&bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1down_f);

            tBDT->Branch("lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_DOWN_AT",&lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_DOWN_AT",&lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_DOWN_AT",&lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_DOWN_AT",&lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_El_Reco_DOWN_AT",&lepSF_SF_El_Reco_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_El_ID_LooseAndBLayerLH_DOWN_AT",&lepSF_SF_El_ID_LooseAndBLayerLH_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_El_ID_TightLH_DOWN_AT",&lepSF_SF_El_ID_TightLH_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_El_Iso_FCLoose_DOWN_AT",&lepSF_SF_El_Iso_FCLoose_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_El_PLVLoose_DOWN_AT",&lepSF_SF_El_PLVLoose_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_El_PLVTight_DOWN_AT",&lepSF_SF_El_PLVTight_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_TTVA_STAT_DOWN_AT",&lepSF_SF_Mu_TTVA_STAT_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_TTVA_SYST_DOWN_AT",&lepSF_SF_Mu_TTVA_SYST_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Loose_STAT_DOWN_AT",&lepSF_SF_Mu_ID_Loose_STAT_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Medium_STAT_DOWN_AT",&lepSF_SF_Mu_ID_Medium_STAT_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Loose_SYST_DOWN_AT",&lepSF_SF_Mu_ID_Loose_SYST_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Medium_SYST_DOWN_AT",&lepSF_SF_Mu_ID_Medium_SYST_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Loose_STAT_LOWPT_DOWN_AT",&lepSF_SF_Mu_ID_Loose_STAT_LOWPT_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Medium_STAT_LOWPT_DOWN_AT",&lepSF_SF_Mu_ID_Medium_STAT_LOWPT_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Loose_SYST_LOWPT_DOWN_AT",&lepSF_SF_Mu_ID_Loose_SYST_LOWPT_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Medium_SYST_LOWPT_DOWN_AT",&lepSF_SF_Mu_ID_Medium_SYST_LOWPT_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_Iso_FCLoose_SYST_DOWN_AT",&lepSF_SF_Mu_Iso_FCLoose_SYST_DOWN_AT_f);
            // tBDT->Branch("lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_DOWN_AT",&lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_PLVLoose_DOWN_AT",&lepSF_SF_Mu_PLVLoose_DOWN_AT_f);
            tBDT->Branch("lepSF_SF_Mu_PLVTight_DOWN_AT",&lepSF_SF_Mu_PLVTight_DOWN_AT_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1down_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1down_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1down_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1down_f);
            
            // up
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1up_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1up_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1up_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1up_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1up_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1up_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1up_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1up_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1up_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1up",&bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1up_f);
            tBDT->Branch("bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1up",&bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1up_f);

            tBDT->Branch("lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_UP_AT",&lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_UP_AT_f);
            tBDT->Branch("lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_UP_AT",&lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_UP_AT_f);
            tBDT->Branch("lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_UP_AT",&lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_UP_AT_f);
            tBDT->Branch("lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_UP_AT",&lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_UP_AT_f);
            tBDT->Branch("lepSF_SF_El_Reco_UP_AT",&lepSF_SF_El_Reco_UP_AT_f);
            tBDT->Branch("lepSF_SF_El_ID_LooseAndBLayerLH_UP_AT",&lepSF_SF_El_ID_LooseAndBLayerLH_UP_AT_f);
            tBDT->Branch("lepSF_SF_El_ID_TightLH_UP_AT",&lepSF_SF_El_ID_TightLH_UP_AT_f);
            tBDT->Branch("lepSF_SF_El_Iso_FCLoose_UP_AT",&lepSF_SF_El_Iso_FCLoose_UP_AT_f);
            tBDT->Branch("lepSF_SF_El_PLVLoose_UP_AT",&lepSF_SF_El_PLVLoose_UP_AT_f);
            tBDT->Branch("lepSF_SF_El_PLVTight_UP_AT",&lepSF_SF_El_PLVTight_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_TTVA_STAT_UP_AT",&lepSF_SF_Mu_TTVA_STAT_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_TTVA_SYST_UP_AT",&lepSF_SF_Mu_TTVA_SYST_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Loose_STAT_UP_AT",&lepSF_SF_Mu_ID_Loose_STAT_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Medium_STAT_UP_AT",&lepSF_SF_Mu_ID_Medium_STAT_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Loose_SYST_UP_AT",&lepSF_SF_Mu_ID_Loose_SYST_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Medium_SYST_UP_AT",&lepSF_SF_Mu_ID_Medium_SYST_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Loose_STAT_LOWPT_UP_AT",&lepSF_SF_Mu_ID_Loose_STAT_LOWPT_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Medium_STAT_LOWPT_UP_AT",&lepSF_SF_Mu_ID_Medium_STAT_LOWPT_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Loose_SYST_LOWPT_UP_AT",&lepSF_SF_Mu_ID_Loose_SYST_LOWPT_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_ID_Medium_SYST_LOWPT_UP_AT",&lepSF_SF_Mu_ID_Medium_SYST_LOWPT_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_Iso_FCLoose_SYST_UP_AT",&lepSF_SF_Mu_Iso_FCLoose_SYST_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_UP_AT",&lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_UP_AT_f);
            // tBDT->Branch("lepSF_SF_Mu_PLVLoose_UP_AT",&lepSF_SF_Mu_PLVLoose_UP_AT_f);
            tBDT->Branch("lepSF_SF_Mu_PLVTight_UP_AT",&lepSF_SF_Mu_PLVTight_UP_AT_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1up_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1up_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1up_f);
            tBDT->Branch("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1up_f);
          }

          f = new TFile(file_path + ntuple,"read" );

          cout<<file_path + ntuple<<endl;
          tree = (TTree*)f->Get(syst_name);
          // file_path = region[k];

          // reader_s->BookMVA(TString(region[k]) + process_list[i] + syst_name + "_BDT method", "dataset/weights/nominal_BDT.weights.xml" );

          histBDT = new TH1F( syst_name + "_" + TString(region[k]) + "BDT",process_list[i],nbin,range[0],range[1]);
          // histBDT = new TH1F( syst_name + "_" + TString(region[k]) + "BDT",process_list[i],16,bin_edges);

          tree->SetBranchAddress( "lep_Pt_0", &lep_Pt_0_s);
          tree->SetBranchAddress( "lep_Pt_1", &lep_Pt_1_s);
          tree->SetBranchAddress( "lep_Pt_2", &lep_Pt_2_s);
          tree->SetBranchAddress( "lep_Pt_3", &lep_Pt_3_s);
          tree->SetBranchAddress( "lep_Etcone30_0", &lep_Etcone30_0_s);
          tree->SetBranchAddress( "lep_Etcone30_1", &lep_Etcone30_1_s);
          tree->SetBranchAddress( "lep_Etcone30_2", &lep_Etcone30_2_s);
          tree->SetBranchAddress( "lep_Etcone30_3", &lep_Etcone30_3_s);
          tree->SetBranchAddress( "lep_Eta_0", &lep_Eta_0_s);
          tree->SetBranchAddress( "lep_Eta_1", &lep_Eta_1_s);
          tree->SetBranchAddress( "lep_Eta_2", &lep_Eta_2_s);
          tree->SetBranchAddress( "lep_Eta_3", &lep_Eta_3_s);
          tree->SetBranchAddress( "jet_E_0", &jet_E_0_s);
          tree->SetBranchAddress( "jet_E_1", &jet_E_1_s);
          tree->SetBranchAddress( "jet_Pt_0", &jet_Pt_0_s);
          tree->SetBranchAddress( "jet_Pt_1", &jet_Pt_1_s);
          tree->SetBranchAddress( "jet_Phi_0", &jet_Phi_0_s);
          tree->SetBranchAddress( "jet_Phi_1", &jet_Phi_1_s);
          tree->SetBranchAddress( "m_12", &m_12_s);
          tree->SetBranchAddress( "m_34", &m_34_s);
          tree->SetBranchAddress( "m_4l", &m_4l_s);
          tree->SetBranchAddress( "m_jj", &m_jj_s);
          // tree->SetBranchAddress( "p_4l", &p_4l_s);
          tree->SetBranchAddress( "p_jj", &p_jj_s);
          tree->SetBranchAddress( "HT", &HT_s);
          // tree->SetBranchAddress( "HT_lep", &HT_lep_s);
          // tree->SetBranchAddress( "cs_pairs", &cs_pairs_s);
          //tree->SetBranchAddress( "cs_Z_pair", &cs_Z_pair_s);
          // tree->SetBranchAddress( "cs_lep_12", &cs_lep_12_s);
          // tree->SetBranchAddress( "cs_lep_34", &cs_lep_34_s);
          tree->SetBranchAddress( "met_met", &met_met_s);
          tree->SetBranchAddress( "Dphi_met_jets", &Dphi_met_jets_s);
          tree->SetBranchAddress( "njets", &njets_s);
          tree->SetBranchAddress( "nbjets", &nbjets_s);
          tree->SetBranchAddress( "mcWeight", &weight_s);
          tree->SetBranchAddress( "dsid", &dsid);

          if(ntuple.Contains("nominal"))
          {
            // down
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1down);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1down);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1down);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1down);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1down);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1down);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1down);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1down);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1down",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1down);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1down",&bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1down);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1down",&bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1down);

            tree->SetBranchAddress("lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_DOWN_AT",&lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_DOWN_AT",&lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_DOWN_AT",&lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_DOWN_AT",&lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_El_Reco_DOWN_AT",&lepSF_SF_El_Reco_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_El_ID_LooseAndBLayerLH_DOWN_AT",&lepSF_SF_El_ID_LooseAndBLayerLH_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_El_ID_TightLH_DOWN_AT",&lepSF_SF_El_ID_TightLH_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_El_Iso_FCLoose_DOWN_AT",&lepSF_SF_El_Iso_FCLoose_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_El_PLVLoose_DOWN_AT",&lepSF_SF_El_PLVLoose_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_El_PLVTight_DOWN_AT",&lepSF_SF_El_PLVTight_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_TTVA_STAT_DOWN_AT",&lepSF_SF_Mu_TTVA_STAT_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_TTVA_SYST_DOWN_AT",&lepSF_SF_Mu_TTVA_SYST_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Loose_STAT_DOWN_AT",&lepSF_SF_Mu_ID_Loose_STAT_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Medium_STAT_DOWN_AT",&lepSF_SF_Mu_ID_Medium_STAT_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Loose_SYST_DOWN_AT",&lepSF_SF_Mu_ID_Loose_SYST_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Medium_SYST_DOWN_AT",&lepSF_SF_Mu_ID_Medium_SYST_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Loose_STAT_LOWPT_DOWN_AT",&lepSF_SF_Mu_ID_Loose_STAT_LOWPT_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Medium_STAT_LOWPT_DOWN_AT",&lepSF_SF_Mu_ID_Medium_STAT_LOWPT_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Loose_SYST_LOWPT_DOWN_AT",&lepSF_SF_Mu_ID_Loose_SYST_LOWPT_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Medium_SYST_LOWPT_DOWN_AT",&lepSF_SF_Mu_ID_Medium_SYST_LOWPT_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_Iso_FCLoose_SYST_DOWN_AT",&lepSF_SF_Mu_Iso_FCLoose_SYST_DOWN_AT);
            // tree->SetBranchAddress("lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_DOWN_AT",&lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_PLVLoose_DOWN_AT",&lepSF_SF_Mu_PLVLoose_DOWN_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_PLVTight_DOWN_AT",&lepSF_SF_Mu_PLVTight_DOWN_AT);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1down);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1down);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1down);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1down",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1down);
            
            // up
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1up);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1up);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1up);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1up);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1up);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1up);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1up);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1up);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1up",&bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1up);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1up",&bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1up);
            tree->SetBranchAddress("bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1up",&bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1up);

            tree->SetBranchAddress("lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_UP_AT",&lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_UP_AT);
            tree->SetBranchAddress("lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_UP_AT",&lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_UP_AT);
            tree->SetBranchAddress("lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_UP_AT",&lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_UP_AT);
            tree->SetBranchAddress("lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_UP_AT",&lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_UP_AT);
            tree->SetBranchAddress("lepSF_SF_El_Reco_UP_AT",&lepSF_SF_El_Reco_UP_AT);
            tree->SetBranchAddress("lepSF_SF_El_ID_LooseAndBLayerLH_UP_AT",&lepSF_SF_El_ID_LooseAndBLayerLH_UP_AT);
            tree->SetBranchAddress("lepSF_SF_El_ID_TightLH_UP_AT",&lepSF_SF_El_ID_TightLH_UP_AT);
            tree->SetBranchAddress("lepSF_SF_El_Iso_FCLoose_UP_AT",&lepSF_SF_El_Iso_FCLoose_UP_AT);
            tree->SetBranchAddress("lepSF_SF_El_PLVLoose_UP_AT",&lepSF_SF_El_PLVLoose_UP_AT);
            tree->SetBranchAddress("lepSF_SF_El_PLVTight_UP_AT",&lepSF_SF_El_PLVTight_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_TTVA_STAT_UP_AT",&lepSF_SF_Mu_TTVA_STAT_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_TTVA_SYST_UP_AT",&lepSF_SF_Mu_TTVA_SYST_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Loose_STAT_UP_AT",&lepSF_SF_Mu_ID_Loose_STAT_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Medium_STAT_UP_AT",&lepSF_SF_Mu_ID_Medium_STAT_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Loose_SYST_UP_AT",&lepSF_SF_Mu_ID_Loose_SYST_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Medium_SYST_UP_AT",&lepSF_SF_Mu_ID_Medium_SYST_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Loose_STAT_LOWPT_UP_AT",&lepSF_SF_Mu_ID_Loose_STAT_LOWPT_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Medium_STAT_LOWPT_UP_AT",&lepSF_SF_Mu_ID_Medium_STAT_LOWPT_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Loose_SYST_LOWPT_UP_AT",&lepSF_SF_Mu_ID_Loose_SYST_LOWPT_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_ID_Medium_SYST_LOWPT_UP_AT",&lepSF_SF_Mu_ID_Medium_SYST_LOWPT_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_Iso_FCLoose_SYST_UP_AT",&lepSF_SF_Mu_Iso_FCLoose_SYST_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_UP_AT",&lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_UP_AT);
            // tree->SetBranchAddress("lepSF_SF_Mu_PLVLoose_UP_AT",&lepSF_SF_Mu_PLVLoose_UP_AT);
            tree->SetBranchAddress("lepSF_SF_Mu_PLVTight_UP_AT",&lepSF_SF_Mu_PLVTight_UP_AT);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1up);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1up);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1up);
            tree->SetBranchAddress("custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1up",&custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1up);
          }
          
          // Output histograms for evaluated value
          // TH1F *histBDT[8];
          // for(i = 0; i < 8; i++)  histBDT[i] = new TH1F(histname,histname,nbin,range[0],range[1]);

          // TH1F *histBDT_data = new TH1F(histname,histname,nbin,range[0],range[1]);
          // histBDT = new TH1F(process_list[i],process_list[i],nbin,range[0],range[1]);

          // Float_t bdt_eval_s, bdt_eval_b, bdt_eval_d;
          // sum = 0;
          for(j = 0; j < tree->GetEntries(); j++)
          {
            tree->GetEntry(j);
            // cout<<j<<endl;

            if(ntuple.Contains("nominal"))
            {
              // down
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1down_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1down;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1down_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1down;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1down_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1down;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1down_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1down;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1down_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1down;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1down_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1down;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1down_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1down;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1down_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1down;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1down_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1down;
              bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1down_f=bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1down;
              bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1down_f=bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1down;

              lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_DOWN_AT_f=lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_DOWN_AT;
              lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_DOWN_AT_f=lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_DOWN_AT;
              lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_DOWN_AT_f=lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_DOWN_AT;
              lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_DOWN_AT_f=lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_DOWN_AT;
              lepSF_SF_El_Reco_DOWN_AT_f=lepSF_SF_El_Reco_DOWN_AT;
              lepSF_SF_El_ID_LooseAndBLayerLH_DOWN_AT_f=lepSF_SF_El_ID_LooseAndBLayerLH_DOWN_AT;
              lepSF_SF_El_ID_TightLH_DOWN_AT_f=lepSF_SF_El_ID_TightLH_DOWN_AT;
              lepSF_SF_El_Iso_FCLoose_DOWN_AT_f=lepSF_SF_El_Iso_FCLoose_DOWN_AT;
              lepSF_SF_El_PLVLoose_DOWN_AT_f=lepSF_SF_El_PLVLoose_DOWN_AT;
              lepSF_SF_El_PLVTight_DOWN_AT_f=lepSF_SF_El_PLVTight_DOWN_AT;
              lepSF_SF_Mu_TTVA_STAT_DOWN_AT_f=lepSF_SF_Mu_TTVA_STAT_DOWN_AT;
              lepSF_SF_Mu_TTVA_SYST_DOWN_AT_f=lepSF_SF_Mu_TTVA_SYST_DOWN_AT;
              lepSF_SF_Mu_ID_Loose_STAT_DOWN_AT_f=lepSF_SF_Mu_ID_Loose_STAT_DOWN_AT;
              lepSF_SF_Mu_ID_Medium_STAT_DOWN_AT_f=lepSF_SF_Mu_ID_Medium_STAT_DOWN_AT;
              lepSF_SF_Mu_ID_Loose_SYST_DOWN_AT_f=lepSF_SF_Mu_ID_Loose_SYST_DOWN_AT;
              lepSF_SF_Mu_ID_Medium_SYST_DOWN_AT_f=lepSF_SF_Mu_ID_Medium_SYST_DOWN_AT;
              lepSF_SF_Mu_ID_Loose_STAT_LOWPT_DOWN_AT_f=lepSF_SF_Mu_ID_Loose_STAT_LOWPT_DOWN_AT;
              lepSF_SF_Mu_ID_Medium_STAT_LOWPT_DOWN_AT_f=lepSF_SF_Mu_ID_Medium_STAT_LOWPT_DOWN_AT;
              lepSF_SF_Mu_ID_Loose_SYST_LOWPT_DOWN_AT_f=lepSF_SF_Mu_ID_Loose_SYST_LOWPT_DOWN_AT;
              lepSF_SF_Mu_ID_Medium_SYST_LOWPT_DOWN_AT_f=lepSF_SF_Mu_ID_Medium_SYST_LOWPT_DOWN_AT;
              lepSF_SF_Mu_Iso_FCLoose_SYST_DOWN_AT_f=lepSF_SF_Mu_Iso_FCLoose_SYST_DOWN_AT;
              // lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_down_AT_f=lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_down_AT;
              lepSF_SF_Mu_PLVLoose_DOWN_AT_f=lepSF_SF_Mu_PLVLoose_DOWN_AT;
              lepSF_SF_Mu_PLVTight_DOWN_AT_f=lepSF_SF_Mu_PLVTight_DOWN_AT;
              custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down_f=custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down;
              custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1down_f=custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1down;
              custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1down_f=custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1down;
              custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down_f=custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down;
              custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1down_f=custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1down;
              custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1down_f=custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1down;
              
              // up
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1up_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_0__1up;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1up_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_1__1up;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1up_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_2__1up;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1up_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_3__1up;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1up_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_4__1up;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1up_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_5__1up;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1up_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_6__1up;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1up_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_7__1up;
              bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1up_f=bTagSF_weight_DL1r_77_FT_EFF_Eigen_B_8__1up;
              bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1up_f=bTagSF_weight_DL1r_77_FT_EFF_extrapolation__1up;
              bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1up_f=bTagSF_weight_DL1r_77_FT_EFF_extrapolation_from_charm__1up;

              lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_UP_AT_f=lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_STAT_UP_AT;
              lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_UP_AT_f=lepSF_SF_El_ChargeMisID_TightLH_FCLoose_STAT_UP_AT;
              lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_UP_AT_f=lepSF_SF_El_ChargeMisID_LooseAndBLayerLH_FCLoose_SYST_UP_AT;
              lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_UP_AT_f=lepSF_SF_El_ChargeMisID_TightLH_FCLoose_SYST_UP_AT;
              lepSF_SF_El_Reco_UP_AT_f=lepSF_SF_El_Reco_UP_AT;
              lepSF_SF_El_ID_LooseAndBLayerLH_UP_AT_f=lepSF_SF_El_ID_LooseAndBLayerLH_UP_AT;
              lepSF_SF_El_ID_TightLH_UP_AT_f=lepSF_SF_El_ID_TightLH_UP_AT;
              lepSF_SF_El_Iso_FCLoose_UP_AT_f=lepSF_SF_El_Iso_FCLoose_UP_AT;
              lepSF_SF_El_PLVLoose_UP_AT_f=lepSF_SF_El_PLVLoose_UP_AT;
              lepSF_SF_El_PLVTight_UP_AT_f=lepSF_SF_El_PLVTight_UP_AT;
              lepSF_SF_Mu_TTVA_STAT_UP_AT_f=lepSF_SF_Mu_TTVA_STAT_UP_AT;
              lepSF_SF_Mu_TTVA_SYST_UP_AT_f=lepSF_SF_Mu_TTVA_SYST_UP_AT;
              lepSF_SF_Mu_ID_Loose_STAT_UP_AT_f=lepSF_SF_Mu_ID_Loose_STAT_UP_AT;
              lepSF_SF_Mu_ID_Medium_STAT_UP_AT_f=lepSF_SF_Mu_ID_Medium_STAT_UP_AT;
              lepSF_SF_Mu_ID_Loose_SYST_UP_AT_f=lepSF_SF_Mu_ID_Loose_SYST_UP_AT;
              lepSF_SF_Mu_ID_Medium_SYST_UP_AT_f=lepSF_SF_Mu_ID_Medium_SYST_UP_AT;
              lepSF_SF_Mu_ID_Loose_STAT_LOWPT_UP_AT_f=lepSF_SF_Mu_ID_Loose_STAT_LOWPT_UP_AT;
              lepSF_SF_Mu_ID_Medium_STAT_LOWPT_UP_AT_f=lepSF_SF_Mu_ID_Medium_STAT_LOWPT_UP_AT;
              lepSF_SF_Mu_ID_Loose_SYST_LOWPT_UP_AT_f=lepSF_SF_Mu_ID_Loose_SYST_LOWPT_UP_AT;
              lepSF_SF_Mu_ID_Medium_SYST_LOWPT_UP_AT_f=lepSF_SF_Mu_ID_Medium_SYST_LOWPT_UP_AT;
              lepSF_SF_Mu_Iso_FCLoose_SYST_UP_AT_f=lepSF_SF_Mu_Iso_FCLoose_SYST_UP_AT;
              lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_UP_AT_f=lepSF_SF_Mu_Iso_FCLoose_STAT_DOWN_ATSF_Mu_PLVLoose_UP_AT;
              // lepSF_SF_Mu_PLVLoose_UP_AT_f=lepSF_SF_Mu_PLVLoose_UP_AT;
              lepSF_SF_Mu_PLVTight_UP_AT_f=lepSF_SF_Mu_PLVTight_UP_AT;

              custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up_f=custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up;
              custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1up_f=custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigStatUncertainty__1up;
              custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1up_f=custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT_MUON_EFF_TrigSystUncertainty__1up;
              custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up_f=custTrigSF_TightElMediumMuID_FCLooseIso_DLT_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up;
              custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1up_f=custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigStatUncertainty__1up;
              custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1up_f=custTrigSF_TightElMediumMuID_FCLooseIso_DLT_MUON_EFF_TrigSystUncertainty__1up;
            }

            // cout<<Dphi_met_jets_s<<endl;
            bdt_eval_s = reader_s->EvaluateMVA("BDT method");

            if(dsid != 342285)
            {
              // if(abs(weight_s) > 0.5) continue;
              if(i != N)  weight_b = weight_s;
              else  weight_b = 1.0;
              histBDT->Fill(bdt_eval_s,weight_b);
              
              tBDT->Fill();
            }
            // else if(i == N-1)  weight_b = 10*weight_s;
            // else  histBDT->Fill(bdt_eval_s);
            // sum = sum + weight_s;

            // ZH reweighting
            else
            {
              if(bdt_eval_s > 0.3)
              {
                for(l = 0; l < 4; l++)
                { 
                  bdt_eval_s = ZH[l][0];
                  weight_b = weight_s/ZH_sum*ZH[l][1];

                  histBDT->Fill(bdt_eval_s,weight_b);
                  tBDT->Fill();
                }
              }
            }
          }
          
          tree->Delete("");
          f->Close();
          
          // if(!(ntuple.Contains("nominal")))  histBDT->Smooth(1);

          // sum_1 = 0;
          // sum_2 = 0;

          // for(j = 1; j < histBDT->GetSize(); j++)
          // {
          //   sum_1 = sum_1 + histBDT->GetBinContent(j);

          //   if(histBDT->GetBinContent(j) < 1e-6)
          //   {
          //     histBDT->SetBinContent(j,1e-6);
          //     histBDT->SetBinError(j,1e-6);
          //   }

          //   sum_2 = sum_2 + histBDT->GetBinContent(j);
          // }
          // cout<<sum_1<<","<<sum_2<<endl;
          
          output = TFile::Open(output_path + region[k] + "/" + process_list[i]+".root","update");
          // cout<<"a"<<endl;
          // histBDT->Write();

          // tBDT->Fill();
          tBDT->Write();
          output->Write();
          output->Close();
        }

        // for(i = 0; i < nbin; i++)
        // {
        //   for(j = 0; j < N+1; j++)
        //     if(histBDT[j]->GetBinContent(i) < 1e-6)
        //     {
        //       histBDT[j]->SetBinContent(i,1e-6);
        //       histBDT[j]->SetBinError(i,1e-6);
        //     }
        // }

        // file_dir = file_dir + region[5];
        // for(i = 0; i < N+1; i++)
        // {
        //   output = new TFile(file_dir + region[k] + "/" + process_list[i]+".root","recreate");
        //   histBDT->Write();
        //   // histBDT->Reset("ICESM");
        //   output->Close();
        // }
      }

    }
  }
  delete reader_s;
  // delete reader_b;
  // delete reader_d;
}
