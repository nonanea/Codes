#include<math.h>

const int N = 6;

void input(TString treename)
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
    Float_t lep_Pt_0_s, lep_Pt_1_s, lep_Pt_2_s, lep_Pt_3_s, lep_Phi_0_s, lep_Phi_1_s, lep_Phi_2_s, lep_Phi_3_s, jet_E_0_s, jet_E_1_s, jet_Pt_0_s, jet_Pt_1_s, jet_Phi_0_s, jet_Phi_1_s, m_12_s, m_34_s, p_4l_s, p_jj_s, m_4l_s, m_jj_s, cs_pairs_s, cs_Z_pair_s, cs_lep_12_s, cs_lep_34_s, met_met_s, HT_s, HT_lep_s, Dphi_met_jets_s, njets_s, nbjets_s;
    Float_t lep_Pt_0_b, lep_Pt_1_b, lep_Pt_2_b, lep_Pt_3_b, lep_Phi_0_b, lep_Phi_1_b, lep_Phi_2_b, lep_Phi_3_b, jet_E_0_b, jet_E_1_b, jet_Pt_0_b, jet_Pt_1_b, jet_Phi_0_b, jet_Phi_1_b, m_12_b, m_34_b, p_4l_b, p_jj_b, m_4l_b, m_jj_b, cs_pairs_b, cs_Z_pair_b, cs_lep_12_b, cs_lep_34_b, met_met_b, HT_b, HT_lep_b, Dphi_met_jets_b, njets_b, nbjets_b;
    Float_t lep_Pt_0_d, lep_Pt_1_d, lep_Pt_2_d, lep_Pt_3_d, lep_Phi_0_d, lep_Phi_1_d, lep_Phi_2_d, lep_Phi_3_d, jet_E_0_d, jet_E_1_d, jet_Pt_0_d, jet_Pt_1_d, jet_Phi_0_d, jet_Phi_1_d, m_12_d, m_34_d, p_4l_d, p_jj_d, m_4l_d, m_jj_d, cs_pairs_d, cs_Z_pair_d, cs_lep_12_d, cs_lep_34_d, met_met_d, HT_d, HT_lep_d, Dphi_met_jets_d, njets_d, nbjets_d;
    Float_t weight_s, weight_b;

    reader_s->AddVariable( "lep_Pt_0:= lep_Pt_0", &lep_Pt_0_s);
    reader_s->AddVariable( "lep_Pt_1:= lep_Pt_1", &lep_Pt_1_s);
    reader_s->AddVariable( "lep_Pt_2:= lep_Pt_2", &lep_Pt_2_s);
    reader_s->AddVariable( "lep_Pt_3:= lep_Pt_3", &lep_Pt_3_s);
    // reader_s->AddVariable( "lep_Phi_0:= lep_Phi_0", &lep_Phi_0_s);
    // reader_s->AddVariable( "lep_Phi_1:= lep_Phi_1", &lep_Phi_1_s);
    // reader_s->AddVariable( "lep_Phi_2:= lep_Phi_2", &lep_Phi_2_s);
    // reader_s->AddVariable( "lep_Phi_3:= lep_Phi_3", &lep_Phi_3_s);
    reader_s->AddVariable( "jet_E_0:= jet_E_0", &jet_E_0_s);
    reader_s->AddVariable( "jet_E_1:= jet_E_1", &jet_E_1_s);
    reader_s->AddVariable( "jet_Pt_0:= jet_Pt_0", &jet_Pt_0_s);
    reader_s->AddVariable( "jet_Pt_1:= jet_Pt_1", &jet_Pt_1_s);
    reader_s->AddVariable( "jet_Phi_0:= jet_Phi_0", &jet_Phi_0_s);
    reader_s->AddVariable( "jet_Phi_1:= jet_Phi_1", &jet_Phi_1_s);
    reader_s->AddVariable( "m_12:= m_12", &m_12_s);
    reader_s->AddVariable( "m_34:= m_34", &m_34_s);
    reader_s->AddVariable( "m_4l:= m_4l", &m_4l_s);
    reader_s->AddVariable( "m_jj:= m_jj", &m_jj_s);
    // reader_s->AddVariable( "p_4l:= p_4l", &p_4l_s);
    reader_s->AddVariable( "p_jj:= p_jj", &p_jj_s);
    // reader_s->AddVariable( "HT:= HT", &HT_s);
    // reader_s->AddVariable( "HT_lep:= HT_lep", &HT_lep_s);
    reader_s->AddVariable( "cs_lep_12:= cs_lep_12", &cs_lep_12_s);
    reader_s->AddVariable( "cs_lep_34:= cs_lep_34", &cs_lep_34_s);
    reader_s->AddVariable( "cs_pairs:= cs_pairs", &cs_pairs_s);
    reader_s->AddVariable( "met_met:= met_met", &met_met_s);
    reader_s->AddVariable( "Dphi_met_jets:= Dphi_met_jets", &Dphi_met_jets_s);
    //reader_s->AddVariable( "cs_Z_pair:= cs_Z_pair", &cs_Z_pair_s);
    reader_s->AddVariable( "njets:= njets", &njets_s);
    reader_s->AddVariable( "nbjets:= nbjets", &nbjets_s);

    // We are going to evaluate only with MLP method
    reader_s->BookMVA( "BDT method", "dataset/weights/TMVAClassification_BDT.weights.xml" );
    // reader_b->BookMVA( "BDT method", "dataset/weights/TMVAClassification_BDT.weights.xml" );
    // reader_d->BookMVA( "BDT method", "dataset/weights/TMVAClassification_BDT.weights.xml" );

    Int_t i, j, k, bkg, nbin = 16;
    Double_t range[2] = {-0.6,0.4};

    TString file_dir;
    // const char* histname;
    const char *process_list[N+1] = {"tt","VV","ttV","Higgs","Z+jets","Signal","data"};
    const char *region[N] = {"tt","ttV","VV+Higgs","Z+jets","SR","VR"};

    // f_s = "Ntuples/signal/" + process_list[i] + "_Preselection.root";
    // f_b = "Ntuples/background/" + process_list[i] + "_Preselection.root";
    // histname = file.ReplaceAll("_Preselection.root","") + "BDT";

    // Prepare input tree for application
    // TFile *input = TFile::Open( filename,"read" );
    TFile *f, *output;
    TTree *tree;
    TH1F *histBDT;
    // for(i = 0; i < N+1; i++)  histBDT[i] = new TH1F("BDT",process_list[i],nbin,range[0],range[1]);

    Float_t bdt_eval_s;
    // double sum;

    // for(k = 0; k == 0; k++)
    for(k = 0; k < N; k++)
    {
      for(i = 0; i < N+1; i++)
      {
        file_dir = "Ntuples/";
        cout<<file_dir + region[k] + "/" + process_list[i] + "/" + treename + "_Preselection.root"<<endl;
        f = TFile::Open(file_dir + region[k] + "/" + process_list[i] + "/" + treename + "_Preselection.root","read" );

        tree = (TTree*)f->Get(treename);
        file_dir = region[k];
        histBDT = new TH1F(file_dir + "BDT",process_list[i],nbin,range[0],range[1]);

        tree->SetBranchAddress( "lep_Pt_0", &lep_Pt_0_s);
        tree->SetBranchAddress( "lep_Pt_1", &lep_Pt_1_s);
        tree->SetBranchAddress( "lep_Pt_2", &lep_Pt_2_s);
        tree->SetBranchAddress( "lep_Pt_3", &lep_Pt_3_s);
        // tree->SetBranchAddress( "lep_Phi_0", &lep_Phi_0_s);
        // tree->SetBranchAddress( "lep_Phi_1", &lep_Phi_1_s);
        // tree->SetBranchAddress( "lep_Phi_2", &lep_Phi_2_s);
        // tree->SetBranchAddress( "lep_Phi_3", &lep_Phi_3_s);
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
        tree->SetBranchAddress( "met_met", &met_met_s);
        // tree->SetBranchAddress( "HT", &HT_s);
        // tree->SetBranchAddress( "HT_lep", &HT_lep_s);
        tree->SetBranchAddress( "cs_pairs", &cs_pairs_s);
        //tree->SetBranchAddress( "cs_Z_pair", &cs_Z_pair_s);
        tree->SetBranchAddress( "cs_lep_12", &cs_lep_12_s);
        tree->SetBranchAddress( "cs_lep_34", &cs_lep_34_s);
        tree->SetBranchAddress( "Dphi_met_jets", &Dphi_met_jets_s);
        tree->SetBranchAddress( "njets", &njets_s);
        tree->SetBranchAddress( "nbjets", &nbjets_s);
        tree->SetBranchAddress( "mcWeight", &weight_s);

        // Output histograms for evaluated value
        // TH1F *histBDT[8];
        // for(i = 0; i < 8; i++)  histBDT[i] = new TH1F(histname,histname,nbin,range[0],range[1]);

        // TH1F *histBDT_data = new TH1F(histname,histname,nbin,range[0],range[1]);
        // histBDT = new TH1F(process_list[i],process_list[i],nbin,range[0],range[1]);

        // Float_t bdt_eval_s, bdt_eval_b, bdt_eval_d;
        // sum = 0;
        for (j = 0; j < tree->GetEntries(); j++) {
          tree->GetEntry(j);

          // cout<<Dphi_met_jets_s<<endl;
          // if(TMath::IsNaN(cs_pairs_s))  cout<<cs_pairs_s<<endl;
          bdt_eval_s = reader_s->EvaluateMVA("BDT method");
          // if(weight_s < 1e-6) weight_s = 1e-6;

          if(i != N)  histBDT->Fill(bdt_eval_s,weight_s);
          else  histBDT->Fill(bdt_eval_s);
          // sum = sum + weight_s;
        }
        tree->Delete("");
        
        for(j = 0; j < nbin; j++)
        {
          if(histBDT->GetBinContent(j) < 1e-6)
          {
            histBDT->SetBinContent(j,1e-6);
            histBDT->SetBinError(j,1e-6);
          }
        }

        file_dir = "Results/";
        output = new TFile(file_dir + region[k] + "/" + process_list[i]+".root","recreate");
        histBDT->Write();
        // histBDT->Reset("ICESM");
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

    delete reader_s;
    // delete reader_b;
    // delete reader_d;
}
