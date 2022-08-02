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
const int N = 7;

void hist()
{
    gErrorIgnoreLevel = kError;

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(.03);

    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 1000, 800);
    c1->SetLeftMargin(0.25);
    c1->SetRightMargin(0.2);
    c1->SetTopMargin(0.2);

    TString plotPath = "Plots/";

    const char *plotname[10] = {"m_jj",
                                "pt_jj",
                                "m_4l",
                                "pt_4l",
                                "m_2l_leading",
                                "m_2l_subleading",
                                "nJets",
                                "nbJets"
                                };
    const char *xaxis_name[10] = {"M_{jj} [GeV]",
                                    "P_{T,jj} [GeV]",
                                    "M_{4l} [GeV]",
                                    "P_{T,4l} [GeV]",
                                    "M_{leading_pair} [GeV]",
                                    "M_{subleading_pair} [GeV]",
                                    "nJets",
                                    "nbJets"
                                    };
    float xrange[10][2] = {{0,800},
                            {0,800},
                            {100,140},
                            {0,800},
                            {0,200},
                            {0,70},
                            {2,10},
                            {1,6}
                            };

    int i,j,k;
    Float_t m_12_s, m_34_s, p_4l_s, p_jj_s, m_4l_s, m_jj_s;
    Float_t m_12_b, m_34_b, p_4l_b, p_jj_b, m_4l_b, m_jj_b;
    Float_t njets_s, nbjets_s, njets_b, nbjets_b;
    Float_t weight_s, weight_b;

    TH1D *hList[10][3];
    THStack *hs[10];
    TTree *t_s;
    TFile *f;
    
    // t_b->SetBranchAddress("m_12",&m_12_b);
    // t_b->SetBranchAddress("m_34",&m_34_b);
    // t_b->SetBranchAddress("m_4l",&m_4l_b);
    // t_b->SetBranchAddress("m_jj",&m_jj_b);
    // t_b->SetBranchAddress("p_4l",&p_4l_b);
    // t_b->SetBranchAddress("p_jj",&p_jj_b);
    // t_b->SetBranchAddress("njets",&njets_b);
    // t_b->SetBranchAddress("nbjets",&nbjets_b);
    // t_b->SetBranchAddress("weight",&weight_b);

    // const char* f_name[3] = {"CategoryReduction_JET_JER_EffectiveNP_6__1down","CategoryReduction_JET_JER_EffectiveNP_6__1up","nominal"};
    const char* process[N] = {"tt","ttV","VV","Higgs","Z+jets","signal","data"}

    // TFile *f = new TFile("TMVA_output.root","read");
    for(i = 0; i < N; i++)
    {
        hList[0][i] = new TH1D("m_jj", "m_jj", 20, 0, 800);
        hList[1][i] = new TH1D("p_jj", "p_jj", 20, 0, 800);
        hList[2][i] = new TH1D("m_4l", "m_4l", 20, 100, 140);
        hList[3][i] = new TH1D("p_4l", "p_4l", 20, 0, 800);
        hList[4][i] = new TH1D("m_2l_leading", "m_2l_leading", 20, 0, 200);
        hList[5][i] = new TH1D("m_2l_subleading", "m_2l_subleaidng", 20, 0, 70);
        hList[6][i] = new TH1D("njets", "nbjets", 8, 2, 10);
        hList[7][i] = new TH1D("nbjets", "njets", 5, 1, 6);
    }

    float scale;
    TH1* frame;
    TPaveText *pt;
    TText *text;
    TLegend *leg;

    gStyle->SetOptStat("");
    for(i = 0; i < 8; i++)
    {   
        leg = new TLegend(0.45,0.75,0.30,0.68);
        // leg->SetFillColor(0);
        // leg->SetLineColor(0);
        leg->SetTextSize(0.02);
        for(k = 0; k < N; k++)
        {
            cout<<TString(process[k])+"_Preselection.root"<<endl;
            f = new TFile(TString("Data/"+process[k])+"/nominal_Preselection.root","read");
            t_s = (TTree*)f->Get("nominal");

            t_s->SetBranchAddress("m_12",&m_12_s);
            t_s->SetBranchAddress("m_34",&m_34_s);
            t_s->SetBranchAddress("m_4l",&m_4l_s);
            t_s->SetBranchAddress("m_jj",&m_jj_s);
            t_s->SetBranchAddress("p_4l",&p_4l_s);
            t_s->SetBranchAddress("p_jj",&p_jj_s);
            t_s->SetBranchAddress("njets",&njets_s);
            t_s->SetBranchAddress("nbjets",&nbjets_s);
            t_s->SetBranchAddress("mcWeight",&weight_s);

            for(j = 0; j < t_s->GetEntries(); j++)
            {
                t_s->GetEntry(j);
                double variables_s[10] = {m_jj_s,p_jj_s,m_4l_s,p_4l_s,m_12_s/1000,m_34_s/1000,njets_s,nbjets_s};

                hList[i][k]->Fill(variables_s[i],weight_s);
            }
            t_s->Delete("");

            if(k == 0)
            {
                scale = (hList[i][k]->GetMaximum())*3.0/2;
                frame = gPad->DrawFrame(xrange[i][0],0,xrange[i][1],scale);
                gPad->Update();
                frame->GetXaxis()->SetNdivisions(505, kTRUE);
                frame->GetYaxis()->SetTitle("Events");
                frame->GetXaxis()->SetTitle(xaxis_name[i]);
                
                // pt=new TPaveText(0.38,0.78,0.28,0.58,"NDC");
                // pt->SetBorderSize(0);
                // pt->SetFillColor(0);
                // pt->SetFillStyle(0);
                // pt->SetTextAlign(12);
                // pt->SetTextSize(0.03);

                // text=pt->AddText("#scale[1.3]{#it{ATLAS} internal}");
                // text=pt->AddText("HH#rightarrow b#bar{b}ZZ#rightarrow b#bar{b}4l");
                // text=pt->AddText("13 TeV, 58.5 fb^{-1}");

                // pt->Draw("same")；
            }

            // hList[i][k]->SetLineColor(k+1);
            // if(TString(f_name[k]) == "nominal")  hList[i][k]->SetLineStyle(kDashed);
            // leg->AddEntry(hList[i][k],f_name[k],"f");
            // leg->Draw("same");

            // hList[i][k]->Draw("same hist");
        }

        rp = new TRatioPlot(hList[i][k]);
        rp->Draw();

        plotPath = plotPath + plotname[i];

        c1->SaveAs(plotPath+"_dist.png");
        c1->SaveAs(plotPath+"_dist.svg");

        c1->Clear();
    }
}
