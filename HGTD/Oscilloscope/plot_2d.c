#include<fstream>
#include <iostream>
#include <vector>

using namespace std;
const TString plotPath = "Plots/";

void plot_2d(TString dataPath)
{   
    gErrorIgnoreLevel = kError;

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(.03);

    TCanvas *c1 = new TCanvas("c1","A Simple 2D Graph",200,100,750,500);
    c1->SetRightMargin(0.32);


    TString txtFile;
    const char  *entry;
    void *dataDir = gSystem->OpenDirectory(dataPath);

    // double x[5] = {0};
    // double y[5] = {0};
    // double z;

    int i, j, n;
    double x,y,z,A,x0,y0,xrange[2]={50,300},yrange[2]={0,50},xbins,ybins;
    
    TH2 *h_2d;

    xbins = 150;
    ybins = 50;
    // h_2d = new TH2D("plotname","plotname",xbins,xrange[0],xrange[1],ybins,yrange[0],yrange[1]);
    h_2d = new TH2D("plotname","plotname",200,0,500,10,0,10*2.5*15);
    
    while((entry = (char *)gSystem->GetDirEntry(dataDir)))
    {
        txtFile = entry;
        
        if(txtFile.EndsWith(".txt") || txtFile.EndsWith(".csv") )
        {
            cout<<txtFile<<endl;
            vector<double>  x_list,y_list,z_list;
            
            TString txtFileFullPath = dataPath + txtFile;

            TTree *tree = new TTree("graph","graph");
            // tree->ReadFile(txtFileFullPath, "x/D:y/D", ',');
            tree->ReadFile(txtFileFullPath,"x/D:y/D:z/D:A/D",',');
            tree->SetBranchAddress("x",&x);
            tree->SetBranchAddress("y",&y);
            tree->SetBranchAddress("A",&z);

            tree->GetEntry(0);

            x0 = x;
            y0 = y;

            for(n = 0; n < tree->GetEntries(); n++)
            {
                tree->GetEntry(n);
                
                // if(x != x0 || y != y0)
                if(x != xbins)
                {
                x_list.push_back(x-x0);
                y_list.push_back((y-y0)/15);
                z_list.push_back(z*1000+15);
                // z_list.push_back(-z*1000);
                // cout<<x<<","<<y<<","<<-z<<endl;

                // x0 = x;
                xbins = x;
                }
            }

            // h_2d = new TH2D("plotname","plotname",x_list.size(),xrange[0],xrange[1],y_list.size(),yrange[0],yrange[1]);

            for(i = 0; i < x_list.size(); i++)
            {
                // for(j = 0; j < sqrt(y_list.size()); j++)
                // {
                    if(i%10 == 0 ) cout<<x_list[i]<<","<<y_list[i]<<","<<z_list[i]<<endl;
                    // if(z_list[i] > 37)  
                    // A = 0.6;
                    // else    A = 0.1;
                    h_2d->SetBinContent(x_list[i],y_list[i],z_list[i]);
                    // h_2d->SetBinContent(int(x_list[i]/4),y_list[i],A);
                    // if(z_list[i*sqrt(y_list.size())+j] > 2.5) h_2d->SetBinContent(j+1,i+1,z_list[i*sqrt(y_list.size())+j]);
                    // if(z_list[i*sqrt(y_list.size())+j] > 2.5) h_2d->SetBinContent(i+1,j+1,1);
                    // else h_2d->SetBinContent(i+1,j+1,0);
                    // if(i/100 == 0 )  cout<<z_list[i]<<endl;
                // }
            }
            // cout<<i<<endl;
        }
    }

    // Int_t MyPalette[100];
    // Double_t r[]    = {0., 0.0, 1.0, 1.0, 1.0};
    // Double_t g[]    = {0., 0.0, 0.0, 1.0, 1.0};
    // Double_t b[]    = {0., 1.0, 0.0, 0.0, 1.0};
    // Double_t stop[] = {0., .25, .50, .75, 1.0};
    // Int_t FI = TColor::CreateGradientColorTable(5, stop, r, g, b, 50);
    // for (int i=0;i<50;i++) MyPalette[i] = FI+i;
    
    // gStyle->SetOptStat();

    h_2d->GetXaxis()->SetNdivisions(5, kTRUE);
    h_2d->GetYaxis()->SetNdivisions(5, kTRUE);
    
    // h_2d->GetXaxis()->SetRangeUser(0, 300);
    // h_2d->GetYaxis()->SetRangeUser(0.5, 4);
    // h_2d->GetZaxis()->SetRangeUser(0, 50);
    // h_2d->GetZaxis()->SetRangeUser(0, 1);
    // h_2d->SetContour(99);

    // h_2d->GetXaxis()->SetTitle("Bias [V]");
    // h_2d->GetYaxis()->SetTitle("Charge [fC]");
    // h_2d->GetZaxis()->SetTitle("Time Resolution [ps]");

    h_2d->GetXaxis()->SetTitle("X [um]");
    h_2d->GetYaxis()->SetTitle("Y [um]");
    h_2d->GetZaxis()->SetTitle("Ampl [a.u.]");
    
    h_2d->Draw("COLZ");
    // h_2d->Draw("CONT");

    c1->SaveAs(plotPath + txtFile.ReplaceAll(".txt","") + "_2d.png");
    c1->SaveAs(plotPath + txtFile.ReplaceAll(".txt","") + "_2d.svg");
    c1->SaveAs(plotPath + txtFile.ReplaceAll(".txt","") + "_2d.pdf");
    // c1->SaveAs(plotPath + "VCT_2d.png");
    // c1->SaveAs(plotPath + "VCT_2d.svg");
    // c1->SaveAs(plotPath + "VCT_2d.pdf");
}
