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

    TCanvas *c1 = new TCanvas("c1","A Simple 2D Graph",200,100,600,500);
    c1->SetRightMargin(0.18);


    TString txtFile;
    const char  *entry;
    void *dataDir = gSystem->OpenDirectory(dataPath);

    // double x[5] = {0};
    // double y[5] = {0};
    // double z;

    int i, j, n;
    double x,y,z,A,x0,y0,xrange[2]={0,6},yrange[2]={0,6};
    
    TH2D *h_2d;

    while((entry = (char *)gSystem->GetDirEntry(dataDir)))
    {
        txtFile = entry;
        
        if(txtFile.EndsWith(".txt") || txtFile.EndsWith(".csv") )
        {
            vector<double>  x_list,y_list,z_list;
            
            TString txtFileFullPath = dataPath + txtFile;

            TTree *tree = new TTree("graph","graph");
            // tree->ReadFile(txtFileFullPath, "x/D:y/D", ',');
            tree->ReadFile(txtFileFullPath,"x/D:y/D:z/D:A/D",',');
            tree->SetBranchAddress("x",&x);
            tree->SetBranchAddress("y",&y);
            tree->SetBranchAddress("A",&z);

            x0 = -1;
            y0 = -1;

            for(n = 0; n < tree->GetEntries(); n++)
            {
                tree->GetEntry(n);
                
                if(x != x0 || y != y0)
                {
                x_list.push_back(x);
                y_list.push_back(y);
                z_list.push_back(-z*1000);
                // cout<<x<<","<<y<<","<<-z<<endl;

                x0 = x;
                y0 = y;
                }
            }
            h_2d = new TH2D("plotname","plotname",sqrt(x_list.size()),xrange[0],xrange[1],sqrt(y_list.size()),yrange[0],yrange[1]);

            for(i = 0; i < sqrt(x_list.size()); i++)
            {
                for(j = 0; j < sqrt(y_list.size()); j++)
                {
                    // cout<<x_list[j]<<","<<y_list[i]<<","<<z_list[i*sqrt(y_list.size())+j]<<endl;
                    // h_2d->SetBinContent(x_list[i],y_list[j],z_list[i*sqrt(y_list.size())+j]);
                    if(z_list[i*sqrt(y_list.size())+j] > 2.5) h_2d->SetBinContent(j+1,i+1,z_list[i*sqrt(y_list.size())+j]);
                    // if(z_list[i*sqrt(y_list.size())+j] > 2.5) h_2d->SetBinContent(i+1,j+1,1);
                    else h_2d->SetBinContent(i+1,j+1,0);
                    // if(i == 200 && j < 50 )  cout<<j<<","<<z_list[i*sqrt(y_list.size())+j]<<",0,0"<<endl;
                }
            }
        }
    }

    gStyle->SetOptStat();

    // h_2d->GetXaxis()->SetNdivisions(4, kTRUE);
    // h_2d->GetZaxis()->SetNdivisions(1, kTRUE);
    
    // h_2d->GetXaxis()->SetRangeUser(1.5, 5);
    // h_2d->GetYaxis()->SetRangeUser(0.5, 4);
    // h_2d->GetYaxis()->SetRangeUser(0.5, 4);

    h_2d->GetXaxis()->SetTitle("X [mm]");
    h_2d->GetYaxis()->SetTitle("Y [mm]");
    h_2d->GetZaxis()->SetTitle("Ampl [mV]");
    
    h_2d->Draw("COLZ");

    c1->SaveAs(plotPath + txtFile.ReplaceAll(".txt","") + "_2d.png");
    c1->SaveAs(plotPath + txtFile.ReplaceAll(".txt","") + "_2d.svg");
    c1->SaveAs(plotPath + txtFile.ReplaceAll(".txt","") + "_2d.pdf");
}