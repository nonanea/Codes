#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

const TString plotStorePath = "Plots/";

// vector<int> colorList{ kGreen, kGreen, kBlue, kBlue,kBlue, kRed, kRed, kRed, kBlue, kGreen, kGreen, kGreen, kGreen+1, kGreen+2,kViolet, kViolet, kCyan, kRed, kYellow, kPink, kOrange, kTeal, kMagenta, kAzure, kTeal + 1, kPink + 1, kMagenta + 1, kAzure + 1, kCyan + 1, kYellow + 1};
vector<int> colorList{kRed, kBlue, kTeal, kGreen, kViolet, kMagenta, kAzure, kCyan, kPink, kYellow, kBlue + 1, kRed + 1, kGreen + 1, kTeal + 1, kPink + 1, kMagenta + 1, kViolet + 1, kAzure + 1, kCyan + 1, kYellow + 1};

double interpolation(vector<double> x, vector<double> y, double x0);

void curve(TString dataPath)
{
    gErrorIgnoreLevel = kError;

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(.03);

    TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,1000,800);
    //c1->SetTopMargin(0.18);
 
    //c1->SetFillColor(42);
    //c1->SetGrid();
    //c1->GetFrame()->SetFillColor(21);
    //c1->GetFrame()->SetBorderSize(12);

    double x = 0, y = 0, x_err = 0,y_err = 0,xrange[2]={99999,0},yrange[2]={99999,0};
    float x_max = -1e10, y_max = -1e10, x_min = 1e10, y_min = 1e10;
    double Q = 4.5;
    int n = 0, i = 0;
    TString plotname = dataPath;
    plotname.ReplaceAll("/","");

    void *dataDir = gSystem->OpenDirectory(dataPath);
    const char  *entry;
    TPaveText *pt;
    TText *text;
    TString txtFile;
    const char* type;
    TLegend *leg;

    TGraphErrors *gr[20];
    TString legend[20];

    while((entry = (char *)gSystem->GetDirEntry(dataDir)))
    {
        txtFile = entry;
        
        if(txtFile.EndsWith(".txt") || txtFile.EndsWith(".csv") )
        {
            vector<double>  x_list,x_err_list,y_list,y_err_list;
            
            TString txtFileFullPath = dataPath + txtFile;

            TTree *tree = new TTree("graph","graph");
            // tree->ReadFile(txtFileFullPath, "x/D:y/D", ',');
            tree->ReadFile(txtFileFullPath,"x/D:y/D:x_err/D:y_err/D",',');
            tree->SetBranchAddress("x",&x);
            tree->SetBranchAddress("y",&y);
            tree->SetBranchAddress("x_err",&x_err);
            tree->SetBranchAddress("y_err",&y_err);

            tree->GetEntry(0);
            double x0 = x;

            // for(n = (tree->GetEntries())*0.1; n < tree->GetEntries(); n++)
            for(n = 0; n < tree->GetEntries(); n++)
            // for(n = 0; n < (tree->GetEntries())*0.2; n++)
            {
                tree->GetEntry(n);
                // cout<<x<<","<<y<<","<<x_err<<","<<y_err<<endl;
                if(txtFile.Contains("B3"))
                {
                    x_list.push_back(x*1e-9);
                    y_list.push_back(-y);
                }
                else if(txtFile.Contains("B7"))
                {
                    x_list.push_back(x);
                    y_list.push_back(y*20);
                }
                else   
                {
                    x_list.push_back(x);
                    y_list.push_back(y);
                }
                // if(i == 0)  y_list.push_back(y*500/1850);
                // else    y_list.push_back(y*500/2741);
                x_err_list.push_back(x_err);
                y_err_list.push_back(y_err);
                // x_err_list.push_back(0);
                // y_err_list.push_back(0);
            }
            
            tree->Delete("");
            
            // sort(x_list.begin(),x_list.end());

            gr[i] = new TGraphErrors(n,&(x_list[0]),&(y_list[0]),&(x_err_list[0]),&(y_err_list[0]));
            
            if(dataPath.Contains("Charge_Bias"))
            {
                cout<<"Charge:"<<Q<<" fC"<<endl;
                cout<<txtFile<<",bias:"<<interpolation(y_list,x_list,Q)<<" V"<<endl;
                // gr[i]->SetPoint(n,interpolation(y_list,x_list,Q),Q);
            }

            //max = gr[i]->GetHistogram()->GetMaximum();
            if(x_max < TMath::MaxElement(gr[i]->GetN(),gr[i]->GetX()))  x_max = TMath::MaxElement(gr[i]->GetN(),gr[i]->GetX());
            if(y_max < TMath::MaxElement(gr[i]->GetN(),gr[i]->GetY()))  y_max = TMath::MaxElement(gr[i]->GetN(),gr[i]->GetY());
            if(x_min > TMath::MinElement(gr[i]->GetN(),gr[i]->GetX()))  x_min = TMath::MinElement(gr[i]->GetN(),gr[i]->GetX());
            if(y_min > TMath::MinElement(gr[i]->GetN(),gr[i]->GetY()))  y_min = TMath::MinElement(gr[i]->GetN(),gr[i]->GetY());
            //x_max = TMath::MaxElement(gr[i]->GetN(),gr[i]->GetX());
            //y_max = TMath::MaxElement(gr[i]->GetN(),gr[i]->GetY());
            //cout<<x_max<<","<<y_max<<endl;

            txtFile = txtFile.ReplaceAll(".txt","");
            txtFile = txtFile.ReplaceAll(".csv","");
            legend[i] = txtFile;


            i++;
            // if(i > 4)   break;
        }
    }

    n = i;

    if(dataPath.Contains("Charge_Bias"))    
    {
        plotname = "Charge_Bias/" + plotname;
        // leg = new TLegend(0.2,0.9,0.4,0.9-n*0.03);
        leg = new TLegend(0.65,0.9,0.85,0.9-n*0.03);
        // leg = new TLegend(0.2,0.78,0.4,0.68);
        // gr[0]->SetTitle(";Bias[V];Charge[a.u.]");
        gr[0]->SetTitle(";Bias [V];Charge [fC]");
    }
    else if(dataPath.Contains("TRes_Bias"))   
    {
        plotname = "TRes_Bias/" + plotname;
        // leg = new TLegend(0.64,0.78,0.84,0.68);
        leg = new TLegend(0.2,0.9,0.4,0.9-n*0.03);
        // leg = new TLegend(0.65,0.9,0.85,0.9-n*0.03);
        gr[0]->SetTitle(";Bias [V];Time Resolution [ps]");
    }

    else if(dataPath.Contains("TRes_Charge"))
    {
        plotname = "TRes_Charge/" + plotname;
        leg = new TLegend(0.65,0.9,0.85,0.9-n*0.03);
        gr[0]->SetTitle(";Charge [fC];Time Resolution [ps]");
    }
    
    else
    {
        plotname = "Others/" + plotname;
        leg = new TLegend(0.2,0.9,0.4,0.9-n*0.03);
        // leg = new TLegend(0.64,0.9,0.84,0.9-n*0.03);
        // gr[0]->SetTitle(";Set;Time  Resolution [ps]");
        // gr[0]->SetTitle(";Time [ns];Ampl [V]");
        gr[0]->SetTitle(";Bias [V];Noise [mV]");
        // gr[0]->SetTitle(";Time [s];Leakage [A]");
        // gr[0]->SetTitle(";Bias [V];Charge [fC] ");
        // gr[0]->SetTitle(";#sigma_{DUT} [ps];#sqrt{V}(#sigma_{DUT})(Ref) [ps]");
    }

    gr[0]->GetYaxis()->SetRangeUser(y_min-(y_max-y_min)/10, y_max+(y_max-y_min)/10);
    // gr[0]->GetYaxis()->SetRangeUser(20, y_max+(y_max-y_min)/10);
    // gr[0]->GetYaxis()->SetRangeUser(20, 35);
    gr[0]->GetXaxis()->SetLimits(x_min-(x_max-x_min)/10, x_max+(x_max-x_min)/10);
    // gr[0]->GetXaxis()->SetLimits(150, 230);
    gStyle->SetEndErrorSize(5);
    
    int color = 0;
    for(i = 0; i < n; i++)
    {   
        if(legend[i].Contains("W7")) color = 3;
        else if(legend[i].Contains("W8")) color = 1;
        else if(legend[i].Contains("W9")) color = 2;
        else if(legend[i].Contains("W10")) color = 0;
        else if(legend[i].Contains("W11")) color = 0;
        else if(legend[i].Contains("W5")) color = 5;
        else if(legend[i].Contains("HPK")) color = 6;
        else    color = i;

        gr[i]->SetMarkerColor(colorList[color]);
        gr[i]->SetLineColor(colorList[color]);
        // gr[i]->SetMarkerColor(colorList[i]);
        // gr[i]->SetLineColor(colorList[i]);
        gr[i]->SetMarkerStyle(20 + i);
        gr[i]->SetMarkerSize(1);
        // gr[i]->SetLineWidth(2);
        if(legend[i].Contains(",M"))    gr[i]->SetLineStyle(10);
        if(i == 0)  gr[i]->Draw("ALP");
        // else if(i == 4 || i ==5)   gr[i]->Draw("LP");
        else    gr[i]->Draw("LP");

        type = legend[i];
        leg->SetTextSize(0.03);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        // leg->SetTextColor(colorList[i]);
        leg->AddEntry(gr[i],type,"LP");
        leg->Draw("same");
        
        pt=new TPaveText(0.6,0.9,0.7,0.85,"NDC");
        pt->SetBorderSize(0);
        // pt->SetFillColor(0);
        pt->SetFillStyle(0);
        pt->SetTextAlign(12);
        pt->SetTextSize(0.07);
        
        text=pt->AddText("#it{HPK Type 3.1}");
        // text=pt->AddText("#it{USTC W5}");
        // pt->Draw();
    }

    c1->SaveAs(plotStorePath + plotname + ".png");
    c1->SaveAs(plotStorePath + plotname + ".svg");
    // c1->SaveAs(plotStorePath + plotname + ".c");
    c1->SaveAs(plotStorePath + plotname + ".pdf");
 
    c1->Clear();
}

// double interpolation(vector<double> y_list, vector<double> x_list, double x)
// {
//     int i, j, n = x_list.size();
//     double temp = 1.0, res = 0.0;
//     double N[20][20] = {0};

//     for(i = 0; i < n; i++)  N[0][i] = y_list[i];
//     for(i = 1; i < n; i++)
//     {
//         for(j = i; j < n; j++)  N[i][j] = (N[i-1][j]-N[i-1][j-1])/(x_list[j]-x_list[j-1]);
//     }

//     // res = N[0][0];
//     for(i = 0; i < n; i++)
//     {
//         temp = N[i][i];
//         for(j = 1; j < i; j++)  temp=temp*(x - x_list[j-1]);
//         res+=temp;
//     }

//     return res;
// }

double interpolation(vector<double> x, vector<double> y, double x0)
{
	int i, j, n = x.size();
	double l[20], res = 0.0;

	for(i = 0; i < n; i++) l[i] = 1;
	
    for(i = 0;i < n; i++)
    {
		for(j = 0; j < n; j++)
        {
			if( j == i)
				l[i] = l[i];
			else
				l[i] = l[i]*(x0-x[j])/(x[i]-x[j]);
		}
	}
	
    for(i = 0; i < n; i++)   res = res + l[i]*y[i];
	
    return res;
}