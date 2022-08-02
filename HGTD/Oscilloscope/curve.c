#include <fstream>
#include <iostream>
#include <vector>
#include "TH1.h"
#include "TF1.h"
#include <sstream>
#include <iomanip>

using namespace std;

const TString plotStorePath = "Plots/";
bool fit = false;
int row = 3;
int column = 4;

// vector<int> colorList{ kGreen, kGreen, kBlue, kBlue,kBlue, kRed, kRed, kRed, kBlue, kGreen, kGreen, kGreen, kGreen+1, kGreen+2,kViolet, kViolet, kCyan, kRed, kYellow, kPink, kOrange, kTeal, kMagenta, kAzure, kTeal + 1, kPink + 1, kMagenta + 1, kAzure + 1, kCyan + 1, kYellow + 1};
vector<int> colorList{kRed, kBlue, kGreen + 2, kViolet, kOrange + 4, kTeal, kMagenta, kAzure, kCyan, kYellow, kPink, kGreen, kBlue + 1, kRed + 1, kTeal + 1, kPink + 1, kMagenta + 1, kViolet + 1, kAzure + 1, kBlack, kCyan + 1};

double interpolation(vector<double> x, vector<double> y, double x0);

void curve(TString dataPath, bool fit = false)
{
    gErrorIgnoreLevel = kError;

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(.03);

    TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,1000,800);
    c1->SetTopMargin(0.05);
    c1->SetRightMargin(0.05);
 
    //c1->SetFillColor(42);
    //c1->SetGrid();
    //c1->GetFrame()->SetFillColor(21);
    //c1->GetFrame()->SetBorderSize(12);

    double x = 0, y = 0, x_err = 0, y_err = 0, xrange[2]={99999,0}, yrange[2]={99999,0}, a = 0, b = 0, c = 0, d = 0;
    double x0, y0;
    float x_max = -1e10, y_max = -1e10, x_min = 1e10, y_min = 1e10;
    double Q = 4;
    int n = 0, i = 0, j = 0;
    double k1[5] = {37.2,29,73.8,39.9,38.3}, k2[5] = {0.584,0.434,0.725,0.572,0.185};
    TString plotname = dataPath;
    plotname.ReplaceAll("/","");

    void *dataDir = gSystem->OpenDirectory(dataPath);
    const char  *entry;
    TPaveText *pt;
    TText *text;
    TString txtFile;
    const char* type;

    TF1 *f1, *f2, *f3, *f4, *f5, *f6;
    f1 = new TF1("f1","expo",0,1000);
    f1->SetParameter(0,a);
    f1->SetParameter(1,b);
    // f1->SetParameter(2,1);

    f2 = new TF1("f2","[2]-expo",0,10);
    // TF1 *f2 = new TF1("f2","[0]+[1]*x",0,10);
    f2->SetParameter(0,a);
    f2->SetParameter(1,b);
    f2->SetParameter(2,c);
    // f2->SetParameter(2,d);

    f3 = new TF1("f3","[0]+[1]*x",180,220);
    // f3 = new TF1("f3","[0]*x",0,100);
    // TF1 *f2 = new TF1("f2","[0]+[1]*x",0,10);
    // f3->SetParameter(0,a);
    // f3->SetParameter(1,b);
    f3->SetParameters(50,-0.001);

    // f4 = new TF1("f4","[3]*exp(-[2]/(1+(x-[0])*(x-[0])/[1]/[1]))",-9600,-8050);
    f4 = new TF1("f4","[0]*ROOT::Math::erfc((x-[1])*[2])+[3]",150,250);
    // TF1 *f2 = new TF1("f2","[0]+[1]*x",0,10);
    f4->SetParameters(-1.0,150.0,1,5);
    // f4->SetParameter(0,-8800);
    // f4->SetParameter(1,1);
    // f4->SetParameter(2,1);
    // f4->SetParameter(3,1);

    f5 = new TF1("f5","[0]*ROOT::Math::erfc(([1]-x)*[2])+[3]*ROOT::Math::erfc((x-[4])*[5])+[6]",210,300);
    // TF1 *f2 = new TF1("f2","[0]+[1]*x",0,10);
    // f4->SetParameter(3,d);
    // f5->SetParameter(0,1);
    // f5->SetParameter(1,1);
    // f5->SetParameter(2,10);
    // f5->SetParameter(3,1);
    // f5->SetParameter(4,1);
    // f5->SetParameter(5,10);
    
    f5->SetParameters(40.0/2,270,1,40.0/2,220,1,5);
    // f5->SetParameters(0,0,0,1.0,150.0,1,5);
    // f5->SetParameter(2,1);
    // f5->SetParameter(3,1);

    // f6 = new TF1("f6","[0]*(exp(4*x*[1]/[2])-1)/(exp((2*x*[1]+[1]*[1])/[2])+exp(4*x*[1]/[2])+1)",10,290);
    f6 = new TF1("f6","[0]*(TMath::Gaus(x,[1]+[4],[2])-TMath::Gaus(x,-[1]+[4],[2]))/(TMath::Gaus(x,[4],[2])+TMath::Gaus(x,[1]+[4],[2])+TMath::Gaus(x,-[1]+[4],[2]))+[3]",150,250);
    f6->SetParameters(0.5,75.0,10.0,-0.1,100.0);

    TF1 *f[20];
    TGraph *gr_fit[20];
    TGraphErrors *gr[20];
    TString legend[20];

    // func->SetParameter(1,c);
    // func->SetParameter(2,d);

    while((entry = (char *)gSystem->GetDirEntry(dataDir)))
    {
        txtFile = entry;
        if(txtFile.Contains("g_"))   continue;
        
        if(txtFile.EndsWith(".txt") || txtFile.EndsWith(".csv") )
        // if(txtFile.EndsWith(".txt") )
        {
            cout<<txtFile<<endl;
            vector<double>  x_list,x_err_list,y_list,y_err_list;
            
            TString txtFileFullPath = dataPath + txtFile;

            TTree *tree = new TTree("graph","graph");
            // tree->ReadFile(txtFileFullPath, "x/D:y/D", ',');
            // tree->ReadFile(txtFileFullPath, "x/D:y/D:x_err/D",' ');
            tree->ReadFile(txtFileFullPath, "x/D:y/D:x_err/D:y_err/D",',');
            tree->SetBranchAddress("x",&x);
            // tree->SetBranchAddress("x",&y_err);
            tree->SetBranchAddress("y",&y);
            tree->SetBranchAddress("x_err",&x_err);
            tree->SetBranchAddress("y_err",&y_err);

            tree->GetEntry(0);

            x0 = x;
            // y0 = 10000000;
            y0 = 0;

            // for(n = (tree->GetEntries())*0.1; n < tree->GetEntries(); n++)
            for(n = 0; n < tree->GetEntries(); n++)
            // for(n = 0; n < (tree->GetEntries())*0.2; n++)
            {
                tree->GetEntry(n);
                // if(y_err > 189) break;
                // if(n*2.5 > 250 || n*2.5 < 100) continue;
                // if(x_err == 1010 )
                // {
                // cout<<x<<","<<y<<","<<x_err<<","<<y_err<<endl;
                    // x_list.push_back(y_err);
                    // y_list.push_back((y0-y)*0.5);
                // }
                    // if(x != x0)
                    // if(TString(to_string(x)).Contains("140720"))  
                    // {
                    //     x0 = 0;
                    //     continue;
                    // }
                    x_list.push_back(x);
                // if(txtFile.Contains("rise"))
                // {
                //     y_list.push_back(y*4);
                // }
                // else   y_list.push_back(y/x_err);
                // {
                    // cout<<x_err<<","<<-x0<<endl;
                    // x_list.push_back(y_err);
                    // cout<<(x-x0)*2.5<<","<<y0/10.0<<endl;
                    // x_list.push_back((x-760)*2.5);
                    y_list.push_back(y);
                    // if(txtFile.Contains("Jitter"))
                    // y_list.push_back(x_err);
                    
                    // x0 = 0;
                    // else
                    // if(txtFile.Contains("USTC"))    y_list.push_back(y*1.1);
                    // else    y_list.push_back(y);
                    // x_err_list.push_back(0);
                    // y_err_list.push_back(0);
                    // x0 = x;
                    // y0 = y;
                // }
                // y0 = y0 - y;
                //     y0 = y0 - y*350;
                // x0 = x;
                // }
                // if(i == 0)  y_list.push_back(y*500/1850);
                // else    y_list.push_back(y*500/2741);
                x_err_list.push_back(x_err);
                y_err_list.push_back(y_err);
            }
            
            tree->Delete("");
            
            // sort(x_list.begin(),x_list.end());

            gr[i] = new TGraphErrors(x_list.size(),&(x_list[0]),&(y_list[0]),&(x_err_list[0]),&(y_err_list[0]));

            //max = gr[i]->GetHistogram()->GetMaximum();
            if(x_max < TMath::MaxElement(gr[i]->GetN(),gr[i]->GetX()))  x_max = TMath::MaxElement(gr[i]->GetN(),gr[i]->GetX());
            if(y_max < TMath::MaxElement(gr[i]->GetN(),gr[i]->GetY()))  y_max = TMath::MaxElement(gr[i]->GetN(),gr[i]->GetY());
            if(x_min > TMath::MinElement(gr[i]->GetN(),gr[i]->GetX()))  x_min = TMath::MinElement(gr[i]->GetN(),gr[i]->GetX());
            if(y_min > TMath::MinElement(gr[i]->GetN(),gr[i]->GetY()))  y_min = TMath::MinElement(gr[i]->GetN(),gr[i]->GetY());
            //x_max = TMath::MaxElement(gr[i]->GetN(),gr[i]->GetX());
            //y_max = TMath::MaxElement(gr[i]->GetN(),gr[i]->GetY());
            //cout<<x_max<<","<<y_max<<endl;

            if(dataPath.Contains("Charge_Bias"))
            {
                cout<<"Charge:"<<Q<<" fC"<<endl;
                gr[i]->Fit(f1,"QCN");
                // cout<<a<<endl;
                f[i] = new TF1("f","expo",0,1000);
                f[i]->SetParameter(0,f1->GetParameter(0));
                f[i]->SetParameter(1,f1->GetParameter(1));
                // f[i]->SetParameter(2,f1->GetParameter(2));
                cout<<txtFile<<",interpolation bias:"<<interpolation(y_list,x_list,Q)<<" V,"<<"exp fit bias:"<<f[i]->GetX(Q)<<" V"<<endl;
                // cout<<txtFile<<",bias:"<< ROOT::Math::Interpolator.Eval(y_list,x_list,Q)<<" V"<<endl;
                // for(n = 0; n < 20; n++)    
                // {
                //     cout<<n<<","<<interpolation(y_list,x_list,Q-2.+0.1*n)<<","<<Q-2.+0.1*n<<endl;
                //     gr_fit[i] = new TGraph();
                //     gr_fit[i]->SetPoint(n,interpolation(y_list,x_list,Q-2.+0.1*n),Q-2.+0.1*n);
                // }
                // gr[i]->SetPoint(n,interpolation(y_list,x_list,Q),Q);
            }

            if(dataPath.Contains("Flux"))
            {
                gr[i]->Fit(f2,"N");
                f[i] = new TF1("f","[2]-expo",0,10);
                f[i]->SetParameter(0,f2->GetParameter(0));
                f[i]->SetParameter(1,f2->GetParameter(1));
                f[i]->SetParameter(2,f2->GetParameter(2));
            }

            // if(dataPath.Contains("Vgl"))
            if(dataPath.Contains("Fitting"))
            // if(txtFile.Contains("fra"))
            {
                gr[i]->Fit(f3,"CN");
                f[i] = new TF1("f","[0]+[1]*x",0,400);
                // f[i] = new TF1("f","[0]*x",-50,50);
                f[i]->SetParameter(0,f3->GetParameter(0));
                f[i]->SetParameter(1,f3->GetParameter(1));
                // cout<<f3->GetParameter(0)<<endl;
            }

            // if(dataPath.Contains("Focus"))
            if(txtFile.Contains("fra"))
            {
                gr[i]->Fit(f4,"CN");
                f[i] = new TF1("f","[0]*ROOT::Math::erfc((x-[1])*[2])+[3]",0,400);
                // f[i] = new TF1("f","[3]*exp(-[2]/(1+(x-[0])*(x-[0])/[1]/[1]))",-10000,-7500);
                // gr[i]->Fit(f6,"CN");
                // f[i] = new TF1("f","[0]*(TMath::Gaus(x,[1]+[4],[2])-TMath::Gaus(x,-[1]+[4],[2]))/(TMath::Gaus(x,[4],[2])+TMath::Gaus(x,[1]+[4],[2])+TMath::Gaus(x,-[1]+[4],[2]))+[3]",0,300);
                f[i]->SetParameter(0,f4->GetParameter(0));
                f[i]->SetParameter(1,f4->GetParameter(1));
                f[i]->SetParameter(2,f4->GetParameter(2));
                f[i]->SetParameter(3,f4->GetParameter(3));
                // f[i]->SetParameter(4,f6->GetParameter(4));
            }

            // if(dataPath.Contains("Profile"))
            if(txtFile.Contains("Profile"))
            {
                gr[i]->Fit(f5,"CNR");
                // cout<<a<<endl;
                f[i] = new TF1("f","[0]*ROOT::Math::erfc(([1]-x)*[2])+[3]*ROOT::Math::erfc((x-[4])*[5])+[6]",200,400);
                f[i]->SetParameter(0,f5->GetParameter(0));
                f[i]->SetParameter(1,f5->GetParameter(1));
                f[i]->SetParameter(2,f5->GetParameter(2));
                f[i]->SetParameter(3,f5->GetParameter(3));
                f[i]->SetParameter(4,f5->GetParameter(4));
                f[i]->SetParameter(5,f5->GetParameter(5));
                f[i]->SetParameter(6,f5->GetParameter(6));
                cout<<f5->GetParameter(1)-f5->GetParameter(4)<<endl;
            }

            // f[i] = new TF1("f","[0]*exp(-[1]*x)",-5,10);
            // f[i]->SetParameter(0,k1[i]);
            // f[i]->SetParameter(1,k2[i]);

            txtFile = txtFile.ReplaceAll(".txt","");
            txtFile = txtFile.ReplaceAll(".csv","");
            legend[i] = txtFile;

            i++;
            // if(i > 4)   break;
        }
    }

    n = i;

    TLegend *leg[row];

    // leg[0] = new TLegend(0.2,0.9,0.4,0.66);
    // leg[1] = new TLegend(0.4,0.9,0.6,0.66);
    // leg[2] = new TLegend(0.6,0.9,0.8,0.66);

    for(i = 0; i < row; i++)    leg[i] = new TLegend(0.2+i*0.2,0.9,0.4+i*0.2,0.9-column*0.03);
    
    // for(i = 0; i < row; i++)    leg[i] = new TLegend(0.45+i*0.2,0.81,0.65+i*0.2,0.81-column*0.063);

    if(dataPath.Contains("Charge_Bias"))    
    {
        plotname = "Charge_Bias/" + plotname;
        gr[0]->SetTitle(";Bias [V];Charge [fC]");
    }
    else if(dataPath.Contains("TRes_Bias"))   
    {
        plotname = "TRes_Bias/" + plotname;
        gr[0]->SetTitle(";Bias [V];Time Resolution [ps]");
    }

    else if(dataPath.Contains("TRes_Charge"))
    {
        plotname = "TRes_Charge/" + plotname;
        gr[0]->SetTitle(";Charge [fC];Time Resolution [ps]");
    }
    
    else
    {
        plotname = "Others/" + plotname;
        // leg = new TLegend(0.2,0.9,0.4,0.9-n*0.03);
        // leg = new TLegend(0.64,0.9,0.84,0.9-n*0.03);
        // gr[0]->SetTitle(";Injected charge [fC];Pulse integral [mV*ns]");
        // gr[0]->SetTitle(";Fluence [#times10^{15} cm^{-2}]; k");
        // gr[0]->SetTitle(";Charge [fC];FWHM [fC]");
        // gr[0]->SetTitle(";Y [um]; Ampl [a.u.]");
        gr[0]->SetTitle(";Y [um]; Fraction");
        // gr[0]->SetTitle(";Nominal IP [um];Effective IP [um]");
        // gr[0]->SetTitle(";Bias [V];Leakage [A]");
        // gr[0]->SetTitle(";Bias [V]; RMS [mV]");
        // gr[0]->SetTitle(";False positive rate;True positive rate");
        // gr[0]->SetTitle(";Time [ns];Voltage [a.u.]");
        // gr[0]->SetTitle(";Z' Mass [GeV];Expected limit on xsec [fb]");
        // gr[0]->SetTitle(";#sigma_{DUT} [ps];#sqrt{V}(#sigma_{DUT})(Ref) [ps]");
    }

    gr[0]->GetYaxis()->SetRangeUser(y_min-(y_max-y_min)/10, y_max+(y_max-y_min)/2);
    // gr[0]->GetYaxis()->SetRangeUser(y_min-(y_max-y_min)/10, 50);
    // gr[0]->GetYaxis()->SetRangeUser(20, y_max+(y_max-y_min)/10);
    // gr[0]->GetYaxis()->SetRangeUser(1e-3, 5);
    // gr[0]->GetYaxis()->SetRangeUser(1e-1, 100);
    gr[0]->GetXaxis()->SetLimits(x_min-(x_max-x_min)/10, x_max+(x_max-x_min)/10);
    // gr[0]->GetXaxis()->SetLimits(1e-1, x_max+(x_max-x_min)/10);
    // gr[0]->GetXaxis()->SetLimits(x_min-(x_max-x_min)/10, 10);
    // gr[0]->GetXaxis()->SetLimits(0, 1000);
    gr[0]->GetXaxis()->SetNdivisions(505, kTRUE);
    gStyle->SetEndErrorSize(5);

    // c1->SetLogy();
    // c1->SetLogx();
    int color = 0;
    for(i = 0; i < n; i++)
    {   
        if(legend[i].Contains("W12")) color = 3;
        else if(legend[i].Contains("W19")) color = 0;
        else if(legend[i].Contains("W21")) color = 2;
        else if(legend[i].Contains("W17")) color = 1;
        else if(legend[i].Contains("W9")) color = 4;
        else if(legend[i].Contains("W10")) color = 2;
        else if(legend[i].Contains("HPK")) color = 6;
        else if(legend[i].Contains("W7")) color = 3;
        else if(legend[i].Contains("W8")) color = 1;
        else if(legend[i].Contains("W11")) color = 0;
        // else if(legend[i].Contains("-10")) color = 0;
        // else if(legend[i].Contains("-20")) color = 1;
        // else if(legend[i].Contains("30")) color = 2;
        // else if(legend[i].Contains(",10")) color = 3;
        // else if(legend[i].Contains(",20")) color = 4;
        else    color = 0;
        
        if(legend[i].Contains("JSI")) color = 6;
        if(legend[i].Contains("USTCB")) color = 0;
        if(legend[i].Contains("UCSCB")) color = 1;

        // gr[i]->SetMarkerColor(kBlue);
        // gr[i]->SetLineColor(kBlue);
        // gr[i]->SetMarkerColor(colorList[color]);
        // gr[i]->SetLineColor(colorList[color]);
        gr[i]->SetMarkerColor(colorList[i]);
        gr[i]->SetLineColor(colorList[i]);
        // gr[i]->SetLineStyle(10 - i);
        gr[i]->SetMarkerStyle(20 + i);
        // gr[i]->SetMarkerStyle(20 + color);
        // gr[i]->SetMarkerSize(1 - 0.2*i);
        gr[i]->SetMarkerSize(2);
        gr[i]->SetLineWidth(2);
        if(legend[i].Contains(",M"))    gr[i]->SetLineStyle(10);
        if(i == 0)  gr[i]->Draw("AP");
        // else if(i == 4 || i ==5)   gr[i]->Draw("LP");
        else    gr[i]->Draw("P");

        // if(dataPath.Contains("Charge_Bias")||dataPath.Contains("Flux")||dataPath.Contains("Vgl"))

        // string A0 = "c(W" + to_string(i+7) + ")=" + to_string(k2[i]*10).substr(0,4) + "#times10^{-16}";
        // cout<<"here"<<endl;
        type = legend[i];
        // cout<<type<<end;
        // cout<<TString("Charge="+to_string(f3->GetParameter(0))+"+"+to_string(f3->GetParameter(1))+"*Ampl")<<endl;
        // leg->SetTextColor(colorList[i]);
        leg[i/column]->SetTextSize(0.03);
        leg[i/column]->SetFillStyle(0);
        leg[i/column]->SetTextFont(42);
        leg[i/column]->AddEntry(gr[i],type,"P");
        // leg[i/column]->AddEntry(gr[i],TString(" "),"P");
        // leg[i/column]->AddEntry(gr[i],TString("Charge="+to_string(f3->GetParameter(0))+"+"+to_string(f3->GetParameter(1))+"*Ampl"),"P");
        leg[i/column]->Draw("same");
        
        pt=new TPaveText(0.5,0.95-i*0.15,0.7,0.65,"NDC");
        // pt=new TPaveText(0.6,0.9,0.7,0.85,"NDC");
        pt->SetBorderSize(0);
        // pt->SetTextColor(kBlack);
        pt->SetTextColor(colorList[i]);
        pt->SetFillStyle(0);
        pt->SetTextAlign(13);
        pt->SetTextSize(0.05);
        
        // string A0 =  "k"+to_string(i)+"="+ to_string(f[i]->GetParameter(0)).substr(0,4);
        // string A0 =  "y = "+ to_string(f3->GetParameter(0)) + "+" + to_string(f3->GetParameter(1)) + "*x";
        // string A1 =  "B = "+ to_string(f4->GetParameter(1)) + "#pm" + to_string(f4->GetParError(1));
        // string A2 =  "C = "+ to_string(f4->GetParameter(2)) + "#pm" + to_string(f4->GetParError(2));
        // string A3 =  "D = "+ to_string(f4->GetParameter(3)) + "#pm" + to_string(f4->GetParError(3));

        // text=pt->AddText("#it{HPK Type 3.1}");
        // text=pt->AddText("#it{USTC W5}");
        // text=pt->AddText("#font[12]{Ampl} = D*exp(-#frac{C}{1+(#frac{#font[12]{Z}-A}{B})^{2}})");
        // text=pt->AddText("");
        // text=pt->AddText(TString(A1));
        // text=pt->AddText(TString(A2));
        // text=pt->AddText(TString(A3));
        // pt->Draw("same");
        if(fit)
        {
            // string A0 = "c(W" + to_string(i+7) + ")=" + to_string(k2[i]*10).substr(0,4) + "#times10^{-16}";
            // text=pt->AddText(TString(A0));
            // f[i]->SetLineColor(colorList[color]);
            // f[i]->SetLineColor(colorList[i]);
            f[i]->SetLineColor(kBlack);
            f[i]->SetLineWidth(2);
            f[i]->Draw("same");
            // TGraph *g = (TGraph*)f[i]->DrawDerivative("same");
            
            // vector<double> z_list;
            // for(j = 0; j < g->GetMaxSize(); j++)    z_list.push_back(-1.0/(g->GetY()[j]));
            // // cout<<z_list[50]<<endl;
            // g = new TGraph(g->GetMaxSize(),g->GetX(),&(z_list[0]));
            // g->SetLineColor(colorList[i]);
            // g->Draw("same");
            
            // cout<<g->GetMaxSize()<<endl;
            // f[i]->DrawDerivative("same");
            // cout<<"fitting"<<endl;
        }
        // pt->Draw("same");
    }

    c1->SaveAs(plotStorePath + plotname + ".png");
    c1->SaveAs(plotStorePath + plotname + ".svg");
    // c1->SaveAs(plotStorePath + plotname + ".eps");
    c1->SaveAs(plotStorePath + plotname + ".pdf");
 
    c1->Clear();
}

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