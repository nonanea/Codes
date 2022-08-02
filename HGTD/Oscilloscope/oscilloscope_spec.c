// Script for Peak/Charge/Tot.. Spectrum visualiztion and analysis
// Made by USTC HGTD Group, September 17,2019

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TPad.h>
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

//const TString dataPath = "Data/";
const TString plotStorePath = "Plots/";
const TString ext = ".txt";
TLegend *legend;
int draw_legend(TString text);
void PedestalAndNoise(vector<double> ylist,double *ped,double *noise);
void Peak(vector<double> xlist, vector<double> Alist, double *Tmax, double *Amax, int *Max, double threshold);
double Charge(vector<double> xlist, vector<double> Alist, int Max, double threshold);
double TOT(vector<double> xlist, vector<double> Alist, double threshold);
double Time, Ampl;
bool multi = true;
vector<int> colorList{kRed, kOrange-3,kGreen +2, kBlue, kViolet,kTeal-1, kMagenta, kViolet, kAzure, kCyan, 
                    kBlue + 1, kRed + 1, kGreen + 1, kTeal + 1, kPink + 1, kMagenta + 1, kViolet + 1, kAzure + 1, kCyan + 1, kYellow + 1};


void oscilloscope_spec(TString dataPath, bool multi = true)
{
    
    gErrorIgnoreLevel = kError;
        
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(.03);
    TText *t = new TText();

    TString plotType = dataPath;
    TString plotname = dataPath;
    plotname.ReplaceAll("/",".");
    dataPath = "Data/" + dataPath;
    void *plotDir = gSystem->OpenDirectory(dataPath);
    
    // Check Plot directory exist or not, if not, create it
    if (!gSystem->OpenDirectory(plotStorePath + plotType))
        gSystem->MakeDirectory(plotStorePath + plotType);

    //cout<<dataPath<<endl;
    TString txtFile,subDir;
    const char *entry, *subEntry;
    //TString plotInfo (dataPath);
    TString plotTarget ="Amp";// Charge / TOT
    
    int i = 0;
    int k = 0;
    int count = 20;
    double threshold = 0;
    //vector<double> Tlist,Alist;
    int n = 0, a=0, nDir=-1, b=0, Max = 0;
    double ped = 0, noise = 0, Tmax = 0, Amax = 0;
    
    vector<TString> dirList;
    while ((entry = (char *)gSystem->GetDirEntry(plotDir)))
    {
        if(TString(entry).Contains("."))    continue;
        //if(!TString(entry).Contains("HV120"))continue;
        dirList.push_back(entry);
    } 
    sort(dirList.begin(), dirList.end());
    //Sort accroding to Directory Name
    
    vector<TH1D*> hPeakList,hChargeList,hChargeCalList,hTOTList;
    
    ofstream myfile;
    myfile.open(plotStorePath + plotType + "Results.csv");
    
    for(int it=0;it<dirList.size();it++)
    {
        entry = dirList[it];     
        void *subDir = gSystem->OpenDirectory(dataPath+entry);    
        TString subDir_str = entry;
        nDir++;
        TH1D *hpeak = new TH1D(Form("peak of " + plotname+ "%d" ,a++), "peak of " + plotname, 100, 0, 100);
        TH1D *hCharge = new TH1D(Form("Charge of " + plotname+ "%d" ,a++), "charge of " + plotname, 100, 0, 600);
        TH1D *hChargeCal = new TH1D(Form("Charge of " + plotname+ "%d" ,a++), "charge of " + plotname, 100, 0, 60);
        TH1D *hTOT = new TH1D(Form("TOT of " + plotname+ "%d" ,a++), "TOT of " + plotname, 50, 0, 10e-9);
        TGraph *gr;

        while ((subEntry = (char *)gSystem->GetDirEntry(subDir)))
        {        
            //if(i++>10)break;
            //cout<<i<<endl;
            txtFile = subEntry;

            vector<double> xlist,ylist,Tlist,Alist;
            if (!txtFile.EndsWith("txt")) continue;
            if (txtFile.Contains("C2"))
            {
                // cout<<txtFile<<endl;
                TString txtFileFullPath = dataPath + entry + "/" +txtFile;
                TString plotInfo = txtFile;
                plotInfo = plotInfo.ReplaceAll(".txt", "");
                TTree *tree = new TTree("signal", "signal of oscillscope");
                tree->ReadFile(txtFileFullPath, "Time/D:Ampl/D", ',');
                //cout<<"txtFile:"<<txtFileFullPath<<endl;
                // Start To read file into TTree
                tree->SetBranchAddress("Time", &Time);
                tree->SetBranchAddress("Ampl", &Ampl);
                gr = new TGraph();
                for (n = 0; n < tree->GetEntries(); n++)
                {
                    tree->GetEntry(n);
                    // cout << Time << "  " << Ampl << " " << endl;
                    xlist.push_back(Time);
                    ylist.push_back(Ampl);
                    // Alist.push_back(Ampl);
                }

                 PedestalAndNoise(ylist,&ped,&noise);
                //  cout<<"ped="<<ped<<endl;
                // cout<<"noise="<<noise<<endl;

                //deduct the pedestal
                for( n = 0; n < ylist.size(); n++){
                    Alist.push_back(ylist[n]-ped);
                    gr->SetPoint(n, xlist[n], Alist[n]);
                }
                
                //get the peak information
                Peak(xlist,Alist,&Tmax,&Amax,&Max,threshold);
                // cout<<"Max="<<Max<<endl;
                //cout<<"Tmax="<<Tmax<<endl;
                double charge =  Charge(xlist,Alist,Max,threshold);
                //cout<<"charge:"<<charge<<endl;
                
                //Cut Here
                //if(Amax*1000 > 500) continue;
                
                hpeak->Fill(Amax*1000);
                //hpeak->Fill(charge*300);
                //Unit[V*ns / 1000 = V*ps]
                hCharge->Fill(charge*1000);
                hChargeCal->Fill(charge*1000.0/7.3733);
                hTOT->Fill(TOT(xlist,Alist,threshold));
                myfile<<TString(entry)<<","<<Amax*1000<<","<<charge*1000<<endl;
            }
        }
     
        hPeakList.push_back(hpeak);
        hChargeList.push_back(hCharge);
        hChargeCalList.push_back(hChargeCal);
        hTOTList.push_back(hTOT);
    }
    
    TString title = plotname;
    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 1000, 800);
//     c1->SetLeftMargin(0.25);
//     c1->SetRightMargin(0.2);
//     c1->SetTopMargin(0.2);
    for(int ii=0;ii<hPeakList.size();ii++)
    { 
        TH1D *hpeak = hPeakList[ii];
        bool doFit = false;
     
        gStyle->SetOptStat("n");

        /*
            DRAW AMP
        */
     
        hpeak->GetXaxis()->SetTitle("Amplitude [mV]");
        hpeak->GetYaxis()->SetTitle("Entries");
        int nentries = hpeak->GetEntries();
        hpeak->Scale(1.0/hpeak->Integral());

        hpeak->SetLineColor(colorList[ii]);
        hpeak->Draw("same hist");
        hpeak->SetAxisRange(0., 0.15,"Y");
        if(dirList[ii].Contains("fC")){
            hpeak->SetAxisRange(0., 1.5,"Y");
            doFit=false;
        }
        double Constant=0,MPV=0,Sigma=0,MPVErr=0; 
        if(doFit){
            TFitResultPtr r = hpeak->Fit("landau","S","C",0,400);
            Constant = r->Value(0);
            MPV = r->Value(1);
            Sigma = r->Value(2);
            MPVErr = r->Error(1);
            cout<<"Constant:"<<Constant<<"  MPV:"<<MPV<<"    Sigma:"<<Sigma<<endl;
            TF1* landau = new TF1("f1","landau",0,700);
            landau->SetParameter(0,Constant);
            landau->SetParameter(1,MPV);
            landau->SetParameter(2,Sigma);
            landau->SetLineColor(colorList[ii]);
            landau->Draw("same");
        }
        legend = new TLegend(0.4, 0.95-ii*0.03 , 0.7, 0.7-ii*0.03 );
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetTextSize(0.02);
        legend->SetTextColor(colorList[ii]);
        legend->SetFillStyle(0);
        legend->SetHeader(Form(dirList[ii]+", MPV:%3.2f,MPV_Err:%2.2f%%,Mean:%.2f, Entries:%d",MPV,100*MPVErr/MPV,hpeak->GetMean(),nentries), "C");
        legend->Draw();
        
        cout<<"draw Peak for" << dirList[ii]  << " ii:"<<ii<<endl;       

    }
    //draw_legend(title);
    gErrorIgnoreLevel = kPrint;
    c1->SaveAs(plotStorePath + plotType + "peak_Dist.png");
    c1->SaveAs(plotStorePath + plotType + "peak_Dist.pdf");
    delete c1;
    
    c1 = new TCanvas("c1", "Graph Draw Options", 1000, 800);
    for(int ii=0;ii<hChargeList.size();ii++)
    { 
        TH1D *hcharge = hChargeList[ii];     
        /*
            DRAW CHARGE
        */
        double Constant=0,MPV=0,Sigma=0,MPVErr=0;     
        hcharge->GetXaxis()->SetTitle("Voltage.Integral [mV*ns]"); //V*ps == mV*ns
        hcharge->GetYaxis()->SetTitle("Entries");
        hcharge->Scale(1.0/hcharge->Integral());
        hcharge->SetLineColor(colorList[ii]);
        hcharge->Draw("same hist");
        if(int(hcharge->GetEntries())==1)
            hcharge->SetAxisRange(0., 1.5,"Y");
        else
            hcharge->SetAxisRange(0., 0.15,"Y");
        
        
        legend = new TLegend(0.4, 0.95-ii*0.03 , 0.7, 0.7-ii*0.03 );
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetTextSize(0.02);
        legend->SetTextColor(colorList[ii]);
        legend->SetFillStyle(0);
        legend->SetHeader(Form(dirList[ii]+
           ", MPV:%3.2f,    MPV_Err:%2.2f%%,    Mean:%.2f,          Entries:%d",
           MPV,             100*MPVErr/MPV,     hcharge->GetMean(), int(hcharge->GetEntries())), "C");     
        legend->Draw();
          
     //double center = hpeak->GetBinCenter(hpeak->GetMaximumBin());
     //hpeak->GetXaxis()->SetTitle("Prop.to Voltage.Integral [V*s]");
     
        // myfile<<hcharge->GetMean()<<endl;

        cout<<"draw Charge for" << dirList[ii]  << " ii:"<<ii<<endl;
    }
    //draw_legend(title);
    gErrorIgnoreLevel = kPrint;
    c1->SaveAs(plotStorePath + plotType + "charge_Dist.png");
    c1->SaveAs(plotStorePath + plotType + "charge_Dist.pdf");
    delete c1;
    
    c1 = new TCanvas("c1", "Graph Draw Options", 1000, 800);
    for(int ii=0;ii<hChargeCalList.size();ii++)
    { 
        TH1D *hcharge = hChargeCalList[ii];     
        /*
            DRAW CHARGE
        */
        double Constant=0,MPV=0,Sigma=0,MPVErr=0;     
        hcharge->GetXaxis()->SetTitle("Collected Charge [fC]"); //V*ps == mV*ns
        hcharge->GetYaxis()->SetTitle("Entries");
        hcharge->Scale(1.0/hcharge->Integral());
        hcharge->SetLineColor(colorList[ii]);
        hcharge->Draw("same hist");
        if(int(hcharge->GetEntries())==1)
            hcharge->SetAxisRange(0., 1.5,"Y");
        else    
            hcharge->SetAxisRange(0., 0.15,"Y");
        
        bool doFit = false;
        if(dirList[ii].Contains("fC")){
            hcharge->SetAxisRange(0., 1.5,"Y");
            doFit = false;
        }
        if(doFit){
            TFitResultPtr r = hcharge->Fit("landau","S","C",0,400);
            Constant = r->Value(0);
            MPV = r->Value(1);
            Sigma = r->Value(2);
            MPVErr = r->Error(1);
            cout<<"Constant:"<<Constant<<"  MPV:"<<MPV<<"    Sigma:"<<Sigma<<endl;
            TF1* landau = new TF1("f1","landau",0,700);
            landau->SetParameter(0,Constant);
            landau->SetParameter(1,MPV);
            landau->SetParameter(2,Sigma);
            landau->SetLineColor(colorList[ii]);
            landau->SetLineWidth(1);
            landau->Draw("same");
        }
                
        
        legend = new TLegend(0.4, 0.95-ii*0.03 , 0.7, 0.7-ii*0.03 );
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetTextSize(0.02);
        legend->SetTextColor(colorList[ii]);
        legend->SetFillStyle(0);
        legend->SetHeader(Form(dirList[ii]+
           ", MPV:%3.2f,    MPV_Err:%2.2f%%,    Mean:%.2f,          Entries:%d",
           MPV,             100*MPVErr/MPV,     hcharge->GetMean(), int(hcharge->GetEntries())), "C");     
        legend->Draw();
          
     //double center = hpeak->GetBinCenter(hpeak->GetMaximumBin());
     //hpeak->GetXaxis()->SetTitle("Prop.to Voltage.Integral [V*s]");
     

        cout<<"draw Charge(Cal) for" << dirList[ii] << " ii:"<<ii<<endl;
    }
    //draw_legend(title);
    gErrorIgnoreLevel = kPrint;
    c1->SaveAs(plotStorePath + plotType + "charge_Dist_Cal.png");
    c1->SaveAs(plotStorePath + plotType + "charge_Dist_Cal.pdf");
    delete c1;  
    
    c1 = new TCanvas("c1", "Graph Draw Options", 1000, 800);
    for(int ii=0;ii<hTOTList.size();ii++)
    { 
        TH1D *hcharge = hTOTList[ii];     
        /*
            DRAW CHARGE
        */
        double Constant=0,MPV=0,Sigma=0,MPVErr=0;     
        hcharge->GetXaxis()->SetTitle("TOT [ns]"); //V*ps == mV*ns
        hcharge->GetYaxis()->SetTitle("Entries");
        hcharge->Scale(1.0/hcharge->Integral());
        hcharge->SetLineColor(colorList[ii]);
        hcharge->Draw("same hist");
        if(int(hcharge->GetEntries())==1)
            hcharge->SetAxisRange(0., 1.5,"Y");
        else    
            hcharge->SetAxisRange(0., 0.15,"Y");
        
        legend = new TLegend(0.4, 0.95-ii*0.03 , 0.7, 0.7-ii*0.03 );
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetTextSize(0.02);
        legend->SetTextColor(colorList[ii]);
        legend->SetFillStyle(0);
        legend->SetHeader(Form(dirList[ii]+
           ", MPV:%3.2f,    MPV_Err:%2.2f%%,    Mean:%.2f,          Entries:%d",
           MPV,             100*MPVErr/MPV,     hcharge->GetMean(), int(hcharge->GetEntries())), "C");     
        legend->Draw();
          
     //double center = hpeak->GetBinCenter(hpeak->GetMaximumBin());
     //hpeak->GetXaxis()->SetTitle("Prop.to Voltage.Integral [V*s]");
     

        cout<<"draw TOT for" << dirList[ii] << " ii:"<<ii<<endl;
    }
    //draw_legend(title);
    gErrorIgnoreLevel = kPrint;
    c1->SaveAs(plotStorePath + plotType + "TOT_Dist_Cal.png");
    c1->SaveAs(plotStorePath + plotType + "TOT_Dist_Cal.pdf");
    delete c1;  

}

int draw_legend(TString text)
{
    //cout << "drawing title" << endl;
    legend = new TLegend(0.3, 0.9, 0.75, 0.8);
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);
    legend->SetHeader(text, "C");
    legend->Draw();
    return 0;
}

void PedestalAndNoise(vector<double> ylist,double *ped,double *noise)
{
    int i;
    double mean = 0,stdev = 0;
    // cout<<ylist.size()<<endl;
    for(i = 0; i < ylist.size()/4; i++){
        mean+=ylist[i];
    }
    mean=mean/(i+1);
    // cout<<mean<<endl;

    for(i = 0; i < ylist.size()/4; i++){
        stdev+=(ylist[i]-mean)*(ylist[i]-mean);
    }
    stdev=sqrt(stdev/i);
    *ped=mean;
    *noise=stdev;
}


void Peak(vector<double> xlist, vector<double> Alist, double *Tmax, double *Amax, int *Max, double threshold)
{
    int k = 0;
    double amax = threshold;

    // m = xlist.size()/2;
    //cout<<"m="<<m<<endl;
    if (xlist.size() != 0)
    {
        for(k = 0; k < xlist.size(); k++)    
        {
            if(Alist[k] > amax)
            {
                amax = Alist[k];
                *Tmax = xlist[k];
                *Amax = Alist[k];
                *Max = k;
                //break;
            }
        } 
    }
}

double Charge(vector<double> xlist, vector<double> Alist, int Max, double threshold)
{
   
    int k = 0, length = 0;
    double tot = 0, integral_voltage = 0, delta_T = 0;
    double window = 0;
    length = xlist.size();
    window = xlist[length-1] - xlist[0];
    delta_T = window*1e9/(length-1);//ns

    cout<<"delta_T="<<delta_T<<endl;
    if (length != 0)
    {
        for(k = Max - 5/delta_T; k < Max + 5/delta_T + 1; k++)
        {
            // if(Alist[k] > threshold)
            // {
                tot +=2e-10;
                integral_voltage+=Alist[k]*delta_T; // 0.2ns per sample(5GHz)
                //cout<<"xlist[k]"<<xlist[k]<<"   dt:"<<xlist[k+1]-xlist[k]<<endl;
            // }
        } 
    }
    //Unit: V*ns
    return integral_voltage ;   
}

double TOT(vector<double> xlist, vector<double> Alist, double threshold){

    
    int k = 0, m = 0;
    double tot = 0, integral_voltage = 0 ;
    m = xlist.size()/2;
    //cout<<"m="<<m<<endl;
    if (m != 0)
    {
        for(k = 0; k < xlist.size(); k++)
        {
            if(Alist[k] > threshold)
            {
                tot +=2e-10;
                //integral_voltage+=Alist[k]*0.2; // 0.2ns per sample(5GHz)
                //cout<<"xlist[k]"<<xlist[k]<<"   dt:"<<xlist[k+1]-xlist[k]<<endl;
            }
        } 
    }
    //Unit: V*ns
    return tot ;   
}



int main(int argc, char **argv)
{
    if (argc == 1 || argc > 2)
        cout << "Usage: " << argv[0] << " dataPath (Image saved in ./Plots/)";
    else
    {
        TString dataPath(argv[1]);
        if (!dataPath.EndsWith("/"))
            dataPath = dataPath + "/";
        // string dataPath (argv[1]);
        cout << dataPath << endl;
        oscilloscope_spec(dataPath,multi);
    }
    return 0;
}
