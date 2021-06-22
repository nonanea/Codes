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
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"

using namespace std;

const TString plotStorePath = "Plots/";

//bool dual_board = true;
double threshold_1, threshold_2, t_cut, t_delta;
// double transimpedance = 1.0;
// double transimpedance = 4.3226;//UCSCv1.1
double transimpedance = 4.6929;//UCSCv1.4
// double transimpedance = 16.571;//USTCv2
// double transimpedance = 18.698;//USTCv2
// double res_ref = 37.7;//T3.1
// double res_ref_err = 0.83;
double res_ref = 31.;//T1.1
double res_ref_err = 1.0;
const int N = 20;

double langaufun(double *x, double *par)
{
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

    // Numeric constants
    double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    double mpshift  = -0.22278298;       // Landau maximum location

    // Control constants
    double np = 200.0;      // number of convolution steps
    double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

    // Variables
    double xx;
    double mpc;
    double fland;
    double sum = 0.0;
    double xlow,xupp;
    double step;
    double i;

    // MP shift correction
    mpc = par[1] - mpshift * par[0]; 

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];

    step = (xupp-xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);

        xx = xupp - (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }
    return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1D *his, double *fitrange, double *startvalues, double *fitparams, double *fiterrors, double *ChiSqr, int *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

    int i;

    TF1 *ffit = new TF1("langau",langaufun,fitrange[0],fitrange[1],4);
    ffit->SetParameters(startvalues);
    ffit->SetParNames("Width","MP","Area","GSigma");

    his->Fit("langau","QRB0");   // fit within specified range, use ParLimits, do not plot

    ffit->GetParameters(fitparams);    // obtain fit parameters
    for (i=0; i<4; i++) {
        fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
    }
    ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
    NDF[0] = ffit->GetNDF();           // obtain ndf

    return (ffit);              // return fit function
}


int langaupro(double *params, double &maxx, double &FWHM)
{
   // Seaches for the location (x value) at the maximum of the 
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   double p,x,fy,fxr,fxl;
   double step;
   double l,lold;
   int i = 0;
   int MAXCALLS = 10000;

   // Search for maximum
   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = langaufun(&x,params);
 
      if (l < lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;

   // Search for right x location of fy
   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;

   // Search for left x location of fy
   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-3);

   fxl = x;

   FWHM = fxr - fxl;
   return (0);
}

void file_read(TString filePath, double threshold_1 = 0.0, double threshold_2 = 0.0, double t_cut = -2.5, double t_delta = 1.4)
{
    gErrorIgnoreLevel = kError;

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(.03);
    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 1000, 800);
    
    //TText *t = new TText();
    c1->SetLeftMargin(0.25);
    c1->SetRightMargin(0.2);
    c1->SetTopMargin(0.2);

    TString dirname = filePath;
    dirname.ReplaceAll("/","");
    filePath = plotStorePath + filePath;
    
    cout<<dirname<<endl;
    
    // Check Plot directory exist or not, if not, create it
    if (!gSystem->OpenDirectory(filePath + dirname))
        gSystem->MakeDirectory(filePath + dirname);

    cout<<filePath+dirname<<endl;
    
    //string target_dir, target_dir_tail;

    //target_dir = filePath;
    //target_dir_tail = target_dir.substr(6);

    //int position = target_dir_tail.find_first_of("/");
    
    //target_dir = target_dir.substr(0,6 + position);
    //cout<<position<<","<<target_dir_tail<<","<<target_dir<<endl;

    ofstream myfile;
    myfile.open (filePath + dirname + "/results.txt");
    myfile<<filePath + dirname<<endl;

    TFile *f= new TFile(filePath+dirname+".root","read");
    // TTree *t[4];
    // t[0] = (TTree*) f->Get("t1");
    // t[1] = (TTree*) f->Get("t2");
    // t[2] = (TTree*) f->Get("t3");
    // t[3] = (TTree*) f->Get("t4");

    TTree *t;
    t = (TTree*) f->Get("tt");
    
    TFitResultPtr r;
    TPaveText *pt;
    TText *text;
    TGraph* gr;
    gr = new TGraph();

    int dual_board = 1, event = 0;
    int n = 0, i = 0;
    Float_t Amax[4], Tmax[4], noise[4], noise_amp[i], charge[4], charge_hw[4], charge_dw[4], CTD[4], CFD[4], ZCD[4], TOA_F[4], TOT[4], nentries, ped[4], risetime_ave[4], risetime[4], risetime_F[4], slope[4], linearity, delta_amax, jitter[4], ped_1[4], ped_2[4], slope_F[4], intercept_F[4], id[4];
    //double Amax1, Tmax1, charge1, ped1, risetime1, slope1, noise1, CTD1, CFD1, ZCD1, TOF_CTD, TOF_ZCD, TOA_Ref;
    double risetime_ave_1, risetime_ave_2, jitter_ave_1, jitter_ave_2;
    double TOF_CFD_12, TOF_CFD_13, TOF_CFD_23;
    double sigma_12,sigma_13,sigma_23,sigma_1,sigma_2,sigma_3, sigma_23_err;
    double MPV = 0, MPV_err = 0, mean = 0, mean_err = 0, sigma = 0, sigma_err = 0;
    double stat_err = 0, ref_err = 0, cali_err = 0;
    double TOF_corr;
    double simple_risetime[4];
    //double a = 0, b = 0, c = 0;
    //double a = -0.0087, b = 3.2514, c = -186.6;
    double a = -0.0009, b = 1.7022, c = -686.77;
    double EoE = 0, TOF_mean = 3948.56;

    TF1 *fit_langau;
    double chisqr, langau_Peak,langau_FWHM;
    int ndf;

    int x[10000];
    double y[10000];

    // for(i = 0; i < 4; i++)
    // {
    //     t[i]->SetBranchAddress("Amax",&Amax[i]);   
    //     t[i]->SetBranchAddress("Tmax",&Tmax[i]);
    //     t[i]->SetBranchAddress("pedestal",&ped[i]);
    //     // t[i]->SetBranchAddress("ped_1",&ped_1[i]);
    //     // t[i]->SetBranchAddress("ped_2",&ped_2[i]);
    //     t[i]->SetBranchAddress("Noise",&noise[i]);
    //     t[i]->SetBranchAddress("Noise_Amp",&noise_amp[i]);
    //     t[i]->SetBranchAddress("Charge",&charge[i]);
    //     t[i]->SetBranchAddress("Charge_HW",&charge_hw[i]);
    //     t[i]->SetBranchAddress("Charge_DW",&charge_dw[i]);
    //     t[i]->SetBranchAddress("CTD",&CTD[i]);
    //     t[i]->SetBranchAddress("CFD",&CFD[i]);
    //     // t[i]->SetBranchAddress("ZCD",&ZCD[i]);
    //     // t[i]->SetBranchAddress("TOA_F",&TOA_F[i]);
    //     t[i]->SetBranchAddress("risetime_ave",&risetime_ave[i]);
    //     t[i]->SetBranchAddress("risetime",&risetime[i]);
    //     // t[i]->SetBranchAddress("risetime_F",&risetime_F[i]);
    //     // t[i]->SetBranchAddress("slope",&slope[i]);
    //     t[i]->SetBranchAddress("ID",&id[i]);
    //     // t[i]->SetBranchAddress("slope_F",&slope_F[i]);
    //     // t[i]->SetBranchAddress("intercept_F",&intercept_F[i]);
    // }
    
    t->SetBranchAddress("Amax",Amax);   
    t->SetBranchAddress("Tmax",Tmax);
    t->SetBranchAddress("pedestal",ped);
    // t->SetBranchAddress("ped_1",ped_1);
    // t->SetBranchAddress("ped_2",ped_2);
    t->SetBranchAddress("Noise",noise);
    t->SetBranchAddress("Noise_Amp",noise_amp);
    t->SetBranchAddress("Charge",charge);
    t->SetBranchAddress("Charge_HW",charge_hw);
    t->SetBranchAddress("Charge_DW",charge_dw);
    t->SetBranchAddress("CTD",CTD);
    t->SetBranchAddress("CFD",CFD);
    // t->SetBranchAddress("ZCD",ZCD);
    t->SetBranchAddress("TOT",TOT);
    t->SetBranchAddress("risetime_ave",risetime_ave);
    t->SetBranchAddress("risetime",risetime);
    // t->SetBranchAddress("risetime_F",risetime_F);
    // t->SetBranchAddress("slope",slope);
    t->SetBranchAddress("ID",id);
    // t->SetBranchAddress("slope_F",slope_F);
    // t->SetBranchAddress("intercept_F",intercept_F);

    double fitting[N] = {0}, error[N] = {0};
    const char *plotname[N] = {"nentries",//0
                                "amplitude_1",//1   DUT1
                                "pedstal_1",//2
                                "RMS_1",//3
                                "charge_1",//4
                                "risetime_1",//5
                                "slope_1",//6
                                "jitter_1",//7
                                "amplitude_2",//8   DUT2
                                "pedstal_2",//9
                                "RMS_2",//10
                                "charge_2",//11
                                "risetime_2",//12
                                "slope_2",//13
                                "jitter_2",//14
                                "tof_cfd_12",//15
                                "tof_cfd_13",//16
                                "tof_cfd_23",//17
                                //"tof_corr"//18
                                "noise_amplitude_1",//18
                                "noise_amplitude_2"//19
                                };
    const char *xaxis_name[N] = {"nentries",//0
                                "Ampl_1 [mV]",//1   DUT1
                                "Pedestal_1 [mV]",//2
                                "Noise RMS_1 [mV]",//3
                                "Collected charge_1 [fC]",//4
                                "risetime_1 [ps]",//5
                                "slope_1 [mV/ps]",//6
                                "jitter_1 [ps]",//7
                                "Ampl_2 [mV]",//8   DUT2
                                "Pedestal_2 [mV]",//9
                                "Noise RMS_2 [mV]",//10
                                "Collected charge_2 [fC]",//11
                                "risetime_2 [ps]",//12
                                "slope_2 [mV/ps]",//13
                                "jitter_2 [ps]",//14
                                // "TOF_CFD_12 [ps]",//15
                                "Amplitude*TOT [mV*ns]",//15
                                "TOT [ps]",//16
                                "TOF_CFD_23 [ps]",//17
                                //"TOF_CORR [ps]"//18
                                "Noise_Ampl_1 [mV]",//18
                                "Noise_Ampl_2 [mV]"//19
                                };
    int maximum_bin = 0;
    float maximum = 0.0;
    int nbin[N] = {0};
    float xrange[N][2] = {{1,1},//0
                            // {70,150},//1    DUT1    high gain
                            // {200,300},//8    DUT2    very high gain
                            {25,45},//1    DUT1    low gain
                            // {3,8},//1    DUT1  noise
                            // {150,150},//1    DUT1  laser high gain
                            // {50,50},//1    DUT1  laser low gain
                            {2,2},//2
                            {2.,2.},//3
                            // {25,40},//4 high gain
                            // {45,75},//4 very high gain
                            {4,8},//4   low gain
                            // {50,50},//4   high gain laser

                            // {15,15},//4   low gain laser
                            {100,150},//5
                            // {0.5,0.5},//6
                            {0.5,0.5},//6
                            {60,200},//7

                            {200,300},//8    DUT2    high gain
                            // {15,50},//8    DUT2  low gain
                            {2,2},//9
                            {1.5,2.5},//10
                            {15,25},//11    high gain
                            // {5,15},//11  low gain
                            {100,150},//12
                            {0.5,0.5},//13
                            {60,150},//14
                            // {500,500},//15
                            {1000,2000},//15
                            {1000,2000},//16
                            {500,500},//16
                            // {200,200},//17
                            // {200,200},//18
                            {5,15},//18
                            {5,15}//19
                            };
    int value_range[2] = {0};

    TH1D *hList[N];
    for(i = 0; i < N; i++)
    {
        value_range[0] = -500;
        value_range[1] = 800;

        if(i > 14 && i < 19)
        {
            value_range[0] = -25000;
            value_range[1] = 25000;
        }
        
        // if(i == 17) 
        // {
        //     value_range[0] = 1300;
        //     value_range[1] = 1600;
        // }
        nbin[i] = (value_range[1]-value_range[0])/(xrange[i][0]+xrange[i][1])*40;
        hList[i] = new TH1D(plotname[i], plotname[i], nbin[i], value_range[0], value_range[1]);
    }

    const char *d2_plotname[10] = {"Amax1 vs. Tmax1",//0
                                    "Noise1 vs. Noise2",//1
                                    "Noise vs. Charge",//2
                                    "Tmax1 vs. Tmax2",//3
                                    "Amax2 vs. Tmax2",//4
                                    "Collected charge vs. Risetime",//5
                                    "Amax vs. Delta_max",//6
                                    "Amax1 vs. Amax2",//7
                                    "Amax vs. SNR"//8
                                    };
    const char *d2_axis_name[10][2] = {{"Amax1 [mV]","Tmax1 [ns]"},//0
                                        {"RMS1 [mV]","RMS2 [mV]"},//1
                                        {"RMS [mV]","Collected Charge [fC]"},//2
                                        {"Tmax1 [ns]","Tmax2 [ns]"},//3
                                        {"Amax2 [mV]","Tmax2 [ns]"},//4
                                        {"Collected charge [fC]","Risetime [ps]"},//5
                                        {"Amax [mV]","Delta_max [mV]"},//6
                                        {"Amax1 [mV]","Amax2 [mV]"},//7
                                        {"Amax [mV]","SNR"}//8
                                        };
    int d2_nbin[10][2] = {{100,100},{100,100},{100,100},{100,100},{100,100},{100,100},{100,100},{100,100},{100,100}};
    float d2_xrange[10][2] = {
                                {0,40},//0
                                {0,5},//1
                                {0,1},//2
                                // {-18,-12},//3
                                {-1,1},//3
                                {0,600},//4
                                {0,40},//5
                                {0,100},//6
                                {0,200},//7
                                {0,20}//8
                                };
    float d2_yrange[10][2] = {
                                {-2,2},//0
                                // {-20,20},//0
                                {0,5},//1
                                {0,40},//2
                                // {-18,-12},//3
                                {-1,1},//3
                                {-2,2},//4
                                // {-20,20},//4
                                {0,1000},//5
                                {-2,50},//6
                                {0,800},//7
                                {0,20}//8
                                };

    TH2D *hList_2d[10];
    for(i = 0; i < 9; i++)
    {
        hList_2d[i] = new TH2D(d2_plotname[i],d2_plotname[i],d2_nbin[i][0],d2_xrange[i][0],d2_xrange[i][1],d2_nbin[i][1],d2_yrange[i][0],d2_yrange[i][1]);
    }

    TH2D *hped_amp = new TH2D("ped vs. amp","ped vs. amp",100,-10,20,100,0,50);
    TH2D *hped_noise = new TH2D("ped vs. noise","ped vs. noise",100,0,20,100,0,10);
    TH2D *hped_l_r = new TH2D("ped_l vs. ped_r","ped_l vs. ped_r",100,-10,20,100,-10,20);
    TH2D *hA_ctd = new TH2D("Amax vs. TOF_CTD","Amax vs. TOF_CTD",100,50,300,100,-200,400);

    // a = -2.5;
    // double intercept_range[2] = {0};
    // intercept_range[0] = a;
    // intercept_range[1] = a + 1.4;
    // Fill the tree
    for(n = 0; n < t->GetEntries(); n++)
    {
        // for(i = 0; i < 4; i++)   t[i]->GetEntry(n);
        t->GetEntry(n);
        if(Tmax[1] > 1 || Tmax[1] < -1)   continue;
        // cout<<Amax[1]*1000<<endl;
        // if(noise_amp[1]*1000 < 6.2)  continue;
        // if(noise[1]*1000 > 3)  continue;
        // if(id[1] < 14000)    continue;
        // if(Amax[1]/noise[1] < 5)    continue;
        // if(Amax[2]/noise[2] < 5)    continue;
        // cout<<Amax[2]<<endl;
        // if(Tmax[2]-Tmax[1] < t_cut || Tmax[2]-Tmax[1] > t_cut + t_delta)   continue; //1ns width used
        if (Amax[1]*1000 < threshold_1)   continue;
        if (Amax[2]*1000 < threshold_2)   continue;
        // if (Amax[1]*1000 > 150)   continue;
        if (Amax[2]*1000 > 350)   continue;
        //if(delta_amax < 5)  continue;
        //gr->SetPoint(n,n,ped*1000);
        //if (Amax*1000 > 230)   continue;

        risetime_ave_1 = risetime_ave[1];
        risetime_ave_2 = risetime_ave[2];
        
        // cout<<Amax[1]<<endl;
        double variables[N] = {nentries,//0
                                Amax[1]*1000,//1    DUT1
                                ped[1]*1000,//2
                                noise[1]*1000,//3
                                charge_dw[1]*1000/transimpedance,//4
                                risetime[1]*1000,//5
                                // risetime[1]*1000,//5
                                Amax[1]*0.2/risetime[1],//6
                                1000*noise[1]/(Amax[1]*0.2/risetime_ave_1),//7
                                // 1000*noise[1]/slope[1],//7

                                Amax[2]*1000,//8    DUT2
                                ped[2]*1000,//9
                                noise[2]*1000,//10
                                charge[2]*1000/transimpedance,//11
                                risetime[2]*1000,//12
                                Amax[2]*0.2/risetime[2],//13
                                // 1000*noise[2]/(Amax[2]*0.2/risetime_ave_2),//14
                                1000*noise[2]/slope[2],//14
                                // (CFD[0]-CFD[1])*1000,//15
                                Amax[2]*1000*TOT[2],//15
                                // (CFD[0]-CFD[2])*1000,//16
                                1000*TOT[2],//16
                                (CFD[1]-CFD[2])*1000,//17
                                // (TOA_F[1]-TOA_F[2])*1000//17
                                //TOA_Ref*1000,//12
                                noise_amp[1]*1000,//18
                                noise_amp[2]*1000//19
                                };

        // cout<<variables[13]<<endl;
        event = event + 1;
        for(i = 0; i < N; i++)  hList[i]->Fill(variables[i]);

        // TOF_corr = a*TOT*TOT+b*TOT+c;
        // //TOF_corr = a*Amax*Amax*1e6+b*Amax*1e3+c;
        // hList[14]->Fill(TOF_CTD-TOF_corr);

        hList_2d[0]->Fill(Amax[1]*1000,Tmax[1]);
        hList_2d[1]->Fill(noise[1]*1000,noise[2]*1000);
        hList_2d[2]->Fill(noise[1]*1000,charge[1]*1000/transimpedance);
        hList_2d[3]->Fill(Tmax[1],Tmax[2]);
        hList_2d[4]->Fill(Amax[2]*1000,Tmax[2]);
        hList_2d[5]->Fill(charge[1]*1000/transimpedance,risetime[1]*1000);
        hList_2d[6]->Fill(Amax[1]*1000,delta_amax);
        hList_2d[7]->Fill(Amax[1]*1000,Amax[2]*1000);
        hList_2d[8]->Fill(Amax[1]*1000,Amax[1]/noise[1]);

        //EoE = EoE + fabs(TOF_CFD*1000-TOF_mean);
    }

    //EoE = EoE/event;
    //cout<<"Err_of_Err:"<<EoE<<endl;
    //myfile<<"Err_of_Err"<<EoE<<endl;

    gStyle->SetOptStat(1111);
    
    // TF1 *f1 = new TF1("f1","[0]+[1]*x");
    // f1->SetLineColor(kRed);
    // r = hList_2d[3]->Fit(f1,"QS","C");
    // a = r->Value(0);
    // b = r->Value(1);

    for(i = 0; i < 9; i++)
    {
        // Draw Amax vs. Tmax distribution
        hList_2d[i]->GetXaxis()->SetTitle(d2_axis_name[i][0]);
        hList_2d[i]->GetYaxis()->SetTitle(d2_axis_name[i][1]);
        hList_2d[i]->GetXaxis()->SetNdivisions(505, kTRUE);
        hList_2d[i]->GetYaxis()->SetNdivisions(505, kTRUE);
        hList_2d[i]->SetMarkerSize(0.2);
        hList_2d[i]->Draw("COLZ");
        
        if(i == 3)
        {
            TF1 *f2 = new TF1("f2","[0]+[1]*x",-25,25);
            f2->SetLineColor(kRed);
            f2->SetLineWidth(2);
            f2->SetParameter(0,t_cut);
            f2->SetParameter(1,1);
            f2->Draw("same");

            TF1 *f3 = new TF1("f3","[0]+[1]*x",-25,25);
            f3->SetLineColor(kRed);
            f3->SetLineWidth(2);
            f3->SetParameter(0,t_cut + t_delta);
            f3->SetParameter(1,1);
            f3->Draw("same");
            
            pt=new TPaveText(0.7,0.75,0.3,0.6,"NDC");
            pt->SetBorderSize(0);
            pt->SetFillColor(0);
            pt->SetFillStyle(0);
            pt->SetTextAlign(12);
            pt->SetTextSize(0.03);

            text=pt->AddText(Form("intercept: %.4f",a));
            text=pt->AddText(Form("slope: %.4f",b));

            // pt->Draw();
        }

        c1->SaveAs(filePath + dirname + "/" + d2_plotname[i] + "_Dist.png");
        c1->SaveAs(filePath + dirname + "/" + d2_plotname[i] + "_Dist.svg");
        c1->Clear();
    }

    // gStyle->SetOptStat(111);
    // cout<<TMath::MaxElement(hList[4]->GetN(),hList[4]->GetX())<<endl;

    double fr[2], sv[4], fp[4], fpe[4];
    // sv[0]=1.8; sv[1]=100.0; sv[2]=50000.0; sv[3]=3.0;

    for(i = 0; i < N; i++)
    {    
        // cout<<plotname[i]<<endl;
        myfile<<plotname[i]<<endl;

        //hpeak->GetXaxis()->SetRangeUser(-1,1);
        hList[i]->GetXaxis()->SetTitle(xaxis_name[i]);
        hList[i]->GetYaxis()->SetTitle("Entries");
        
        maximum_bin = hList[i]->GetMaximumBin();
        maximum = hList[i]->GetXaxis()->GetBinCenter(maximum_bin);
        hList[i]->GetXaxis()->SetRangeUser(maximum - xrange[i][0], maximum + xrange[i][1]);
        // if(i == 1 || i == 4 || i == 8 || i == 11)   hList[i]->GetXaxis()->SetRangeUser(maximum/3, maximum*3);
        hList[i]->GetXaxis()->SetNdivisions(505, kTRUE);
        //cout<<maximum - xrange[i][0]<<","<<maximum + xrange[i][1]<<endl;
        //hList[i]->SetStats(1);
        hList[i]->Draw();

        if(i == 0)
        {
            pt=new TPaveText(1.2,0.75,0.6,0.6,"NDC");
            pt->SetBorderSize(0);
            pt->SetFillColor(0);
            pt->SetFillStyle(0);
            pt->SetTextAlign(12);
            pt->SetTextSize(0.03);

            text=pt->AddText(Form("Entries: %d",n));

            pt->Draw();
        }

        //amplitude, charge, tot 
        else if(i == 1 || i==3 ||i == 4 || i == 7 || i == 8 || i == 11 || i == 14 || i == 15 || i == 16 || i == 18 || i == 19)
        // else if(i == 8 || i == 15)//For Laser
        {
            // r = hList[i]->Fit("landau", "S", "C", maximum/3, maximum*3);
            // r = hList[i]->Fit("landau", "S", "C", maximum - xrange[i][0], maximum + xrange[i][1]);
            // MPV=r->Value(1);
            // MPV_err=r->ParError(1);
            // sigma=r->Value(2);

            fr[0] = maximum - xrange[i][0];
            fr[1] = maximum + xrange[i][1];
            sv[0] = hList[i]->GetStdDev()/100;
            sv[1] = hList[i]->GetMean();
            sv[2] = 2*(fr[1] - fr[0]);
            sv[3] = sv[0];
            // cout<<i<<","<<sv[1]<<endl;
            
            fit_langau = langaufit(hList[i],fr,sv,fp,fpe,&chisqr,&ndf);
            fit_langau->Draw("lsame");
            langaupro(fp,langau_Peak,langau_FWHM);

            // cout<<langau_FWHM<<endl;
            // cout<<sigma<<endl;
            
            fitting[i] = langau_Peak;
            if(i != 4)  myfile<<"MPV:"<<langau_Peak<<",MPV_ERR_fitting:"<<0<<",MPV_ERR:"<<langau_FWHM/sqrt(event)<<",Sigma:"<<langau_FWHM<<",Entries:"<<event<<endl;
            // myfile<<"MPV:"<<MPV<<",MPV_ERR_fitting:"<<MPV_err<<",MPV_ERR:"<<sigma/sqrt(event)<<",Sigma:"<<sigma<<",Entries:"<<event<<endl;

            else 
            {
                cout<<"MPV:"<<fp[1]<<",MPV ERR:"<<fpe[1]<<endl;;
                stat_err = langau_FWHM/sqrt(event);
                cali_err = langau_Peak*0.00294/0.2131;

                myfile<<"MPV:"<<langau_Peak<<",MPV_ERR_fitting:"<<0<<",MPV_ERR:"<<sqrt(stat_err*stat_err+cali_err*cali_err)<<",Sigma:"<<langau_FWHM<<",Entries:"<<event<<endl;

                cout<<"charge:stat_err:"<<stat_err<<",cali_err:"<<cali_err<<endl;
            }
    
            pt=new TPaveText(1.2,0.75,0.6,0.6,"NDC");
            pt->SetBorderSize(0);
            pt->SetFillColor(0);
            pt->SetFillStyle(0);
            pt->SetTextAlign(12);
            pt->SetTextSize(0.03);

            text=pt->AddText(Form("MPV: %.2f",langau_Peak));
            // text=pt->AddText(Form("MPV_Err: %.2f",MPV_err));
            text=pt->AddText(Form("Sigma: %.2f",langau_FWHM));
            text=pt->AddText(Form("Entries: %d",event));

            pt->Draw();
        }

        else
        {                   
            r = hList[i]->Fit("gaus", "QS", "C", maximum - xrange[i][0], maximum + xrange[i][1]);
            mean=r->Value(1);
            mean_err=r->ParError(1);
            sigma=r->Value(2);
            sigma_err=r->ParError(2);

            if(i == 15) sigma_12 = sigma;
            if(i == 16) sigma_13 = sigma;
            if(i == 17)
            {
                sigma_23 = sigma;
                sigma_23_err = sigma_err;
            }

            fitting[i] = mean;
            myfile<<"Mean:"<<mean<<",Mean_ERR_fitting:"<<mean_err<<",Mean_ERR:"<<sqrt(mean_err*mean_err+sigma*sigma/event)<<",Sigma:"<<sigma<<",Sigma_ERR:"<<sigma_err<<",Entries:"<<event<<endl;
    
            pt=new TPaveText(1.2,0.75,0.6,0.6,"NDC");
            pt->SetBorderSize(0);
            pt->SetFillColor(0);
            pt->SetFillStyle(0);
            pt->SetTextAlign(12);
            pt->SetTextSize(0.03);

            text=pt->AddText(Form("Mean: %.2f",mean));
            text=pt->AddText(Form("Mean_Err: %.2f",mean_err));
            text=pt->AddText(Form("Sigma: %.2f",sigma));
            text=pt->AddText(Form("Sigma_Err: %.2f",sigma_err));
            text=pt->AddText(Form("Entries: %d",event));

            pt->Draw();
        }

        c1->SaveAs(filePath + dirname + "/" + plotname[i] + "_Dist.png");
        c1->SaveAs(filePath + dirname + "/" + plotname[i] + "_Dist.svg");
        c1->Clear(); 
    }

    sigma_1 = sqrt((sigma_12*sigma_12+sigma_13*sigma_13-sigma_23*sigma_23)/2);
    // sigma_2 = sqrt((sigma_12*sigma_12+sigma_23*sigma_23-sigma_13*sigma_13)/2);
    // cout<<sigma_23<<endl;
    // cout<<res_ref<<endl;
    sigma_2 = sqrt(sigma_23*sigma_23-res_ref*res_ref); //  HPK T3.1 single ref
    // sigma_2 = sqrt(sigma_23*sigma_23-27.29*27.29);    //  HPK T1.1 single ref
    // sigma_2 = sigma_23/sqrt(2);    //  same ref
    sigma_3 = sqrt((sigma_23*sigma_23+sigma_13*sigma_13-sigma_12*sigma_12)/2);

    stat_err = sigma_23/sigma_2*sigma_23_err;
    ref_err = res_ref/sigma_2*res_ref_err;

    jitter_ave_1 = 1000*fitting[3]/(fitting[1]*0.8/risetime_ave_1);
    jitter_ave_2 = 1000*fitting[10]/(fitting[8]*0.8/risetime_ave_2);

    cout<<"charge_1:"<<fitting[4]<<",charge_2:"<<fitting[11]<<endl;

    cout<<"sigma_1:"<<sigma_1<<",sigma_2:"<<sigma_2<<",sigma_3:"<<sigma_3<<endl;

    myfile<<"SNR_1"<<endl<<fitting[1]/fitting[3]<<endl;
    myfile<<"SNR_2"<<endl<<fitting[8]/fitting[10]<<endl;

    myfile<<"sigma_1"<<endl<<sigma_1<<endl;
    // myfile<<"sigma_2"<<endl<<"value:"<<sigma_2<<",error:"<<sqrt(2)*sigma_23_err<<endl;    
    myfile<<"sigma_2"<<endl<<"value:"<<sigma_2<<",error:"<<sqrt(stat_err*stat_err+ref_err*ref_err)<<",stat_err:"<<stat_err<<",ref_err:"<<ref_err<<endl;
    myfile<<"sigma_3"<<endl<<"value:"<<sigma_3<<",error:"<<sigma_23/sigma_3*sigma_23_err<<endl;

    cout<<"TRes: stat_err:"<<stat_err<<",ref_err:"<<ref_err<<endl;

    myfile<<"risetime_ave_1"<<endl<<risetime_ave_1<<endl;
    myfile<<"risetime_ave_2"<<endl<<risetime_ave_2<<endl;
    myfile<<"slope_ave_1"<<endl<<fitting[1]*0.2/risetime_ave_1<<endl;
    myfile<<"slope_ave_2"<<endl<<fitting[8]*0.2/risetime_ave_2<<endl;
    myfile<<"jitter_ave_1"<<endl<<jitter_ave_1<<endl;
    // myfile<<fitting[7]<<endl;
    myfile<<"jitter_ave_2"<<endl<<jitter_ave_2<<endl;
    myfile<<"landau_1"<<endl<<sqrt(sigma_2*sigma_2-jitter_ave_1*jitter_ave_1)<<endl;
    myfile<<"landau_2"<<endl<<sqrt(sigma_3*sigma_3-jitter_ave_2*jitter_ave_2)<<endl;

    cout<<"total number:"<<n<<",events number:"<<event<<",background rejection:"<<1-event*1.0/n<<endl;
}

int main(int argc, char **argv)
{
    if (argc == 1 || argc > 2)
        cout << "Usage: " << argv[0] << " filePath";
    else
    {
        TString filePath(argv[1]);
        if (!filePath.EndsWith(".root"))
            filePath = filePath + ".root";
        // string dataPath (argv[1]);
        cout << filePath << endl;
        file_read(filePath,threshold_1);
    }
    return 0;
}