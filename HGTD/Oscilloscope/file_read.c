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
// double transimpedance = 0.4;
// double transimpedance = 4.3226;//UCSCv1.1
// double transimpedance = 4.6929;//UCSCv1.4
// double transimpedance = 16.262;//USTCv4 B5 cold
double transimpedance = 18.19;//USTCv4 B1
// double transimpedance = 26.42;//USTCv2
// double res_ref = 37.7;//T3.1
// double res_ref_err = 0.83;
const int N = 34;

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
    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 900, 800);
    
    //TText *t = new TText();
    // c1->SetLeftMargin(0.25);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.05);

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
    myfile.open(filePath + dirname + "/results.txt");
    myfile<<filePath + dirname<<endl;

    // ofstream temp;
    // temp.open("tmp.csv");
    // temp<<"Amplitude,Charge,RMS,delta_TOA"<<endl;

    double res_ref = 30.;//T1.1
    double res_ref_err = 1.0;

    if(filePath.Contains("/20/"))   //20 degree
    {    
        res_ref = 29.34;//T1.1
        res_ref_err = 0.83;
    }
    else if(filePath.Contains("/M30/")) //-30 degree
    {
        res_ref = 33.93;//T1.1
        res_ref_err = 0.81;
    }
    
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
    Float_t Amax[4], Tmax[4], noise[4], tail_noise[4], noise_amp[i], charge[4], charge_hw[4], charge_dw[4], CTD[4], CFD[4], ZCD[4], TOA_F[4], TOT[4], nentries, ped[4], risetime_ave[4], risetime[4], risetime_F[4], slope[4], linearity, delta_amax, jitter[4], ped_1[4], ped_2[4], slope_F[4], intercept_F[4], id[4];
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
    
    t->SetBranchAddress("Amax",Amax);   
    t->SetBranchAddress("Tmax",Tmax);
    t->SetBranchAddress("pedestal",ped);
    // t->SetBranchAddress("ped_1",ped_1);
    // t->SetBranchAddress("ped_2",ped_2);
    t->SetBranchAddress("Noise",noise);
    t->SetBranchAddress("Tail_Noise",tail_noise);
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
    const char *plotname[N] = {"nentries",      //0

                                "amplitude_1",  //1   Ch2
                                "pedstal_1",    //2
                                "RMS_1",        //3
                                "charge_1",     //4
                                "risetime_1",   //5
                                "slope_1",      //6
                                "jitter_1",     //7

                                "amplitude_2",  //8   Ch3
                                "pedstal_2",    //9
                                "RMS_2",        //10
                                "charge_2",     //11
                                "risetime_2",   //12
                                "slope_2",      //13
                                "jitter_2",     //14
                                
                                "amplitude_3",  //15   Ch0
                                "pedstal_3",    //16
                                "RMS_3",        //17
                                "charge_3",     //18
                                "risetime_3",   //19
                                "slope_3",      //20
                                "jitter_3",     //21
                                
                                "amplitude_4",  //22   Ch4
                                "pedstal_4",    //23
                                "RMS_4",        //24
                                "charge_4",     //25
                                "risetime_4",   //26
                                "slope_4",      //27
                                "jitter_4",     //28

                                "tof_cfd_12",   //29
                                "tof_cfd_13",   //30
                                "tof_cfd_23",   //31
                                //"tof_corr"    //32
                                "noise_amplitude_1",//33
                                "noise_amplitude_2"//34
                                };

    const char *xaxis_name[N] = {"nentries",//0
                                // "Ampl_DUT [mV]",            //1   Ch2
                                // "Pedestal_DUT [mV]",        //2
                                // "Noise RMS_DUT [mV]",       //3
                                // "Collected charge_DUT [fC]",//4
                                // "risetime_DUT [ps]",        //5
                                // "slope_DUT [mV/ps]",        //6
                                // "jitter_DUT [ps]",          //7
                                // "Ampl [mV]",            //1   Ch2
                                "Amplitude [mV]",            //1   Ch2
                                "Pedestal [mV]",        //2
                                "RMS Noise [mV]",       //3
                                "Collected charge [fC]",//4
                                "risetime [ps]",        //5
                                "slope [mV/ps]",        //6
                                "jitter [ps]",          //7

                                "Ampl_Ref [mV]",            //8   Ch3
                                "Pedestal_Ref [mV]",        //9
                                "Noise RMS_Ref [mV]",       //10
                                "Collected charge_Ref [fC]",//11
                                "risetime_Ref [ps]",        //12
                                "slope_Ref [mV/ps]",        //13
                                "jitter_Ref [ps]",          //14
                                
                                "Ampl_DUT_2 [mV]",            //15   Ch0
                                "Pedestal_DUT_2 [mV]",        //16
                                "Noise RMS_DUT_2 [mV]",       //17
                                "Collected charge_DUT_2 [fC]",//18
                                "risetime_DUT_2 [ps]",        //19
                                "slope_DUT_2 [mV/ps]",        //20
                                "jitter_DUT_2 [ps]",          //21
                                
                                "Ampl_DUT_3 [mV]",            //22  Ch4
                                "Pedestal_DUT_3 [mV]",        //23
                                "Noise RMS_DUT_3 [mV]",       //24
                                "Collected charge_DUT_3 [fC]",//25
                                "risetime_DUT_3 [ps]",        //26
                                "slope_DUT_3 [mV/ps]",        //27
                                "jitter_DUT_3 [ps]",          //28

                                // "TOF_CFD_12 [ps]",//15
                                "Amplitude*TOT [mV*ns]",//29
                                "TOT [ps]",//30
                                // "TOF_CFD_23 [ps]",//17
                                "#DeltaTOA [ps]",//31
                                //"TOF_CORR [ps]"//18
                                "Noise_Ampl_DUT [mV]",//32
                                "Noise_Ampl_Ref [mV]"//33
                                };
    int maximum_bin = 0;
    float maximum = 0.0;
    int nbin[N] = {0};
    float xrange[N][2] = {{1,1},//0
                            // {80,150},       //1    Ch2    high gain
                            {280,510},      //1    Ch2    high gain
                            // {30,45},     //1    Ch2    low gain
                            // {3,5},       //1    Ch2  noise
                            // {150,150},   //1    Ch2  laser high gain
                            // {20,20},     //1    Ch2  laser low gain
                            {2,2},          //2
                            // {1.,1.},     //3
                            {3.,3.},      //3 high noise
                            {15,35},     //4 high gain
                            // {45,75},        //4 very high gain
                            // {10,25},      //4   low gain
                            // {500,500},   //4   high gain laser 
                            // {10,10},        //4   low gain laser
                            {100,150},      //5
                            // {0.5,0.5},   //6
                            {0.5,0.5},      //6
                            {5,8},          //7
                            // {10,20},     //7

                            {280,510},      //8    Ch3    high gain
                            // {80,150},    //8    Ch3  low gain
                            // {150,150},        //8    Ch3  laser low gain
                            {2,2},          //9
                            // {1.,1.},     //10
                            {3.,3.},      //10 high noise
                            {15,35},        //11    high gain
                            // {5,15},      //11  low gain
                            // {10,10},      //11  low gain
                            {100,150},      //12
                            {0.5,0.5},      //13
                            {20,30},        //14

                            // {250,400},      //15   Ch0  high gain
                            // {60,140},    //15   Ch0  low gain
                            {150,150},        //15    Ch0  laser low gain
                            {2,2},          //16
                            // {1.,1.},     //10
                            {2.5,2.5},      //17 high noise
                            // {15,25},        //18    high gain
                            // {5,15},      //18  low gain
                            {10,10},      //18  low gain
                            {100,150},      //19
                            {0.5,0.5},      //20
                            {20,30},        //21

                            // {250,400},      //22    Ch4    high gain
                            // {60,140},    //22    Ch4  low gain
                            {200,200},        //22    Ch4  laser low gain
                            {2,2},          //23
                            // {1.,1.},     //10
                            {2.5,2.5},      //24 high noise
                            {25,25},        //25    high gain
                            // {5,15},      //25  low gain
                            {100,150},      //26
                            {0.5,0.5},      //27
                            {20,30},        //28

                            {200,200},      //29
                            // {1000,2000}, //25
                            // {1500,3500}, //26
                            {500,500},      //30
                            // {100,100},   //27
                            {500,500},      //31
                            // {5,15},      //28
                            {3,5},          //32
                            {5,15}          //33
                            };
    int value_range[2] = {0};

    TH1D *hList[N];
    for(i = 0; i < N; i++)
    {
        // value_range[0] = -500;
        // value_range[1] = 500;
        value_range[0] = -3500;
        value_range[1] = 3500;

        if(i > 28 && i < 33)    //  Timing measurement
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

    const char *d2_plotname[10] = {"Amax1_vs_Tmax1",//0
                                    "Noise1_vs_Noise2",//1
                                    "Head_Noise_vs_Tail_Noise",//2
                                    "Tmax1_vs_Tmax2",//3
                                    "Amax2_vs_Tmax2",//4
                                    "Collected charge_vs_Risetime",//5
                                    "Amax_vs_Delta_max",//6
                                    "Amax1_vs_Amax2",//7
                                    "Amax_vs_SNR"//8
                                    };
    const char *d2_axis_name[10][2] = {{"Amax_DUT [mV]","Tmax_DUT [ns]"},//0
                                        {"RMS1 [mV]","RMS2 [mV]"},//1
                                        {"Head_RMS [mV]","Tail_RMS [mV]"},//2
                                        {"Tmax_DUT [ns]","Tmax_Ref [ns]"},//3
                                        {"Amax2 [mV]","Tmax2 [ns]"},//4
                                        {"Collected charge [fC]","Risetime [ps]"},//5
                                        {"Amax [mV]","Delta_max [mV]"},//6
                                        {"Amax1 [mV]","Amax2 [mV]"},//7
                                        {"Amax [mV]","SNR"}//8
                                        };
    int d2_nbin[10][2] = {{100,100},{100,100},{100,100},{100,100},{100,100},{100,100},{100,100},{100,100},{100,100}};
    float d2_xrange[10][2] = {
                                {0,800},//0
                                {0,5},//1
                                {0,100.},//2
                                // {-18,-12},//3
                                {-2,2},//3
                                {0,1000},//4
                                {0,500},//5
                                {0,100},//6
                                {0,2000},//7
                                {0,20}//8
                                };
    float d2_yrange[10][2] = {
                                {-2,2},//0
                                // {-25,25},//0
                                // {-20,20},//0
                                {0,5},//1
                                {0,100.},//2
                                // {-18,-12},//3
                                {-2,2},//3
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

    TH2D *hped_amp = new TH2D("ped_vs_amp","ped_vs_amp",100,-10,20,100,0,50);
    TH2D *hped_noise = new TH2D("ped_vs_noise","ped_vs_noise",100,0,20,100,0,10);
    TH2D *hped_l_r = new TH2D("ped_l_vs_ped_r","ped_l_vs_ped_r",100,-10,20,100,-10,20);
    TH2D *hA_ctd = new TH2D("Amax_vs_TOF_CTD","Amax_vs_TOF_CTD",100,50,300,100,-200,400);

    // a = -2.5;
    // double intercept_range[2] = {0};
    // intercept_range[0] = a;
    // intercept_range[1] = a + 1.4;
    // Fill the tree
    for(n = 0; n < t->GetEntries(); n++)
    {
        // for(i = 0; i < 4; i++)   t[i]->GetEntry(n);
        t->GetEntry(n);
        // if(Tmax[1] > -15 || Tmax[1] < -100)   continue;
        // cout<<Amax[1]*1000<<endl;
        // if(noise_amp[1]*1000 < threshold_1)  continue;
        // if(noise[1]*1000 > 3)  continue;
        if(tail_noise[1]*1000 > 20)  continue;
        // if(id[1] < 14000)    continue;
        // if(Amax[1]/noise[1] < 5)    continue;
        // if(Amax[2]/noise[2] < 5)    continue;
        // cout<<Amax[2]<<endl;
        // if(Tmax[2]-Tmax[1] < t_cut || Tmax[2]-Tmax[1] > t_cut + t_delta)   continue; //1ns width used
        if (Amax[1]*1000 < threshold_1)   continue;
        if (Amax[2]*1000 < threshold_2)   continue;
        if (Amax[1]*1000 > 900)   continue;
        if (Amax[2]*1000 > 900)   continue;
        // if (Amax[0]*1000 > 1000)   continue;
        // if (Amax[3]*1000 > 1000)   continue;
        //if(delta_amax < 5)  continue;
        // if ((CFD[1]-CFD[2])*1000 < 100)   continue;
        // if (1000*TOT[1] > 4500)   continue;

        // risetime_ave_1 = risetime_ave[1];
        // risetime_ave_2 = risetime_ave[2];
        
        // temp<<Amax[1]*1000<<","<<charge[1]*1000/transimpedance<<","<<noise[1]*1000<<","<<(CFD[1]-CFD[2])*1000<<endl;

        // cout<<Amax[1]<<endl;
        double variables[N] = {nentries,//0
                                Amax[1]*1000,//1    Ch2
                                ped[1]*1000,//2
                                noise[1]*1000,//3
                                charge[1]*1000/transimpedance,//4
                                risetime[1]*1000,//5
                                // risetime[1]*1000,//5
                                Amax[1]*0.2/risetime[1],//6
                                1000*noise[1]/(Amax[1]*0.2/risetime_ave[1]),//7
                                // 1000*noise[1]/slope[1],//7

                                Amax[2]*1000,//8    Ch3
                                ped[2]*1000,//9
                                noise[2]*1000,//10
                                charge[2]*1000/transimpedance,//11
                                risetime[2]*1000,//12
                                Amax[2]*0.2/risetime[2],//13
                                1000*noise[2]/(Amax[2]*0.2/risetime_ave[2]),//14
                                
                                Amax[0]*1000,//15    Ch0
                                ped[0]*1000,//16
                                noise[0]*1000,//17
                                charge_dw[0]*1000/transimpedance,//18
                                risetime[0]*1000,//19
                                // risetime[1]*1000,//5
                                Amax[0]*0.2/risetime[0],//20
                                1000*noise[0]/(Amax[0]*0.2/risetime_ave[0]),//21
                                
                                Amax[3]*1000,//22    Ch4
                                ped[3]*1000,//23
                                noise[3]*1000,//24
                                charge_dw[3]*1000/transimpedance,//25
                                risetime[3]*1000,//26
                                // risetime[1]*1000,//5
                                Amax[3]*0.2/risetime[3],//27
                                1000*noise[3]/(Amax[3]*0.2/risetime_ave[3]),//28
                                // 1000*noise[1]/slope[1],//7

                                // 1000*noise[2]/slope[2],//14
                                // (CFD[0]-CFD[1])*1000,//15
                                Amax[2]*1000*TOT[2],//29
                                // (CFD[0]-CFD[2])*1000,//16
                                1000*TOT[1],//30
                                // CFD[3]*1000,//17
                                (CFD[1]-CFD[2])*1000,//31
                                // (TOA_F[1]-TOA_F[2])*1000//17
                                //TOA_Ref*1000,//12
                                noise_amp[1]*1000,//32
                                noise_amp[2]*1000//33
                                };

        // cout<<variables[13]<<endl;
        event = event + 1;
        for(i = 0; i < N; i++)  hList[i]->Fill(variables[i]);

        // TOF_corr = a*TOT*TOT+b*TOT+c;
        // //TOF_corr = a*Amax*Amax*1e6+b*Amax*1e3+c;
        // hList[14]->Fill(TOF_CTD-TOF_corr);

        hList_2d[0]->Fill(Amax[1]*1000,Tmax[1]);
        // hList_2d[0]->Fill(Amax[1]*1000,TOT[1]);
        hList_2d[1]->Fill(noise[1]*1000,noise[2]*1000);
        hList_2d[2]->Fill(noise[1]*1000,tail_noise[1]*1000);
        hList_2d[3]->Fill(Tmax[1],Tmax[2]);
        hList_2d[4]->Fill(Amax[2]*1000,Tmax[2]);
        hList_2d[5]->Fill(charge[1]*1000/transimpedance,risetime[1]*1000);
        hList_2d[6]->Fill(Amax[1]*1000,delta_amax);
        hList_2d[7]->Fill(Amax[1]*1000,Amax[0]*1000);
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
        // Draw Amax vs Tmax distribution
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
            
            pt=new TPaveText(0.8,0.85,0.4,0.7,"NDC");
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
        // c1->SaveAs(filePath + dirname + "/" + d2_plotname[i] + "_Dist.svg");
        // c1->SaveAs(filePath + dirname + "/" + d2_plotname[i] + "_Dist.pdf");
        c1->Clear();
    }

    // gStyle->SetOptStat(111);
    // cout<<TMath::MaxElement(hList[4]->GetN(),hList[4]->GetX())<<endl;

    double fr[2], sv[4], fp[4], fpe[4];
    // sv[0]=1.8; sv[1]=100.0; sv[2]=50000.0; sv[3]=3.0;

    for(i = 0; i < N; i++)
    {    
        // cout<<plotname[i]<<endl;
        if(!hList[i]->Integral())   hList[i]->SetBinContent(1,1);
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

        pt=new TPaveText(1.25,0.9,0.65,0.65,"NDC");
        pt->SetBorderSize(0);
        pt->SetFillColor(0);
        pt->SetFillStyle(0);
        pt->SetTextAlign(12);
        pt->SetTextSize(0.045);

        // if(i == 0 || i == 18)
        if(i == 0)
        {
            if(i==0)    text=pt->AddText(Form("Entries: %d",n));
            else    text=pt->AddText(Form("Entries: %d",event));
        }

        //amplitude, charge, tot 
        else if(i == 1 || i == 4 || i == 7 || 
        i == 8 || i == 11 || i == 14 ||
        i == 15 || i == 18 || i == 21 ||
        i == 22 || i == 25 || i == 28||
        i == 12 || i == 33 || i == 34)
        // else if(i == 29)//For Laser
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

            text=pt->AddText(Form("MPV: %.2f",langau_Peak));
            // text=pt->AddText(Form("MPV_Err: %.2f",MPV_err));
            text=pt->AddText(Form("FWHM: %.2f",langau_FWHM));
            text=pt->AddText(Form("Entries: %d",event));
        }

        else
        {                   
            r = hList[i]->Fit("gaus", "QS", "C", maximum - xrange[i][0], maximum + xrange[i][1]);
            mean=r->Value(1);
            mean_err=r->ParError(1);
            sigma=r->Value(2);
            sigma_err=r->ParError(2);

            if(i == 29) sigma_12 = sigma;
            if(i == 30) sigma_13 = sigma;
            if(i == 31)
            {
                sigma_23 = sigma;
                sigma_23_err = sigma_err;
            }

            fitting[i] = mean;
            // myfile<<"Mean:"<<mean<<",Mean_ERR_fitting:"<<mean_err<<",Mean_ERR:"<<sqrt(mean_err*mean_err+sigma*sigma/event)<<",Sigma:"<<sigma<<",Sigma_ERR:"<<sigma_err<<",Entries:"<<event<<endl;
                
            myfile<<"Mean:"<<hList[i]->GetMean()<<",Mean_ERR_fitting:"<<0<<",Mean_ERR:"<<0<<",Sigma:"<<hList[i]->GetRMS()<<",Entries:"<<event<<endl;

            text=pt->AddText(Form("Mean: %.2f",mean));
            text=pt->AddText(Form("Mean_Err: %.2f",mean_err));
            text=pt->AddText(Form("Sigma: %.2f",sigma));
            text=pt->AddText(Form("Sigma_Err: %.2f",sigma_err));
            text=pt->AddText(Form("Entries: %d",event));

        }

        pt->Draw();
        c1->SaveAs(filePath + dirname + "/" + plotname[i] + "_Dist.png");
        c1->SaveAs(filePath + dirname + "/" + plotname[i] + "_Dist.svg");
        // c1->SaveAs(filePath + dirname + "/" + plotname[i] + "_Dist.pdf");
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

    jitter_ave_1 = 1000*fitting[3]/(fitting[1]*0.2/risetime_ave[1]);
    jitter_ave_2 = 1000*fitting[10]/(fitting[8]*0.2/risetime_ave[2]);

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
    myfile<<"slope_ave_1"<<endl<<fitting[1]*0.2/risetime_ave[1]<<endl;
    myfile<<"slope_ave_2"<<endl<<fitting[8]*0.2/risetime_ave[2]<<endl;
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