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
void FitPedestal(vector<double> xlist, vector<double> ylist, int length, double *a, double *b, double *c, int j);
void PedestalAndNoise(vector<double> ylist, double *ped, double *noise, double *noise_amp);
double Peak(vector<double> xlist, vector<double> Alist,double *Tmax, double *Amax, int *Max, double threshold, double init_f, double fin_f);
double leftTime(vector<double> xlist, vector<double> Alist, double threshold, int Max);
double rightTime(vector<double> xlist, vector<double> Alist, double threshold, int Max);
double Charge(vector<double> xlist, vector<double> Alist, int Max, int left, int right, double threshold);
double Time, Ampl;
double transImpedance = 1.0;
bool average_only = false;
bool multi = false;
bool draw = false;
vector<int> colorList{ kRed, kBlue, kGreen, kViolet, kMagenta, kTeal, kPink, kAzure, kCyan, kYellow,
                    kBlue + 1, kRed + 1, kGreen + 1, kTeal + 1, kPink + 1, kMagenta + 1, kViolet + 1, kAzure + 1, kCyan + 1, kYellow + 1};

void oscilloscope_draw(TString dataPath, bool average_only = false, bool draw = false, bool multi = false)
{
    gErrorIgnoreLevel = kError;

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(.03);

    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 1000, 800);
    //TText *t = new TText();
    c1->SetLeftMargin(0.25);
    c1->SetRightMargin(0.2);
    c1->SetTopMargin(0.2);

    TString plotType = dataPath;
    TString plotname = dataPath;
    plotname.ReplaceAll("/","");
    dataPath = "Data/" + dataPath;
    void *plotDir = gSystem->OpenDirectory(dataPath);
    
    // Check Plot directory exist or not, if not, create it
    if (!gSystem->OpenDirectory(plotStorePath + plotType))
        gSystem->MakeDirectory(plotStorePath + plotType);

    //cout<<dataPath<<endl;
    TString txtFile;
    const char *entry, *entry_temp;

    //Some parameters
    double threshold = 0.0;
    double trig_th = 0.02, th = 0.04, fra = 0.5, ZCD_fra = 0.5, shift = 1;
    //if(plotType.Contains("Board1_Ch2")) threshold = 0.02;
    
    int i = 0, j = 0, nentries = 0, dual_board = 0, ch_num = 0, fit_ped = 0, doZCD = 0;
    int N = 20, Ns = 500;
    int n = 0, length = 9999, Max[4] = {0};
    double a, b, c, d;
    double sum[4][9999] = {0}, t[4][9999] = {0};
    double ped[4] = {0}, noise[4] = {0}, noise_amp[4] = {0}, charge[4] = {0}, charge_hw[4] = {0}, charge_dw[4] = {0}, Tmax[4] = {0}, Amax[4] = {0}, CTD[4] = {0}, CFD[4] = {0}, ZCD[4] = {0}, TOA[4] = {0}, TOT[4] = {0}, rise_start = 0, rise_end = 0, risetime_ave[4] = {0}, risetime[4] = {0}, slope[4] = {0}, linearity[4] = {0}, delta_amax[4] = {0}, ped_1[4] = {0}, ped_2[4] = {0}, ave_risetime[4] = {0}, noise_1[4], noise_2[4], noise_amp_1[4], noise_amp_2[4], id[4] = {0};
    //double ped1 = 0, noise1 = 0, noise_amp1 = 0, charge1 = 0, Tmax1 = 0, Amax1 = 0, CTD1 = 0, CFD1 = 0, ZCD1 = 0, TOT1 = 0, risetime1 = 0, slope1 = 0;
    //double TOF_CTD = 0, TOF_CFD = 0, TOF_ZCD = 0;
    double TOA_Ref = 0, CTD_C = 0, init_f = 0.1, fin_f = 0.9;
    double simple_rise_start, simple_rise_end, simple_risetime;

    vector<TString> dirList, list_0, list_1, list_2;
    while ((entry = (char *)gSystem->GetDirEntry(plotDir)))
    {
        if(!TString(entry).Contains(".txt"))    continue;
        if(TString(entry).Contains("C1"))   list_0.push_back(entry);
        if(TString(entry).Contains("C2"))   list_1.push_back(entry);
        if(TString(entry).Contains("C3"))   list_2.push_back(entry);
        dirList.push_back(entry);
    } 
    sort(dirList.begin(), dirList.end());

    // Average waveform
    TString txtFileFullPath;
    vector<int>  t_ref,id_ref;
    int  ch_dut = 1, ch_ref = 2;
    fin_f = 0.4;
    
    // Reference processing
    for(j = 0; j < list_2.size(); j++)
    {
        //if(j > 100)  break;
        //cout<<j<<endl;
        if(list_2.size() < 1)   break;
        txtFile = list_2[j];
        vector<double>  xlist_ave, ylist_ave;

        txtFileFullPath = dataPath + txtFile;
        //cout<<txtFileFullPath<<endl;
        TTree *tree = new TTree("waveforms", "waveforms of oscillscope");
        tree->ReadFile(txtFileFullPath, "Time/D:Ampl/D", ',');
        //tree->ReadFile(txtFileFullPath, "Time/D:Ampl/D", ' ');
        
        // Start To read file into TTree
        tree->SetBranchAddress("Time", &Time);
        tree->SetBranchAddress("Ampl", &Ampl);

        for (n = 0; n < tree->GetEntries(); n++)
        {
            tree->GetEntry(n);
            xlist_ave.push_back(Time*1e9);
            // ylist_ave.push_back(Ampl);
            ylist_ave.push_back(-Ampl);
        }
        
        tree->Delete("");
        
        if(length > xlist_ave.size())   length = xlist_ave.size();

        PedestalAndNoise(ylist_ave,&ped[ch_ref],&noise[ch_ref],&noise_amp[ch_ref]);
        Peak(xlist_ave,ylist_ave,&Tmax[ch_ref],&Amax[ch_ref],&Max[ch_ref],threshold,init_f,fin_f);

        if(Amax[ch_ref]*1000 < 50)  continue;

        n = Max[ch_ref];
        while(ylist_ave[n] > (Amax[ch_ref]+ped[ch_ref])*0.5)  Max[ch_ref] = n--;

        t_ref.push_back(Max[ch_ref]);
        // t_ref.push_back(Max[ch_ref]);
        id_ref.push_back(TString(string(txtFile).substr(string(txtFile).length()-9,5)).Atoi());

        for( n = Max[ch_ref] - length/5; n < Max[ch_ref] + length/5; n++)
        {
            // t[0][n-Max[0]+length/10]=xlist_ave[n]-Tmax[0];
            t[ch_ref][n-Max[ch_ref]+length/5]=xlist_ave[n]-xlist_ave[Max[ch_ref]];
            sum[ch_ref][n-Max[ch_ref]+length/5]=sum[ch_ref][n-Max[ch_ref]+length/5]+ylist_ave[n]-ped[ch_ref];
        }
        // if(i++ > Ns)  break;
    }
    cout<<i<<" signals of the reference"<<endl;
    
    //DUT processing
    i = 0;
    for(j = 0; j < id_ref.size(); j++)
    {
        //if(j > 100)  break;
        //cout<<j<<endl;
        if(list_1.size() < 1)   break;
        txtFile = list_1[id_ref[j]];
        vector<double>  xlist_ave, ylist_ave;

        txtFileFullPath = dataPath + txtFile;
        //cout<<txtFileFullPath<<endl;
        TTree *tree = new TTree("waveforms", "waveforms of oscillscope");
        tree->ReadFile(txtFileFullPath, "Time/D:Ampl/D", ',');
        //tree->ReadFile(txtFileFullPath, "Time/D:Ampl/D", ' ');
        
        // Start To read file into TTree
        tree->SetBranchAddress("Time", &Time);
        tree->SetBranchAddress("Ampl", &Ampl);

        for (n = 0; n < tree->GetEntries(); n++)
        {
            tree->GetEntry(n);
            xlist_ave.push_back(Time*1e9);
            // ylist_ave.push_back(Ampl);
            ylist_ave.push_back(-Ampl);
        }
        
        tree->Delete("");
        
        if(length > xlist_ave.size())   length = xlist_ave.size();
        //cout<<length<<endl;

        PedestalAndNoise(ylist_ave,&ped[ch_dut],&noise[ch_dut],&noise_amp[ch_dut]);
        
        Peak(xlist_ave,ylist_ave,&Tmax[ch_dut],&Amax[ch_dut],&Max[ch_dut],threshold,init_f,fin_f);

        if(Amax[ch_dut]*1000 < 50)  continue;

        // n = Max[ch_dut];
        // while(ylist_ave[n] > 0.5*Amax[ch_dut])  Max[ch_dut] = n--;

        for( n = t_ref[j] - length/5; n < t_ref[j] + length/5; n++)
        // for( n = Max[ch_dut] - length/10; n < Max[ch_dut] + length/10; n++)
        {
            // t[0][n-Max[0]+length/10]=xlist_ave[n]-Tmax[0];
            t[ch_dut][n-t_ref[j]+length/5]=xlist_ave[n]-xlist_ave[t_ref[j]];
            sum[ch_dut][n-t_ref[j]+length/5]=sum[ch_dut][n-t_ref[j]+length/5]+ylist_ave[n]-ped[ch_dut];
            
            // t[ch_dut][n-Max[ch_dut]+length/10]=xlist_ave[n]-xlist_ave[Max[ch_dut]];
            // sum[ch_dut][n-Max[ch_dut]+length/10]=sum[ch_dut][n-Max[ch_dut]+length/10]+ylist_ave[n]-ped[ch_dut];
        }
        i++;
    }
    cout<<i<<" signals of the DUT"<<endl;

    TPaveText *pt;
    TText *text;
    vector<vector<double>> Tlist_ave(4), Alist_ave(4);

    ofstream myfile;
    myfile.open (plotStorePath + plotType + "average.txt");

    TGraph *gr_ave[4];

    for(n = 0; n < 4; n++)
    {
        gr_ave[n] = new TGraph();
        for(i = 0; i < length*2/5; i++)
        {
            Tlist_ave[n].push_back(t[n][i]);
            Alist_ave[n].push_back(sum[n][i]/Ns);
            gr_ave[n]->SetPoint(i,t[n][i],sum[n][i]/Ns);
            if(n == 1)  myfile<<t[n][i]<<","<<sum[n][i]<<endl;     
        }
    }

    int P_left[4] = {0}, P_right[4] = {0};
    fin_f = 0.4;

    for(i = 0; i < 4; i++)
    {
        //Integral range
        Peak(Tlist_ave[i],Alist_ave[i],&Tmax[i],&Amax[i],&Max[i],threshold,0.5,0.9);
        
        n = Max[i];
        while(Alist_ave[i][n] > Amax[i]*0.1)  P_left[i] = n--;
        
        n = Max[i];
        while(Alist_ave[i][n] > Amax[i]*0.1)  P_right[i] = n++;

        cout<<"Length:"<<Alist_ave[i].size()<<",M:"<<Max[i]<<",L:"<<P_left[i]<<",R:"<<P_right[i]<<endl;
        myfile<<"Length:"<<Alist_ave[i].size()<<",M:"<<Max[i]<<",L:"<<P_left[i]<<",R:"<<P_right[i]<<endl;

        //Draw average waveform
        double y_max;
        y_max = TMath::MaxElement(gr_ave[i]->GetN(),gr_ave[i]->GetY());

        gr_ave[i]->SetMarkerStyle(21);
        gr_ave[i]->SetMarkerSize(0.4);
        gr_ave[i]->SetMarkerColor(2);
        gr_ave[i]->SetLineColor(2);

        gr_ave[i]->GetXaxis()->SetTitle("Time [ns]");
        gr_ave[i]->GetYaxis()->SetTitle("Ampl [V]");

        gr_ave[i]->GetXaxis()->SetNdivisions(505, kTRUE);
        gr_ave[i]->GetXaxis()->SetLimits(-5, 10);
        // gr_ave[i]->GetYaxis()->SetRangeUser(-0.008, 0.1);
        gr_ave[i]->GetYaxis()->SetRangeUser(-0.08*y_max, 1.1*y_max);

        gr_ave[i]->Draw("ALP");

        //slope of average waveform
        rise_start = leftTime(Tlist_ave[i],Alist_ave[i],Amax[i]*0.4,Max[i]);
        rise_end = leftTime(Tlist_ave[i],Alist_ave[i],Amax[i]*0.6,Max[i]);
        risetime_ave[i] = rise_end - rise_start;

        //slope of average waveform
        // for(j = 1; j < 5; j++)
        // {
        //     rise_start = leftTime(Tlist_ave[i],Alist_ave[i],Amax[i]*0.1*j,Max[i]);
        //     rise_end = leftTime(Tlist_ave[i],Alist_ave[i],Amax[i]*(1-0.1*j),Max[i]);
        //     risetime_ave[i] = rise_end - rise_start;
            
        //     fa = (1 - 0.2*j)*Amax[i]/risetime_ave[i];
        //     fb = Amax[i]*0.1*j - fa*rise_start;
        //     cout<<Amax[i]*0.1*j<<","<<rise_start<<endl;
        //     cout<<fa<<","<<fb<<endl;

        //     f = new TF1("f","[0]+[1]*x",1,2.5);
        //     f->SetLineColor(colorList[j]);
        //     f->SetLineWidth(1);
        //     f->SetParameter(0,fb);
        //     f->SetParameter(1,fa);
        //     f->Draw("same");
        // }

        slope[i] = Amax[i]*0.2/(rise_end-rise_start);
        cout<<"Amax_20:"<<Amax[i]*0.2<<",Rise time:"<<rise_end-rise_start<<",Average slope:"<<slope[i]<<endl;
        myfile<<"Amax_20:"<<Amax[i]*0.2<<",Rise time:"<<rise_end-rise_start<<",Average slope:"<<slope[i]<<endl;

        pt=new TPaveText(1.2,0.75,0.6,0.6,"NDC");
        pt->SetBorderSize(0);
        pt->SetFillColor(0);
        pt->SetFillStyle(0);
        pt->SetTextAlign(12);
        pt->SetTextSize(0.02);

        text=pt->AddText(Form("Amax_20: %.4f",Amax[i]*0.2*1000));
        text=pt->AddText(Form("Rise time_4060: %.4f",rise_end-rise_start));
        text=pt->AddText(Form("Slope: %.4f",slope[i]*1000));

        pt->Draw();

        c1->Update();
        c1->SaveAs(plotStorePath + plotType + i + "average.png");
        c1->SaveAs(plotStorePath + plotType + i + "average.svg");
        c1->Clear();
    }

    if(average_only)    return;
    // //refresh the directory
    // gSystem->FreeDirectory(plotDir);
    // gSystem->OpenDirectory(dataPath);

    //Draw waveforms & Variables extraction
    TGraph *gr;
    TFile *f=new TFile(plotStorePath + plotType + plotname + ".root","recreate");
    TTree *tt[4];
    tt[0]=new TTree("t1","channel 1");
    tt[1]=new TTree("t2","channel 2");
    tt[2]=new TTree("t3","channel 3");
    tt[3]=new TTree("t4","channel 4");

    for(n = 0; n < 4; n++)
    {
        tt[n]->Branch("Amax",&Amax[n]);
        tt[n]->Branch("Tmax",&Tmax[n]);
        tt[n]->Branch("Charge",&charge[n]);
        tt[n]->Branch("Charge_HW",&charge_hw[n]);
        tt[n]->Branch("Charge_DW",&charge_dw[n]);
        tt[n]->Branch("CTD",&CTD[n]);
        tt[n]->Branch("CFD",&CFD[n]);
        tt[n]->Branch("ZCD",&ZCD[n]); 
        tt[n]->Branch("TOT",&TOT[n]);
        tt[n]->Branch("risetime_ave",&ave_risetime[n]);
        tt[n]->Branch("risetime",&risetime[n]);
        tt[n]->Branch("slope",&slope[n]);
        tt[n]->Branch("pedestal",&ped[n]);
        tt[n]->Branch("ped_1",&ped_1[n]);
        tt[n]->Branch("ped_2",&ped_2[n]);
        tt[n]->Branch("Noise_Amp",&noise_amp[n]);
        tt[n]->Branch("Noise",&noise[n]);
        tt[n]->Branch("ID",&id[n]);
    }

    i = 0;
    j = 0;
    TH2D *hA_T[4];
    int t_range[2] = {-20,-10};

    for(n = 0; n < 4; n++)
    {
        // if(n == 1 || n == 2)
        // {
        //     t_range[0] = -5;
        //     t_range[1] = 5;
        // }
        hA_T[n] = new TH2D("Amax vs. Tmax","Amax vs. Tmax",100,0,200,40,t_range[0],t_range[1]);
    }

    for(j = 0; j < dirList.size(); j++)
    {
        // if(j > 50)  break;
        //cout<<j<<endl;
        txtFile = dirList[j];
        dual_board = 0;
        
        // myfile<<txtFile<<endl;
        TGraph *g_rising;
        double RefT = -99999;
        double MaxTime = -99999;
        vector<double> xlist,ylist,ylist_1,ylist_2,Tlist,Alist,ZCD_Alist;
        
        if(txtFile.Contains("C1"))  ch_num = 0;
        else if(txtFile.Contains("C2"))  ch_num = 1;
        else if(txtFile.Contains("C3"))  ch_num = 2;
        else if(txtFile.Contains("C4"))  ch_num = 3;
        else ch_num = 0;
        
        // if(ch_num == 0) continue;
        txtFileFullPath = dataPath + txtFile;
        //cout<<txtFileFullPath<<endl;
        TString plotInfo = txtFile;
        plotInfo = plotInfo.ReplaceAll(".txt", "");
        TTree *tree = new TTree("waveforms", "waveforms of oscillscope");
        //tree->ReadFile(txtFileFullPath, "Time/D:Ampl/D", ' ');
        tree->ReadFile(txtFileFullPath, "Time/D:Ampl/D", ',');
        
        // Start To read file into TTree
        tree->SetBranchAddress("Time", &Time);
        tree->SetBranchAddress("Ampl", &Ampl);
        gr = new TGraph();

        for (n = 0; n < tree->GetEntries(); n++)
        {
            tree->GetEntry(n);
            //cout << Time << "  " << Ampl << " " << endl;
            xlist.push_back(Time*1e9);
            // ylist.push_back(Ampl);
            ylist.push_back(-Ampl);
        }
        tree->Delete("");

        if(length > xlist.size())   length = xlist.size();

        //DUT
        PedestalAndNoise(ylist,&ped[ch_num],&noise[ch_num],&noise_amp[ch_num]);
        for( n = 0; n < length; n++)
        {
            Alist.push_back(ylist[n]-ped[ch_num]);
            //Tlist.push_back(xlist[n]-TOA_Ref);

            // if(n < length/50) ylist_1.push_back(ylist[n]);
            // if(n > length/50 && n < length/25) ylist_2.push_back(ylist[n]);

            // sum[ch_num][n]=sum[ch_num][n]+Alist[n];
            t[ch_num][n]=xlist[n];
            if(draw)    gr->SetPoint(n, xlist[n], Alist[n]);
        }
        
        // PedestalAndNoise(ylist_1,&ped_1[ch_num],&noise_1[ch_num],&noise_amp_1[ch_num]);
        // PedestalAndNoise(ylist_2,&ped_2[ch_num],&noise_2[ch_num],&noise_amp_2[ch_num]);

        // Get the peak information
        //delta_amax = Peak(xlist,Alist,&Tmax[ch_num],&Amax[ch_num],&Max[ch_num],threshold);
        init_f = 0.1;
        fin_f = 0.9;
        if(ch_num == ch_dut || ch_num == ch_ref)
        {
            init_f = 0.2;
            fin_f = 0.44;
            // init_f = 0.24;
            // fin_f = 0.36;
        }
        // if(ch_num == ch_dut) init_f = 0.2;
        Peak(xlist,Alist,&Tmax[ch_num],&Amax[ch_num],&Max[ch_num],threshold,init_f,fin_f);
        //cout<<"Amax="<<Amax<<endl;
        //cout<<"Max="<<Max<<endl;
        //if(Amax < threshold) continue;

        charge[ch_num] = Charge(xlist,Alist,Max[ch_num],P_left[ch_num],P_right[ch_num],threshold)/transImpedance;
        charge_hw[ch_num] = Charge(xlist,Alist,Max[ch_num],0.8*P_left[ch_num],0.8*P_right[ch_num],threshold)/transImpedance;
        charge_dw[ch_num] = Charge(xlist,Alist,Max[ch_num],2*P_left[ch_num],2*P_right[ch_num],threshold)/transImpedance;
        // if(ch_num == 1) cout<<charge_dw[ch_num]*1e12<<endl;

        CTD[ch_num] = leftTime(xlist,Alist,th,Max[ch_num]);
        //CTD_C = rightTime(xlist,Alist,th,Max);
        //TOT = CTD_C - CTD; 
        CFD[ch_num] = leftTime(xlist,Alist,Amax[ch_num]*fra,Max[ch_num]);
        //cout<<"TOT="<<TOT<<endl

        // simple_rise_start = -99999;
        // simple_rise_end = -99999;
        
        // int start_P = 0, end_P = 0;
        // g_rising = new TGraph;

        // for(n = Max[ch_num]; n > 0; n--)
        // {
        //     if(simple_rise_end < Tlist[n] && Alist[n] < Amax[ch_num]*0.6)
        //     {
        //         simple_rise_end = Tlist[n];
        //         end_P = n;
        //     }

        //     //if(end_P > 0)   g_rising->SetPoint(end_P - n, Tlist[n], Alist[n]*1000);

        //     if(simple_rise_end > Tlist[n] && simple_rise_start < Tlist[n] && Alist[n] < Amax[ch_num]*0.4)
        //     {
        //         simple_rise_start = Tlist[n];
        //         start_P = n;
        //         break;
        //     }
        // }
        // //cout<<"90:"<<simple_rise_end<<",10:"<<simple_rise_start<<endl;
        // simple_risetime[ch_num] = simple_rise_end - simple_rise_start;

        //g_rising->Fit("pol1","Q","", simple_rise_start + 1, simple_rise_end);
        //linearity = g_rising->GetFunction("pol1")->GetChisquare();
        
        //g_rising->GetXaxis()->SetNdivisions(505, kTRUE);
        //g_rising->GetXaxis()->SetLimits(-2 , 2);
        //g_rising->GetYaxis()->SetRangeUser(-8, 25);

        //g_rising->Draw("AP");
        //c1->SaveAs("linearity.png");
        //cout<<linearity<<endl;

        //linearity = (Alist[end_P] - Alist[(end_P+start_P)/2])/(Alist[(end_P+start_P)/2] - Alist[start_P]);

        rise_start = leftTime(xlist,Alist,Amax[ch_num]*0.4,Max[ch_num]);
        rise_end = leftTime(xlist,Alist,Amax[ch_num]*0.6,Max[ch_num]);
        risetime[ch_num] = rise_end - rise_start;
        slope[ch_num] = Amax[ch_num]*0.2/risetime[ch_num];
        TOT[ch_num] = rightTime(xlist,Alist,Amax[ch_num]*0.1,Max[ch_num]) - leftTime(xlist,Alist,Amax[ch_num]*0.1,Max[ch_num]);

        ave_risetime[ch_num] = risetime_ave[ch_num];
        
        id[ch_num] = TString(string(txtFile).substr(string(txtFile).length()-9,5)).Atoi();
        // if(doZCD == 1)
        // {
        //     for(n = 0; n < length-shift-1; n++)  ZCD_Alist.push_back(Alist[n]-ZCD_fra*Alist[n+shift]);
        //     //cout<<ZCD_Alist[n-shift-1]<<endl;   
        //     //gr->SetPoint(n, xlist[n], ZCD_Alist[n-shift-1]);
        //     ZCD = leftTime(xlist,ZCD_Alist,0,Max);
        // }

        // TOF_CTD = CTD1 - CTD;
        // TOF_CFD = CFD1 - CFD;
        // TOF_ZCD = ZCD1 - ZCD;

        //cout<<TOF<<endl;
        tt[ch_num]->Fill();

        //T vs. A
        hA_T[ch_num]->Fill(Amax[ch_num]*1000,Tmax[ch_num]);
        //cout<<Amax*1000<<endl;

        // if(draw && Amax[ch_num]*1000 > 30 && ch_num == 1 && i < 50)
        if(draw && ch_num == 1)
        {            
            i = i + 1;
            //Start to draw
            gr->SetMarkerStyle(21);
            gr->SetMarkerSize(0.4);
            gr->SetMarkerColor(kRed);
            gr->SetLineColor(kRed);

            gr->GetXaxis()->SetTitle("Time [ns]");
            gr->GetYaxis()->SetTitle("Ampl [V]");

            //gr->GetXaxis()->SetLimits(xlist[0], xlist[length-1]);
            gr->GetXaxis()->SetLimits(-20, 20);
            gr->GetXaxis()->SetNdivisions(505, kTRUE);
            gr->GetYaxis()->SetRangeUser(-0.002, 0.20);

            // if (i%N == 0)
                gr->Draw("ALP");
            // else
            //     gr->Draw("LP");

            // Choose plot type
            if(multi)
            {
                i = i + 1;
                if(i%N == 0)
                {
                    c1->Update();
                    c1->SaveAs(plotStorePath + plotType + plotname + "_" + i/N +".png");
                    //c1->SaveAs(plotStorePath + plotType + plotname + "_" + i/N +".svg");
                    c1->Clear();
                }
            }
            else
            {
                c1->Update();
                c1->SaveAs(plotStorePath + plotType + plotInfo + ".png");
                //c1->SaveAs(plotStorePath + plotType + plotInfo + ".svg");
                c1->Clear();
            }
        }
        printf("Processing:%.2f%%\r",j*100.0/dirList.size());
    }
    cout<<"files number:"<<j<<endl;
    //cout<<"multi="<<multi<<endl;

    for(n = 0; n < 4; n++)  tt[n]->Write();

    // Draw Amax vs. Tmax distribution
    for(n = 0; n < 4; n++)
    {
        hA_T[n]->GetXaxis()->SetTitle("Amax [mV]");
        hA_T[n]->GetYaxis()->SetTitle("Tmax [ns]");
        hA_T[n]->GetXaxis()->SetNdivisions(505, kTRUE);
        hA_T[n]->SetMarkerSize(0.2);
        hA_T[n]->Draw("COLZ");

        c1->SaveAs(plotStorePath + plotType + n + "Amax_vs_Tmax.png");
        c1->Clear();
    }
}

void PedestalAndNoise(vector<double> ylist, double *ped, double *noise, double *noise_amp)
{
    int i, length = ylist.size();
    double mean = 0, stdev = 0, nmax = 0;
 
    for(i = 0; i < length*0.2 ; i++)
    {
        mean+=ylist[i];
        if(ylist[i]>nmax)   nmax = ylist[i];
    }
    mean = mean/(length*0.2);

    for(i = 0; i < length*0.2 ; i++) stdev+=(ylist[i]-mean)*(ylist[i]-mean);
    stdev = sqrt(stdev/(length*0.2-1));

    *ped = mean;
    *noise = stdev;
    *noise_amp = nmax - mean;
}

double Peak(vector<double> xlist, vector<double> Alist, double *Tmax, double *Amax, int *Max, double threshold, double init_f, double fin_f)
{
    int k = 0, length = xlist.size();
    double amax = threshold, amax_sub = threshold, delta = 0;
    //cout<<"m="<<m<<endl;

    *Tmax = 0;
    *Amax = 0;
    *Max = 0;

    if (length != 0)
    {
        for(k = length*init_f; k < length*fin_f ; k++)
        {
            // if(Alist[k] > Alist [k - 1] && Alist[k] > Alist[k + 1] && Alist[k] > amax)
            if(Alist[k] > amax)
            {   
                amax_sub = amax;
                amax = Alist[k];
                *Tmax = xlist[k];
                *Amax = Alist[k];
                *Max = k;
            }
        } 
    }
    delta = amax - amax_sub;
    //cout<<amax_sub<<endl;
    
    return delta*1000;
}

double leftTime(vector<double> xlist, vector<double> Alist, double threshold, int Max)
{
    int i = Max;
    double a = 0,b = 0,lefttime = 0;
    
    while(Alist[i]>threshold) i--;

    a=(Alist[i+1]-Alist[i])/(xlist[i+1]-xlist[i]);
    b=Alist[i]-a*xlist[i];
    lefttime=(threshold-b)/a;

    return lefttime;
}

double rightTime(vector<double> xlist, vector<double> Alist, double threshold, int Max)
{
    int i = Max;
    double a = 0,b = 0,righttime = 0;

    while(Alist[i]>threshold) i++;

    a=(Alist[i]-Alist[i-1])/(xlist[i]-xlist[i-1]);
    b=Alist[i]-a*xlist[i];
    righttime=(threshold-b)/a;

    return righttime;
}

double Charge(vector<double> xlist, vector<double> Alist, int Max, int left, int right, double threshold)
{
    int i,width;
    double deltaT,charge;
    width = right - left;
    deltaT = xlist[1] - xlist[0];
    //cout<<deltaT<<endl;

    for(i = Max - 0.5*width; i < Max + width; i++)//1e-9/deltaT = points in 1ns
    // for(i = Max - 10; i < Max + width*2/3; i++)//1e-9/deltaT = points in 1ns
    {
        //if(Alist[i]>threshold) 
            charge+=Alist[i];
    }
    
    charge=charge*deltaT;
    return charge;
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
        oscilloscope_draw(dataPath,draw,multi);
    }
    return 0;
}
