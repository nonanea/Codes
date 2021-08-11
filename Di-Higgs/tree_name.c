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

void tree_name(TString file)
{
    int i;
    vector<TString> s;
    // TString ss;

    // TFile *f = new TFile("user.chihao.25776474._000001.output.root");
    TFile *f = new TFile(file);
    TList *l = (TList*) f->GetListOfKeys();

    for(i = 0; i < l->GetEntries(); i++)
    // for(i = 0; i < 5; i++)
    {
        // ss = (*l->At(i)).GetName();
        s.push_back((*l->At(i)).GetName());
    }

    for(i = 0; i < s.size(); i++) cout<<"\""<<s[i]<<"\","<<endl;
    cout<<s.size()<<endl;
}