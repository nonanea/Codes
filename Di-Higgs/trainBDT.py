#!/usr/bin/env python

# 
# This example is basically the same as $ROOTSYS/tmva/test/TMVAClassification.C
# 

import ROOT
import sys
import os

# in order to start TMVA
from ROOT import TMVA
ROOT.TMVA.Tools.Instance()

# note that it seems to be mandatory to have an
# output file, just passing None to TMVA::Factory(..)
# does not work. Make sure you don't overwrite an
# existing file.

# treename = sys.argv[1:]
# print(treename)
# sig_path = "Ntuples/SR/signal/"
# bkg_path = "Ntuples/SR/background/"

# sig_list = os.listdir(sig_path)
# bkg_list = os.listdir(bkg_path)

sig = "Ntuples/SR/signal/nominal_Preselection.root"
bkg = "Ntuples/SR/background/nominal_Preselection.root"

# open input file, get trees, create output file
# for i in range(len(sig_list)):

# treename = sig_list[i].replace("_Preselection.root","")
treename = "nominal"

# f_s = ROOT.TFile(sig_path + sig_list[i])
f_s = ROOT.TFile(sig)
tree_s = f_s.Get(treename)
# f_b = ROOT.TFile(bkg_path + bkg_list[i])
f_b = ROOT.TFile(bkg)
tree_b = f_b.Get(treename)

n_sig = tree_s.GetEntries()*0.8
n_bkg = tree_b.GetEntries()*0.8

fout = ROOT.TFile("Training/"+treename+"training.root","RECREATE")

# define factory with options
factory = ROOT.TMVA.Factory(
                           # "TMVAClassification",
                           treename,
                           fout,
                           ":".join([    "!V",
                                          "!Silent",
                                          "Color",
                                          "DrawProgressBar",
                                          "Transformations=I;D;P;G,D",
                                          "AnalysisType=Classification"]
                                    ))


dataloader = TMVA.DataLoader("dataset")

# add discriminating variables for training
dataloader.AddVariable("lep_Pt_0","F")
dataloader.AddVariable("lep_Pt_1","F")
dataloader.AddVariable("lep_Pt_2","F")
dataloader.AddVariable("lep_Pt_3","F")
dataloader.AddVariable("lep_Etcone30_0","F")
dataloader.AddVariable("lep_Etcone30_1","F")
dataloader.AddVariable("lep_Etcone30_2","F")
dataloader.AddVariable("lep_Etcone30_3","F")
# dataloader.AddVariable("SumLabel := lep_iso_0+lep_iso_1+lep_iso_2+lep_iso_3","F")
# dataloader.AddVariable("lep_Eta_0","F")
dataloader.AddVariable("lep_Eta_1","F")
dataloader.AddVariable("lep_Eta_2","F")
dataloader.AddVariable("lep_Eta_3","F")
# dataloader.AddVariable("jet_E_0","F")
# dataloader.AddVariable("jet_E_1","F")
dataloader.AddVariable("jet_Pt_0","F")
# dataloader.AddVariable("jet_Pt_1","F")
# dataloader.AddVariable("jet_Phi_0","F")
# dataloader.AddVariable("jet_Phi_1","F")
dataloader.AddVariable("m_12","F")
dataloader.AddVariable("m_34","F")
dataloader.AddVariable("m_4l","F")
dataloader.AddVariable("m_jj","F")
# dataloader.AddVariable("p_4l","F")
dataloader.AddVariable("p_jj","F")
# dataloader.AddVariable("cs_jet","F")
dataloader.AddVariable("HT","F")
# dataloader.AddVariable("HT_lep","F")
# dataloader.AddVariable("cs_lep_12","F")
# dataloader.AddVariable("cs_lep_34","F")
# dataloader.AddVariable("cs_pairs","F")
dataloader.AddVariable("met_met","F")
dataloader.AddVariable("Dphi_met_jets","F")
# dataloader.AddVariable("cs_Z_pair","F")
# dataloader.AddVariable("njets","F")
dataloader.AddVariable("nbjets","F")

# define signal and background trees
dataloader.AddSignalTree(tree_s)
dataloader.AddBackgroundTree(tree_b)

dataloader.SetSignalWeightExpression( "mcWeight" )
dataloader.SetBackgroundWeightExpression( "mcWeight" )

# define additional cuts 
sigCut = ROOT.TCut("1")
bgCut = ROOT.TCut("1")

train_sig = "nTrain_Signal=" + str(n_sig)
train_bkg = "nTrain_Background=" + str(n_bkg)
test_sig = "nTest_Signal=" + str(n_sig)
test_bkg = "nTest_Background=" + str(n_bkg)

# set options for trainings
dataloader.PrepareTrainingAndTestTree(sigCut, 
                                 bgCut, 
                                 ":".join([train_sig,
                                             train_bkg,
                                             test_sig,
                                             test_bkg,
                                             "SplitMode=Random",
                                             "NormMode=NumEvents",
                                             "!V"
                                             ]))

# book and define methods that should be trained

factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDT",
                           ":".join([ "!H",
                                       "!V",
                                       "NTrees=1000",
                                       "MinNodeSize=3%",
                                       "MaxDepth=2",
                                    #    "BoostType=AdaBoost",
                                    #    "AdaBoostBeta=0.5",
                                       "BoostType=Grad",
                                       "Shrinkage=0.3",
                                       "SeparationType=GiniIndex",
                                       "nCuts=10",
                                       "PruneMethod=NoPruning",
                                       ]))

# self-explaining
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

# print(sig_list[i])