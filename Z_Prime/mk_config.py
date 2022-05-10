import os
import sys

# sig = ['005','009','015','019','023','027','031','035','039','045','051','054','060','066','069','075']
sig = ['005','009','015','019','023','027','031','035']
# sig = ['039','045','051','054','060','066','069','075']

sh = open("dofit.sh","a")
sh.seek(0)
sh.truncate()

for mass in sig:
    # if mass != '031':   continue

    # config = "config/Zp"+mass+"_fit.config"
    
    cg = open("config/Zp"+mass+"_fit.config","a")
    cg.seek(0)
    cg.truncate()

    job = "Job:\"Zp"+mass+"_fit\"\n Label: L_mu-L_tau Zp\n CmeLabel: 13 TeV\n LumiLabel: 139 fb^{-1}\n Lumi: 139\n POI: XSec\n %POI: rootgZp\n ReadFrom: NTUP\n NtuplePath: \"/Lab/Z_prime/lllv/Minitrees/Outputs/MVA/low_low/\"\n NtupleName: nominal\n MCweight: weight\n %Selection: mll_1 < 80 && lep0_isoTight == 1 && lep1_isoTight == 0 && lep1_pt > 10 && lep2_pt > 8 && n_jets == 0\n %Selection: DNN_" + str(int(mass)) + " > 0.8\n DebugLevel: 2\n UseGammaPulls: TRUE\n PlotOptions: NORMSIG,NOSIG\n ImageFormat: png,svg\n\n"

    fit = "Fit: \"myFit\"\n FitType: BONLY\n FitRegion: CRSR\n POIAsimov: 1\n LHscanSteps: 500\n UseMinos: XSec\n %UseMinos: gZp\n\n"

    SR = "Region: \"SR\"\n Type: SIGNAL\n %Variable: \"mll_1\",20,0,80\n Variable: \"DNN_" + str(int(mass)) + "\",20,.0,1.0\n VariableTitle: \"m_{ll,1} [GeV]\"\n Label: \"SR\"\n ShortLabel: \"SR\"\n LogScale: TRUE\n DataType: ASIMOV\n\n"

    SR_sample1 = "Sample: \"Background\"\n Type: BACKGROUND\n Title: \"Background\"\n FillColor: 4\n LineColor: 4\n NtupleFile: \"bkg\"\n Regions: SR\n\n"

    signal = "Sample: \"Signal\"\n Type: SIGNAL\n Title: \"Zp\"\n FillColor: 2\n LineColor: 2\n Selection: mass == " + str(int(mass)) + "\n NtupleFile: \"sig\"\n %MCweight: \"weight_g\"\n Regions: SR\n\n"

    # NF1 = "NormFactor:\"Xsec\"\n Title: \"Xsec\"\n Nominal:1\n Min: -10\n Max: 10\n Samples: Signal\n Regions: SR\n\n"

    NF1 = "NormFactor:\"XSec\"\n Title: \"XSec\"\n%NormFactor:\"gZp\"\n %Title: \"gZp\"\n Nominal:1\n Min: -10\n Max: 10\n %Expression: (rootgZp*rootgZp):rootgZp[0,2.]\n Samples: Signal\n Regions: SR\n\n"

    syst1 = "Systematic: \"luminosity\"\n Title: \"luminosity\"\n Type: OVERALL\n OverallUp: 0.017\n OverallDown: -0.017\n Samples: all\n Category: Instrumental\n\n"

    syst2 = "Systematic: \"bkgxsec\"\n Title: \"bkgxsec\"\n Type: OVERALL\n OverallUp: 0.1\n OverallDown: -0.1\n Samples: all\n Category: Instrumental\n\n"

    limit = "Limit: \"limit\"\n LimitType: ASYMPTOTIC\n"

    cg.write(job+fit+SR+SR_sample1+signal+NF1+syst1+syst2)
    cg.close()

    sh.write("trex-fitter -nwdfl config/Zp"+mass+"_fit.config\n")

sh.close()