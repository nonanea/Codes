import os
import sys

sig = ['005','009','015','019','023','027','031','035','039','045','051','054','060','066','069','075']
# sig = ['005','009','015','019','023','027','031','035']
# sig = ['039','045','051','054','060','066','069','075']

res = [
        0.082970911
        ,0.155532186
        ,0.254278777
        ,0.322111209
        ,0.391810168
        ,0.460365453
        ,0.518714803
        ,0.591198869
        ,0.654172172
        ,0.770273591
        ,0.893132374
        ,0.951376734
        ,1.106011789
        ,1.320550183
        ,1.438338998
        ,1.794264318
]

sh = open("dofit.sh","a")
sh.seek(0)
sh.truncate()

poi = "XSec"
# poi = "rootgZp"

comment = "%"

if poi == "rootgZp":    comment = ""
    # POI = poi
    # LHscan = poi

for i,mass in enumerate(sig):
    # if mass != '031':   continue

    # config = "config/Zp"+mass+"_fit.config"
    
    cg = open("config/Zp"+mass+"_fit.config","a")
    cg.seek(0)
    cg.truncate()

    job = "Job:\"Zp"+mass+"_fit\"\n Label: L_mu-L_tau Zp\n CmeLabel: 13 TeV\n LumiLabel: 139 fb^{-1}\n %Lumi: 139\n POI: " + poi + "\n ReadFrom: NTUP\n NtuplePath: \"/Lab/Z_prime/lllv/Minitrees/Outputs/MVA/no_mll/high_39/\"\n NtupleName: nominal\n MCweight: weight\n %Selection: mll_1 < 80 && lep0_isoTight == 1 && lep1_isoTight == 0 && lep1_pt > 10 && lep2_pt > 8 && n_jets == 0\n Selection: DNN_" + str(int(mass)) + " > 0.8\n" + " %Selection: mll_1 > " + str(int(mass)-5*res[i]) + " && mll_1 < " + str(int(mass)+5*res[i]) + "\n" + " DebugLevel: 2\n UseGammaPulls: TRUE\n PlotOptions: NORMSIG,NOSIG\n ImageFormat: png,svg\n\n"

    fit = "Fit: \"myFit\"\n FitType: SPLUSB\n FitRegion: CRSR\n POIAsimov: 0\n %doLHscan: "+ poi + "\n %LHscanSteps: 50\n %UseMinos: " + poi + "\n\n"

    SR = "Region: \"SR\"\n Type: SIGNAL\n Variable: \"mll_1\",20," + str(int(mass)-5*res[i]) + "," + str(int(mass)+5*res[i]) + "\n %Variable: \"DNN_" + str(int(mass)) + "\",50,0.,1.0\n VariableTitle: \"m_{ll,2} [GeV]\"\n Label: \"SR\"\n ShortLabel: \"SR\"\n LogScale: TRUE\n DataType: ASIMOV\n\n"

    # SR = "Region: \"SR\"\n Type: SIGNAL\n %Variable: \"mll_1\",80,0,80\n Variable: \"DNN_" + str(int(mass)) + "\",50,0.,1.0\n VariableTitle: \"m_{ll,2} [GeV]\"\n Label: \"SR\"\n ShortLabel: \"SR\"\n LogScale: TRUE\n DataType: ASIMOV\n\n"

    SR_sample1 = "Sample: \"Background\"\n Type: BACKGROUND\n Title: \"Background\"\n FillColor: 4\n LineColor: 4\n NtupleFile: \"bkg\"\n Regions: SR\n\n"

    signal = "Sample: \"Signal\"\n Type: SIGNAL\n Title: \"Zp\"\n FillColor: 2\n LineColor: 2\n Selection: mass == " + str(int(mass)) + "\n NtupleFile: \"sig\"\n " + comment + "MCweight: \"weight_g\"\n Regions: SR\n\n"

    # NF1 = "NormFactor:\"Xsec\"\n Title: \"Xsec\"\n Nominal:1\n Min: -10\n Max: 10\n Samples: Signal\n Regions: SR\n\n"

    NF1 = "NormFactor:\""+ poi.replace("root","") + "\"\n Title: \""+ poi.replace("root","") + "\"\n Nominal:1\n Min: -2\n Max: 10\n " + comment + "Expression: (rootgZp*rootgZp):rootgZp[0,2.]\n Samples: Signal\n Regions: SR\n\n"

    syst1 = "Systematic: \"luminosity\"\n Title: \"luminosity\"\n Type: OVERALL\n OverallUp: 0.017\n OverallDown: -0.017\n Samples: all\n Category: Instrumental\n\n"

    syst2 = "Systematic: \"bkgxsec\"\n Title: \"bkgxsec\"\n Type: OVERALL\n OverallUp: 0.1\n OverallDown: -0.1\n Samples: all\n Category: Instrumental\n\n"

    limit = "Limit: \"limit\"\n LimitType: ASYMPTOTIC\n"

    cg.write(job+fit+SR+SR_sample1+signal+NF1+syst1+syst2)
    cg.close()

    sh.write("trex-fitter -ndwfl config/Zp"+mass+"_fit.config\n")

sh.close()