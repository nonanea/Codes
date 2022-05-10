import ROOT

# sig = ['005','009','015',
sig = ['005','009','015','019','023','027','031','035','039','045','051','054','060','066','069','075']

for mass in sig:
    fit = "Zp"+mass+"_fit"
    lf_name = "/Limits/asymptotics/myLimit_BLIND_CL95.root"

    lf = ROOT.TFile(fit+lf_name,"read")
    lt = lf.Get("stats")

    # print lt
    for event in lt:
        print int(mass),format(event.exp_upperlimit_minus2,'.2g'),format(event.exp_upperlimit_minus1,'.2g'),format(event.exp_upperlimit,'.2g'),format(event.exp_upperlimit_plus1,'.2g'),format(event.exp_upperlimit_plus2,'.2g')