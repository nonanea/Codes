import ROOT

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gROOT.SetBatch(ROOT.kTRUE)

f = ROOT.TFile.Open("user.chihao.25776474._000001.output.root")
l = f.GetListOfKeys()

print(len(l))

l = l[3:len(l)-2]
# del l[0:3]
# del l[len(l)-1,len(l)]
# l = f.GetTrees()

print(len(l))

for name in l:
    print(name.GetName())