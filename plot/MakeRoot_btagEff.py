import math
import gc
import sys
import ROOT
import numpy as np
import copy
import os
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1;")
ROOT.TH1.AddDirectory(ROOT.kFALSE)
ROOT.gStyle.SetOptStat(0)
from array import array
from ROOT import TColor
from ROOT import TGaxis
from ROOT import THStack
import gc
TGaxis.SetMaxDigits(2)

num=['h2_BTaggingEff_Num_b', 'h2_BTaggingEff_Num_c', 'h2_BTaggingEff_Num_udsg']
denom=['h2_BTaggingEff_Denom_b', 'h2_BTaggingEff_Denom_c', 'h2_BTaggingEff_Denom_udsg']
Name = ['h2_BTaggingEff_b', 'h2_BTaggingEff_c', 'h2_BTaggingEff_udsg']
HistAddress = '/user/rgoldouz/NewAnalysis2020/Analysis/hists/'

Samples = ['mc_2016.root','mc_2017.root','mc_2018.root']
year = ['2016', '2017', '2018']
Hists = []

for numyear, nameyear in enumerate(Samples):
    f = ROOT.TFile.Open(HistAddress + nameyear)
    for numvar, namevar in enumerate(denom):
        h1= f.Get(num[numvar])
        h2= f.Get(denom[numvar])
        h1.Divide(h2)
        h1.SetName(year[numyear] + '_' + Name[numvar])
        Hists.append(h1)
    f.Close()

myfile = ROOT.TFile( 'btagEff.root', 'RECREATE')
for H in Hists:
    H.Write()
