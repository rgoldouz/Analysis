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

def puDraw(hists,Fnames,year,can):
    ratio=[]
    canvas = ROOT.TCanvas(can,can,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    legend = ROOT.TLegend(0.7,0.65,0.9,0.88)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)

    pad1=ROOT.TPad("pad1", "pad1", 0, 0.315, 1, 0.99 , 0)#used for the hist plot
    pad2=ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.305 , 0)#used for the ratio plot
    pad1.Draw()
    pad2.Draw()
    pad2.SetGridy()
    pad2.SetTickx()
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.05)
    pad2.SetTopMargin(0.1)
    pad2.SetBottomMargin(0.4)
    pad2.SetLeftMargin(0.14)
    pad2.SetRightMargin(0.05)
    pad2.SetFillStyle(0)
    pad1.SetFillStyle(0)
    pad1.cd()
    pad1.SetLogx(ROOT.kFALSE)
    pad2.SetLogx(ROOT.kFALSE)
    pad1.SetLogy(ROOT.kFALSE)

    y_min=0
    y_max=1.2*hists[0].GetMaximum()
    hists[0].SetTitle("")
    hists[0].GetYaxis().SetTitle('Probability')
    hists[0].GetXaxis().SetLabelSize(0)
    hists[0].GetYaxis().SetTitleOffset(0.8)
    hists[0].GetYaxis().SetTitleSize(0.07)
    hists[0].GetYaxis().SetLabelSize(0.04)
    hists[0].GetYaxis().SetRangeUser(y_min,y_max)
    hists[0].Draw("hist")
    nomhist = hists[0].Clone()
    for H in range(0,len(hists)):
        hists[H].SetLineWidth(2)
        if 'nominal' in Fnames[H]:
            hists[H].SetLineColor(1)
            legend.AddEntry(hists[H],'Data','L')
            nomhist = hists[H].Clone()
        if 'up' in Fnames[H]:
            hists[H].SetLineColor(16)
            hists[H].SetLineStyle(9)
            legend.AddEntry(hists[H],'Data (+1 #sigma)','L')
            ratio.append(hists[H].Clone())
        if 'down' in Fnames[H]:
            hists[H].SetLineColor(16)
            hists[H].SetLineStyle(6)
            legend.AddEntry(hists[H],'Data (-1 #sigma)','L')
            ratio.append(hists[H].Clone())
        if 'MC' in Fnames[H]:
            hists[H].SetLineColor(2)
            legend.AddEntry(hists[H],'MC','L')
            ratio.append(hists[H].Clone())
        hists[H].Draw("histSAME")
    hists[0].Draw("AXISSAMEY+")
    hists[0].Draw("AXISSAMEX+")

    legend.Draw("same")
    Lumi = '137.42'
    
    if (year == '2016'):
        Lumi = '35.92'
    if (year == '2017'):
        Lumi = '41.53'
    if (year == '2018'):
        Lumi = '59.97'
    label_cms="CMS Preliminary"
    Label_cms = ROOT.TLatex(0.2,0.92,label_cms)
    Label_cms.SetNDC()
    Label_cms.SetTextFont(61)
    Label_cms.Draw()
    Label_lumi = ROOT.TLatex(0.71,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.Draw("same")
    pad1.Update()

    pad2.cd()
    print ratio
    for rH in ratio: 
        rH.SetTitle("")
        rH.SetMarkerStyle(20)
        rH.SetMarkerSize(1.2)
        rH.GetXaxis().SetTitle('True Number of Interactions')
        rH.GetYaxis().CenterTitle()
        rH.GetXaxis().SetMoreLogLabels()
        rH.GetXaxis().SetNoExponent()
        rH.GetXaxis().SetTitleSize(0.04/0.3)
        rH.GetYaxis().SetTitleSize(0.04/0.3)
        rH.GetXaxis().SetTitleFont(42)
        rH.GetYaxis().SetTitleFont(42)
        rH.GetXaxis().SetTickLength(0.05)
        rH.GetYaxis().SetTickLength(0.05)
        rH.GetXaxis().SetLabelSize(0.115)
        rH.GetYaxis().SetLabelSize(0.089)
        rH.GetXaxis().SetLabelOffset(0.02)
        rH.GetYaxis().SetLabelOffset(0.01)
        rH.GetYaxis().SetTitleOffset(0.42)
        rH.GetXaxis().SetTitleOffset(1.1)
        rH.GetYaxis().SetNdivisions(504)
        rH.GetYaxis().SetRangeUser(0.5,2)
        rH.Divide(nomhist)
        for b in range(rH.GetNbinsX()):
            rH.SetBinContent(b+1,1/rH.GetBinContent(b+1))
        rH.SetStats(ROOT.kFALSE)
        rH.GetYaxis().SetTitle('Data/x')
        rH.Draw('histSAME')
        rH.Draw("AXISSAMEY+")
        rH.Draw("AXISSAMEX+")


    canvas.Print(can + ".png")
    del canvas
    gc.collect()



samples = {}

samples["2017_down"]=[    "PileupHistogram-goldenJSON-13tev-2017-66000ub-99bins.root"]
samples["2017_nominal"]=[    "PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root"]
samples["2017_up"]=[    "PileupHistogram-goldenJSON-13tev-2017-72400ub-99bins.root"]

samples["2016_down"]=[    "PileupHistogram-goldenJSON-13tev-2016-66000ub-75bins.root"]
samples["2016_nominal"]=[    "PileupHistogram-goldenJSON-13tev-2016-69200ub-75bins.root"]
samples["2016_up"]=[    "PileupHistogram-goldenJSON-13tev-2016-72400ub-75bins.root"]

samples["2018_down"]=[    "PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root"]
samples["2018_nominal"]=[    "PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root"]
samples["2018_up"]=[    "PileupHistogram-goldenJSON-13tev-2018-72400ub-99bins.root"]

MC2017_bins = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 
27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 
39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 
51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 
63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 
75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 
87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98
]


#https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py
#https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py
MC2017_bins =[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99]

MC2017_value =[
3.39597497605e-05,6.63688402133e-06,1.39533611284e-05,3.64963078209e-05,6.00872171664e-05,9.33932578027e-05,0.000120591524486,0.000128694546198,0.000361697233219,0.000361796847553,0.000702474896113,0.00133766053707,0.00237817050805,0.00389825605651,0.00594546732588,0.00856825906255,0.0116627396044,0.0148793350787,0.0179897368379,0.0208723871946,0.0232564170641,0.0249826433945,0.0262245860346,0.0272704617569,0.0283301107549,0.0294006137386,0.0303026836965,0.0309692426278,0.0308818046328,0.0310566806228,0.0309692426278,0.0310566806228,0.0310566806228,0.0310566806228,0.0307696426944,0.0300103336052,0.0288355370103,0.0273233309106,0.0264343533951,0.0255453758796,0.0235877272306,0.0215627588047,0.0195825559393,0.0177296309658,0.0160560731931,0.0146022004183,0.0134080690078,0.0129586991411,0.0125093292745,0.0124360740539,0.0123547104433,0.0123953922486,0.0124360740539,0.0124360740539,0.0123547104433,0.0124360740539,0.0123387597772,0.0122414455005,0.011705203844,0.0108187105305,0.00963985508986,0.00827210065136,0.00683770076341,0.00545237697118,0.00420456901556,0.00367513566191,0.00314570230825,0.0022917978982,0.00163221454973,0.00114065309494,0.000784838366118,0.000533204105387,0.000358474034915,0.000238881117601,0.0001984254989,0.000157969880198,0.00010375646169,6.77366175538e-05,4.39850477645e-05,2.84298066026e-05,1.83041729561e-05,1.17473542058e-05,7.51982735129e-06,6.16160108867e-06,4.80337482605e-06,3.06235473369e-06,1.94863396999e-06,1.23726800704e-06,7.83538083774e-07,4.94602064224e-07,3.10989480331e-07,1.94628487765e-07,1.57888581037e-07,1.2114867431e-07,7.49518929908e-08,4.6060444984e-08,2.81008884326e-08,1.70121486128e-08,1.02159894812e-08]

#/Neutrino_E-10_gun/RunIISpring15PrePremix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v2-v2/GEN-SIM-DIGI-RAW
#2016_25ns_Moriond17MC_PoissonOOTPU
#https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py
MC2016_bins =[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74]
MC2016_value =[1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05]
    
MC2018_bins =[    
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
    40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
    60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
    80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99
    ]

MC2018_value =[
    4.695341e-10, 1.206213e-06, 1.162593e-06, 6.118058e-06, 1.626767e-05,
    3.508135e-05, 7.12608e-05, 0.0001400641, 0.0002663403, 0.0004867473,
    0.0008469, 0.001394142, 0.002169081, 0.003198514, 0.004491138,
    0.006036423, 0.007806509, 0.00976048, 0.0118498, 0.01402411,
    0.01623639, 0.01844593, 0.02061956, 0.02273221, 0.02476554,
    0.02670494, 0.02853662, 0.03024538, 0.03181323, 0.03321895,
    0.03443884, 0.035448, 0.03622242, 0.03674106, 0.0369877,
    0.03695224, 0.03663157, 0.03602986, 0.03515857, 0.03403612,
    0.0326868, 0.03113936, 0.02942582, 0.02757999, 0.02563551,
    0.02362497, 0.02158003, 0.01953143, 0.01750863, 0.01553934,
    0.01364905, 0.01186035, 0.01019246, 0.008660705, 0.007275915,
    0.006043917, 0.004965276, 0.004035611, 0.003246373, 0.002585932,
    0.002040746, 0.001596402, 0.001238498, 0.0009533139, 0.0007282885,
    0.000552306, 0.0004158005, 0.0003107302, 0.0002304612, 0.0001696012,
    0.0001238161, 8.96531e-05, 6.438087e-05, 4.585302e-05, 3.23949e-05,
    2.271048e-05, 1.580622e-05, 1.09286e-05, 7.512748e-06, 5.140304e-06,
    3.505254e-06, 2.386437e-06, 1.625859e-06, 1.111865e-06, 7.663272e-07,
    5.350694e-07, 3.808318e-07, 2.781785e-07, 2.098661e-07, 1.642811e-07,
    1.312835e-07, 1.081326e-07, 9.141993e-08, 7.890983e-08, 6.91468e-08,
    6.119019e-08, 5.443693e-08, 4.85036e-08, 4.31486e-08, 3.822112e-08
]


gr2016=[]
gr2017=[]
gr2018=[]
Ngr2016=[]
Ngr2017=[]
Ngr2018=[]

weights={}
for key, value in samples.items():
    if '2016' in key:
        W = []
        f = ROOT.TFile.Open(value[0])
        h = f.Get('pileup')
        h.Scale(1/h.Integral()) 
        gr2016.append(h)
        Ngr2016.append(key)
        for b in range(len(MC2016_bins)):
            W.append(h.GetBinContent(b+1)/MC2016_value[b])
        weights[key] = W
    if '2017' in key:
        W = []
        f = ROOT.TFile.Open(value[0])
        h = f.Get('pileup')
        h.Scale(1/h.Integral())
        gr2017.append(h)
        Ngr2017.append(key)
        for b in range(len(MC2017_bins)-1):
            W.append(h.GetBinContent(b+1)/MC2017_value[b])
        weights[key] = W
    if '2018' in key:
        W = []
        f = ROOT.TFile.Open(value[0])
        h = f.Get('pileup')
        h.Scale(1/h.Integral())
        gr2018.append(h)
        Ngr2018.append(key)
        for b in range(len(MC2018_bins)-1):
            W.append(h.GetBinContent(b+1)/MC2018_value[b])
        weights[key] = W

hcopy = gr2016[0].Clone()
for b in range(hcopy.GetNbinsX()):
    hcopy.SetBinContent(b+1,MC2016_value[b])
gr2016.append(hcopy)
Ngr2016.append('MC_2016')

hcopy = gr2017[0].Clone()
for b in range(hcopy.GetNbinsX()):
    hcopy.SetBinContent(b+1,MC2017_value[b])
gr2017.append(hcopy)
Ngr2017.append('MC_2017')

hcopy = gr2018[0].Clone()
for b in range(hcopy.GetNbinsX()):
    hcopy.SetBinContent(b+1,MC2018_value[b])
gr2018.append(hcopy)
Ngr2018.append('MC_2018')

for key, value in samples.items():
    print 'double pu' + key + '['+str(len(weights[key])) +']={'
    print weights[key]
print sum(MC2016_value)
print sum(MC2017_value)
print sum(MC2018_value)

puDraw(gr2016,Ngr2016,'2016','pu2016')
puDraw(gr2017,Ngr2017,'2017','pu2017')
puDraw(gr2018,Ngr2018,'2018','pu2018')


