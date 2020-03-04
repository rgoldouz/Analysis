# This file skims the data and saves the output to ./tmp
# Do not combine files across runs, otherwise you may get inconsistent TTree structures!
# Doing things file by file is the safest way to avoid this problem, and comes at almost
# no extra cost.
# You can copy and paste json sources directly from https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/

import os
import ROOT
path = []

#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_v3/191130_081246/0000/')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_v3_ext1/191130_081323/0000/')

#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_v3_ext1/191130_075804/0000')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_v3_ext2/191130_081130/0000')

#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_TTTo2L2Nu/191130_081559/0000')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WWTo2L2Nu_13TeV-powheg/crab_WWTo2L2Nu/191130_090952/0000')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ZZTo2L2Nu_13TeV_powheg_pythia8/crab_ZZTo2L2Nu/191130_091725/0000')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_antitop_v3_ext1/191130_092526/0000')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_v3_ext1/191130_093520/0000')

#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ZZTo4L_13TeV_powheg_pythia8/crab_ZZTo4L/191130_092020/0000')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo2L2Q/191130_091455/0000/')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_v3_ext1_v1/191201_062001/0000/')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_v3_ext2_v2/191130_113603/0000/')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_v3_v2/191130_115257/0000/')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ/191130_110848/0000')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuN/191130_110706/0000/')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_v3_ext1_v2/191201_071154/0000/')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_v3_ext3_v1/191201_071416/0000/')

#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu_v3_ext1_v2/191201_070704/0000/')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu/191130_110349/0000')
#path.append('/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToQQ/191130_110154/0000/')
#path.append('/pnfs/iihe/cms/store/user/rgoldouz/TOPptSamples/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_TT_Mtt_1000toInf/200226_084805/0000')
#path.append('/pnfs/iihe/cms/store/user/rgoldouz/TOPptSamples/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_TT/200226_084823/0000')
path.append('/pnfs/iihe/cms/store/user/rgoldouz/TOPptSamples/TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_TT_Mtt_700to1000/200226_084843/0000')

nEventsraw = 0
neventsweight = 0
nEventsStored = 0
nEventsiihe = 0
for a in path:
    filenames = os.listdir(a)
    for fname in filenames:
        filename = a + '/' + fname
        print fname 
        if 'fail' in fname:
            continue
        f = ROOT.TFile.Open(filename)
        
        if not f:
            print 'rm -rf '+fname
        tree_in = f.Get('IIHEAnalysis')
        tree_meta = f.Get('meta')
        nEventsiihe += tree_in.GetEntries()
        tree_meta.GetEntry(0)    
        print tree_meta.nEventsRaw
        nEventsraw += tree_meta.nEventsRaw
        nEventsStored += tree_meta.nEventsStored
        neventsweight += tree_meta.mc_nEventsWeighted
        f.Close()
print 'nEventsraw %d   '%(nEventsraw)
print 'neventsweight %d   '%(neventsweight)
print 'nEventsStored %d   '%(nEventsStored)
print 'nEventsiihe %d   '%(nEventsiihe)
