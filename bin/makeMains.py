import sys
import os
import subprocess
import readline
import string

MCsamples = {}
#MCsamples ['DYM10to50'] = ['address', 'data/mc','year', 'run', 'cross section','lumi','Neventsraw']
MCsamples['DYM10to50_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_v3/191130_081246/0000/'], 'mc','2016', '','18610','35.92','78843820']
MCsamples['DYM10to50_2'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_v3_ext1/191130_081323/0000/'], 'mc','2016', '','18610','35.92','78843820']


MCsamples['DYM50_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_v3_ext1/191130_075804/0000/'], 'mc','2016', '','5765.4','35.92','146280395']
MCsamples['DYM50_2'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_v3_ext2/191130_081130/0000/'], 'mc','2016', '','5765.4','35.92','146280395']

MCsamples['TTTo2L2Nu_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_TTTo2L2Nu/191130_081559/0000/'], 'mc','2016', '','87.31','35.92','79140880']
MCsamples['ST_tW_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_antitop_v3_ext1/191130_092526/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_v3_ext1/191130_093520/0000/'], 'mc','2016', '','38.94','35.92','6368069']

MCsamples['WWTo2L2Nu_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WWTo2L2Nu_13TeV-powheg/crab_WWTo2L2Nu/191130_090952/0000/'], 'mc','2016', '','12.178','35.92','1999000']
MCsamples['ZZTo2L2Nu_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ZZTo2L2Nu_13TeV_powheg_pythia8/crab_ZZTo2L2Nu/191130_091725/0000/'], 'mc','2016', '','0.564','35.92','8931750']
MCsamples['ZZTo4L_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ZZTo4L_13TeV_powheg_pythia8/crab_ZZTo4L/191130_092020/0000/'], 'mc','2016', '','1.212','35.92','6669988']
MCsamples['WZTo2L2Q_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo2L2Q/191130_091455/0000/'], 'mc','2016', '','5.595','35.92','15879472']
MCsamples['WZTo3LNu_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_v3_ext1_v1/191201_062001/0000/'], 'mc','2016', '','4.43','35.92','18000000']

MCsamples['WJetsToLNu_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_v3_ext2_v2/191130_113603/0000/'], 'mc','2016', '','61526.7','35.92','86916455']
MCsamples['WJetsToLNu_2'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_v3_v2/191130_115257/0000/'], 'mc','2016', '','61526.7','35.92','86916455']


MCsamples['TTWJetsToQQ_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToQQ/191130_110154/0000/'], 'mc','2016', '','0.4062','35.92','430310']
MCsamples['TTWJetsToLNu_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu_v3_ext1_v2/191201_070704/0000/','/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu/191130_110349/0000/'], 'mc','2016', '','0.2043','35.92','2716249']
MCsamples['TTZToLLNuN_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuN/191130_110706/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_v3_ext1_v2/191201_071154/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_v3_ext3_v1/191201_071416/0000/'], 'mc','2016', '','0.2529','35.92','6420825']
MCsamples['TTZToQQ_1'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ/191130_110848/0000/'], 'mc','2016', '','0.5297','35.92','351164']


rootlib1 = subprocess.check_output("root-config --cflags", shell=True)
rootlib11="".join([s for s in rootlib1.strip().splitlines(True) if s.strip()])
rootlib2 = subprocess.check_output("root-config --glibs", shell=True)
rootlib22="".join([s for s in rootlib2.strip().splitlines(True) if s.strip()])
print 'llll'
print rootlib11
print rootlib22

for key, value in MCsamples.items():
    text = ''
    text += '    TChain* ch    = new TChain("IIHEAnalysis") ;\n'

    for S  in value[0]:
        text += '    ch ->Add("' + S+ '*.root");\n'
 
    text += '    MyAnalysis t1(ch);\n'
    text += '    t1.Loop("ntuples/' + key + '.root", "' + value[1] + '" , "'+ value[2] + '" , "'+ value[3] + '" , '+ value[4] + ' , ' + value[5] + ' , '+ value[6] + ');\n'
    

    SHNAME = key + '.C'
    SHFILE='#include "MyAnalysis.h"\n' +\
    'main(){\n' +\
    text +\
    '}'

    open(SHNAME, 'wt').write(SHFILE)

    os.system('g++ -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2  -I./../include   '+ rootlib11 +' -ldl  -o ' + key + ' ' + SHNAME+ ' ../lib/main.so ' + rootlib22 + '  -lMinuit -lMinuit2 -lTreePlayer -lGenVector')
#    os.system("chmod +x "+SHNAME)
