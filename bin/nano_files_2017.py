import sys
import os
import subprocess
import readline
import string

data2017_samples = {}
mc2017_samples = {}
#data2017_samples ['DYM10to50'] = ['address', 'data/mc','dataset','year', 'run', 'cross section','lumi','Neventsraw']                        

mc2017_samples['2017_DYM10to50'] = [['root://cms-xrd-global.cern.ch////store/user/piedavid/topNanoAOD/v6-1-1/2017/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/TopNanoAODv6-1-1_2017/200613_065556/0000/'], 'mc','','2017', '','18610','41.53','39521230'] ## '78,994,955' ## <-  is event number from miniaod, we must be missing an extension

mc2017_samples['2017_DYM50'] = [['root://cms-xrd-global.cern.ch////store/user/piedavid/topNanoAOD/v6-1-1/2017/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/TopNanoAODv6-1-1_2017/200613_070408/0000/'], 'mc','','2017', '','5765.4','41.53','182217609'] ## '84,597,528' ## now we have ore events ?  182,217,609

## example nanoAOD datset on DAS where I got number of events from
## https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fphys03&input=dataset+%3D+%2FTTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8%2Fpalencia-TopNanoAODv6-1-1_2018-0d1d4920f08f56d048ece029b873a2cc%2FUSER
mc2017_samples['2017_TTTo2L2Nu'] = [['root://cms-xrd-global.cern.ch////store/group/phys_top/topNanoAOD/v6-1-1/2018/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/TopNanoAODv6-1-1_2018/200610_155014/0000/'], 'mc','','2017', '','87.31','41.53', '64310000' ] ##'8926992'] what is going on with number of events?
