import sys
import os
import subprocess
import readline
import string

sys.path.append('/user/rgoldouz/NewAnalysis2020/Analysis/bin')

import Files_2016
import Files_2017
import Files_2018
SAMPLES = {}
mc_2016 = True
data_2016 = True
mc_2017 = True
data_2017 = True
mc_2018 = True
data_2018 = True

if mc_2016:
    SAMPLES.update(Files_2016.mc2016_samples)
if data_2016:
    SAMPLES.update(Files_2016.data2016_samples)
if mc_2017:
    SAMPLES.update(Files_2017.mc2017_samples)
if data_2017:
    SAMPLES.update(Files_2017.data2017_samples)
if mc_2018:
    SAMPLES.update(Files_2018.mc2018_samples)
if data_2018:
    SAMPLES.update(Files_2018.data2018_samples)


addedFilesData = {"2016": [], "2017": [], "2018": []}
addedFilesMc = {"2016": [], "2017": [], "2018": []}

for key, value in SAMPLES.items():
#    if '2016' not in key:
#        continue
    year = value[3]
    nf = 50
    if 'TTTo2L2Nu' in key or 'tw' in key or 'DY' in key or 'TTsys' in key:
        nf = 25
    if value[1]=='data':
        nf = 175
    for idx, S in enumerate(value[0]):
        for subdir, dirs, files in os.walk(S):
            sequance = [files[i:i+nf] for i in range(0,len(files),nf)]
            for num,  seq in enumerate(sequance):
                if os.path.isfile('/user/rgoldouz/NewAnalysis2020/Analysis/hists/' + year + '/' + key +'_' + str(idx) +'_' + str(num) + '.root'):
                    continue
                f = key +'_' + str(idx) +'_' + str(num)
                qsub = "qsub -q localgrid  -o STDOUT/" + f + ".stdout -e STDERR/" + f + ".stderr Jobs/" + f + '.sh'
                exit_status = subprocess.call(qsub, shell=True)
            break
