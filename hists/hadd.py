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
    year = value[3]
    os.system('rm '+ year + '/' +key + '.root')
    nf = 40
    if value[1]=='data':
        addedFilesData[year].append( year + '/' + key + '.root ')
    if ('TTTo2L2Nu' not in key) and ('DY' not in  key) and ('ST_tW' not in  key) and ('WJetsToLNu' not in  key) and (value[1]=='mc'):
        addedFilesMc[year].append(  year + '/' + key + '.root ')
    hadd='hadd ' + year + '/' + key + '.root '
    for idx, S in enumerate(value[0]):
        if value[1]=='data':
            nf = 250
        for subdir, dirs, files in os.walk(S):
            sequance = [files[i:i+nf] for i in range(0,len(files),nf)]
            for num,  seq in enumerate(sequance):
                hadd +=  year + '/' + key +'_' + str(idx) +'_' + str(num) + '.root '

    os.system(hadd)

os.system('rm *_data.root')
os.system('rm *_others.root')
hadddata_2016 ='hadd 2016_data' + '.root ' + ' '.join(addedFilesData['2016'])
hadddata_2017 ='hadd 2017_data' + '.root ' + ' '.join(addedFilesData['2017'])
hadddata_2018 ='hadd 2018_data' + '.root ' + ' '.join(addedFilesData['2018'])

haddmc_2016 ='hadd 2016_others' + '.root ' + ' '.join(addedFilesMc['2016'])
haddmc_2017 ='hadd 2017_others' + '.root ' + ' '.join(addedFilesMc['2017'])
haddmc_2018 ='hadd 2018_others' + '.root ' + ' '.join(addedFilesMc['2018'])


os.system(haddmc_2016)
os.system(haddmc_2017)
os.system(haddmc_2018)

os.system(hadddata_2016)
os.system(hadddata_2017)
os.system(hadddata_2018)

os.system('rm *_DY.root')
os.system('hadd 2016_DY.root 2016/2016_DYM50.root 2016/2016_DYM10to50.root')
os.system('hadd 2017_DY.root 2017/2017_DYM50.root 2017/2017_DYM10to50.root')
os.system('hadd 2018_DY.root 2018/2018_DYM50.root 2018/2018_DYM10to50.root')

Y = ['2016','2017','2018']
Sam = ['TTTo2L2Nu','ST_tW','WJetsToLNu']

for y in Y:
    for s in Sam:
        txt = 'cp ' + y + '/' + y + '_' + s + '.root .'
        os.system(txt)

allSamples = ['TTTo2L2Nu','ST_tW','WJetsToLNu','data', 'others', 'DY']

for s in allSamples:
    os.system('rm All_' + s + '.root ')
    haddall='hadd All_' + s + '.root '
    for y in Y:
        haddall +=  y + '_' + s + '.root '
    os.system(haddall) 
       
