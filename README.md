# Charged Lepton Flavour Violation Analysis  

Author: Reza Goldouzian

Modified by Jack Li for trilepton selection

This frame work can be compiled standalone, although it depends on ROOT&&Boost libraries. Make sure you have "boost" installed when you try to compile this on your own machine. 

## To compile & run 

```sh
git clone -b Trilepton_Selection https://github.com/Jingyan95/Analysis.git TopLFV
cd TopLFV
make all
./RunAll

it takes 1~2 min to generate test.root file where you can find all the interesting distributions.  
```

## To write & submit jobs (Condor) 

```sh
cd TopLFV/bin
python makeJobs.py
python submitJobs.py

Output files can be found at TopLFV/hists/2017/
Make sure replacing the directories that are hard-coded in the two python files above with your own directories on lxplus. 
```


