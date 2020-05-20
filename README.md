# Top LFV Analysis  

Author: Reza Goldouzian

Modified by Northeastern Group for trilepton selection

This framework can be compiled standalone, although it depends on ROOT&&BOOST libraries. Make sure you have "BOOST" installed when you try to compile this on your own machine. 

## I. File Lists

bin/Files_2016.py                   Python file to store addresses of 2016 samples                  
bin/Files_2017.py                   ... 2017 ...
bin/Files_2018.py                   ... 2018 ...
bin/Files_2017_A.py                 Addresses of samples for Top Trilepton analysis
bin/makeGFAL-COPYtxtfiles.py        Python function to make gfal-copy files
bin/GFAL-COPY-ALL.py                Python function to run makeGFAL-COPYtxtfiles.py
bin/makeJobs.py                     Python function for writing condor jobs
bin/submitJobs.py                   Python function for submitting condor jobs
bin/input/                          Duplicate of input/ directory
bin/Jobs/                           Directory where condor job files live
helper/                             Directory where utility functions live
hists/2016/                         Directory to store output file for 2016 data/mc
hists/2017/                         ... 2017 ...
hists/2018/                         ... 2018 ...
hists/hadd.py                       Utility function for merging root files
include/                            Directory where header files live
input/                              Directory where input files live 
plot/                               Directory where plotting functions live
src/BTagCalibrationStandalone.cc    B-tagging calibration file
src/MyAnalysis.cc                   Main analysis file 
src/PU_reWeighting.cc               Pile-Up reweighting file 
src/RoccoR.cc                       Rochester correction file 
src/jet_candidate.cc                Object class for jet
src/lepton_candidate.cc             Object class for lepton
src/main.cc                         Testing file 
Makefile                            Tool for compiling 

## II. To compile & run 

```sh
git clone -b Trilepton_Selection https://github.com/Jingyan95/Analysis.git TopLFV
cd TopLFV
make all

This generates an executable "RunAll" to run src/main.cc which is a testing file with signal samples as input.

./RunAll

It takes 1~2 min to produce test.root file where you can find all the interesting distributions.  
```

## III. To write & submit jobs (Condor) 

The first thing you need to do is to modify these two lines https://github.com/Jingyan95/Analysis/blob/Trilepton_Selection/bin/makeJobs.py#L40-L41 , as well as this https://github.com/Jingyan95/Analysis/blob/Trilepton_Selection/bin/submitJobs.py#L8 . Make sure replacing these with your own paths on lxplus. 

```sh
cd TopLFV/bin
python makeJobs.py
python submitJobs.py

Output files can be found at TopLFV/hists/2017/ 
```

### makeJobs.py

This function will walk through addresses stored in Files_2017_A.py and write two types of files needed for submitting condor jobs under bin/Jobs/. 
The first type of files is .C file which we use to call our analysis code. This is the block https://github.com/Jingyan95/Analysis/blob/Trilepton_Selection/bin/makeJobs.py#L59-L70 where a .C file is written. An example of .C files:

2017_LFVStVecC_0_0.C

```sh
#include "MyAnalysis.h"
main(){
    TChain* ch    = new TChain("IIHEAnalysis") ;
    ch ->Add("/eos/user/a/asparker/TopLFV/2017_LFVStVecC/outfile_1467.root");
    ch ->Add("/eos/user/a/asparker/TopLFV/2017_LFVStVecC/outfile_1476.root");
    MyAnalysis t1(ch);
    t1.Loop("/afs/cern.ch/user/j/jingyan/TopLFV/hists/2017/2017_LFVStVecC_0_0.root", "mc" , "" , "2017" , "" , 0.0512 , 41.53 , 500000);
}
```

The second type of file is .sh file which is an executable where we write some commands to tell the machine what to do once the job is running. These are relevant lines https://github.com/Jingyan95/Analysis/blob/Trilepton_Selection/bin/makeJobs.py#L73-L83 . Here is an example of .sh file:

2017_LFVStVecC_0_0.sh

```sh
#!/bin/bash
cd /afs/cern.ch/user/j/jingyan/TopLFV/bin
g++ -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2  -I./../include   -pthread -std=c++11 -m64 -I/usr/include/root -ldl  -o 2017_LFVStVecC_0_0 Jobs/2017_LFVStVecC_0_0.C ../lib/main.so -L/usr/lib64/root -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic  -lMinuit -lTreePlayer
./2017_LFVStVecC_0_0
FILE=/afs/cern.ch/user/j/jingyan/TopLFV/hists/2017/2017_LFVStVecC_0_0.root
if [ -f "$FILE" ]; then
    rm  2017_LFVStVecC_0_0
fi 
```

As you may have noticed, this file will tell to the machine to go to the right directory (TopLFV/bin), compile the .C file and generate an executable (2017_LFVStVecC_0_0). Then this executable will be executed and deleted after it finishes running. 

### submitJobs.py

Now we have both .C file and .sh files in place, but we are not yet ready to submit jobs to clusters. A configuration file (.sub) is needed which can be written by submitJobs.py function https://github.com/Jingyan95/Analysis/blob/Trilepton_Selection/bin/submitJobs.py#L39-L49 to set up all the configurations. An example of .sub files:

submit01.sub

```sh
universe = vanilla
arguments = "$(argument)"
output = submit01.out
error = submit01.err
log = submit01.log
+JobFlavour = "tomorrow"
queue
```

To actually submit jobs we use the "condor_submit" command to submit the .sub file along with the .sh executable generated by makeJobs.py https://github.com/Jingyan95/Analysis/blob/Trilepton_Selection/bin/submitJobs.py#L64 .
To remove jobs, we can use condor_rm command. For example:

```sh
condor_rm jingyan

This removes all the jobs submitted by user jingyan.

condor_rm 4698914.0

This removes the job with a certain ID.
```

More on Condor can be found at http://information-technology.web.cern.ch/services/fe/lxbatch/howto/quickstart-guide-htcondor


