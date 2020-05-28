# Top LFV Analysis  

Author: Reza Goldouzian

Modified by Northeastern Group for trilepton selection

This framework can be compiled standalone, although it depends on ROOT&&BOOST libraries. Make sure you have "BOOST" installed when you try to compile this on your own machine. 

## I. File Lists

<table border="0">
 <tr>
    <td><b style="font-size:30px">Files</b></td>
    <td><b style="font-size:30px">Description</b></td>
 </tr>
 <tr>
    <td>bin/Files_2016.py</td>
    <td>Python file to store addresses of 2016 samples</td>
 </tr>
 <tr>
    <td>bin/Files_2017.py</td>
    <td>Python file to store addresses of 2017 samples</td>
 </tr>
 <tr>
    <td>bin/Files_2018.py</td>
    <td>Python file to store addresses of 2018 samples</td>
 </tr>
 <tr>
    <td>bin/Files_2017_A.py</td>
    <td>Addresses of samples for Top Trilepton analysis</td>
 </tr>
 <tr>
    <td>bin/makeGFAL-COPYtxtfiles.py</td>
    <td>Python function to make gfal-copy files</td>
 </tr>
 <tr>
    <td>bin/GFAL-COPY-ALL.py</td>
    <td>Python function to run makeGFAL-COPYtxtfiles.py</td>
 </tr>
 <tr>
    <td>bin/makeJobs.py</td>
    <td>Python function for writing condor jobs</td>
 </tr>
 <tr>
    <td>bin/submitJobs.py</td>
    <td>Python function for submitting condor jobs</td>
 </tr>
 <tr>
    <td>bin/input/</td>
    <td>Duplicate of input/ directory</td>
 </tr>
 <tr>
    <td>bin/Jobs/</td>
    <td>Directory where condor job files live</td>
 </tr>
 <tr>
    <td>helper/</td>
    <td>Directory where utility functions live</td>
 </tr>
 <tr>
    <td>hists/2016/</td>
    <td>Directory to store output file for 2016 data/mc</td>
 </tr>
 <tr>
    <td>hists/2017/</td>
    <td>Directory to store output file for 2017 data/mc</td>
 </tr>
 <tr>
    <td>hists/2018/</td>
    <td>Directory to store output file for 2018 data/mc</td>
 </tr>
 <tr>
    <td>hists/hadd.py</td>
    <td>Utility function for merging root files</td>
 </tr>
 <tr>
    <td>include/</td>
    <td>Directory where header files live</td>
 </tr>
 <tr>
    <td>input/</td>
    <td>Directory where input files live</td>
 </tr>
 <tr>
    <td>plot/</td>
    <td>Directory where plotting functions live</td>
 </tr>
 <tr>
    <td>src/BTagCalibrationStandalone.cc</td>
    <td>B-tagging calibration file</td>
 </tr>
 <tr>
    <td>src/MyAnalysis.cc</td>
    <td>Main analysis file</td>
 </tr>
 <tr>
    <td>src/PU_reWeighting.cc</td>
    <td>Pile-Up reweighting file</td>
 </tr>
 <tr>
    <td>src/RoccoR.cc</td>
    <td>Rochester correction file</td>
 </tr>
 <tr>
    <td>src/jet_candidate.cc</td>
    <td>Object class for jet</td>
 </tr>
 <tr>
    <td>src/lepton_candidate.cc</td>
    <td>Object class for lepton</td>
 </tr>
 <tr>
    <td>src/main.cc</td>
    <td>Testing file</td>
 </tr>
 <tr>
    <td>Makefile</td>
    <td>Tool for compiling</td>
 </tr>
</table>

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

The first thing you need to do is to modify the --l option below, replacing these with your own paths on lxplus.The --n option tells the name in a given sample you want to run, the one below will give you the signal samples or 2017 would give you all 2017 files, see Files_2017_A.py for the keys to the dictionary in order to choose a different --n option to suit your needs. 

```sh
cd TopLFV/bin
python makeJobs.py --l /afs/cern.ch/user/a/asparker/public/LFVTopCode_MyFork/TopLFV/ --n SMEFTfr
python submitJobs.py --l /afs/cern.ch/user/a/asparker/public/LFVTopCode_MyFork/TopLFV/ --n SMEFTfr



Output files can be found at TopLFV/hists/2017/ 
```

### makeJobs.py

This function will walk through addresses stored in Files_2017_A.py and write two types of files needed for submitting condor jobs under bin/Jobs/. 
The first type of files is .C file which we use to call our analysis code. This is the block https://github.com/Jingyan95/Analysis/blob/Trilepton_Selection/bin/makeJobs.py#L78-L89 where a .C file is written. An example of .C files:

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

The second type of file is .sh file which is an executable where we write some commands to tell the machine what to do once the job is running. These are relevant lines https://github.com/Jingyan95/Analysis/blob/Trilepton_Selection/bin/makeJobs.py#L92-L102 . Here is an example of .sh file:

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

Now we have both .C file and .sh files in place, but we are not yet ready to submit jobs to clusters. A configuration file (.sub) is needed which can be written by submitJobs.py function https://github.com/Jingyan95/Analysis/blob/Trilepton_Selection/bin/submitJobs.py#L55-L65 to set up all the configurations. An example of .sub files:

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

To actually submit jobs we use the "condor_submit" command to submit the .sub file along with the .sh executable generated by makeJobs.py https://github.com/Jingyan95/Analysis/blob/Trilepton_Selection/bin/submitJobs.py#L80 . We can use condor_q to monitor the jobs. To remove jobs, we can use condor_rm command. For example:

```sh
condor_rm jingyan

This removes all the jobs submitted by user jingyan.

condor_rm 4698914.0

This removes the job with a certain ID.
```

More on Condor can be found at http://information-technology.web.cern.ch/services/fe/lxbatch/howto/quickstart-guide-htcondor


