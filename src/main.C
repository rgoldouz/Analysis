#include "MyAnalysis.h"

main(){
    TChain* ch_DY    = new TChain("IIHEAnalysis") ;
//    ch_DY   ->Add("/group/HEEP/SYNCH/MC2016_ttTolnu.root") ;
    ch_DY   ->Add("/group/HEEP/TTTo2L2Nu_2017.root");    
    MyAnalysis t1(ch_DY);
    t1.Loop("/user/rgoldouz/NewAnalysis2020/Analysis/hists/test.root", "mc" , "","2017" , "" , 5765.4 , 35.92 , 146280395);
}




