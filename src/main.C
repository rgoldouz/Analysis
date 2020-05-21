#include "MyAnalysis.h"
main(){
    TChain* ch    = new TChain("IIHEAnalysis") ;
    ch ->Add("/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_vector_emutu/outfile_2562.root");
    MyAnalysis t1(ch);
    t1.Loop("test.root", "mc" , "" , "2018" , "" , 2.06 , 41.53 , 494000);
}
