#include "../include/MyAnalysis.h"
int main(){
    TChain* ch    = new TChain("IIHEAnalysis") ;
    ch ->Add("/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_TT_clequ1_emutu/outfile_9073.root");
    MyAnalysis * t1 = new MyAnalysis(ch);

//    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runD/191203_070404/0000/outfile_85.root");
//    MyAnalysis t1(ch);
//    t1.Loop("2016_D_MuonEG", "2016_D_MuonEG_0_0.root", "data" , "MuonEG" , "2016" , "D" , 1 , 1 , 1);
    t1->Loop("LFVTt.root","LFVTt",  "mc" , "" , "2016" , "" , 87.31 , 35.92 , 67312164);
//    ch ->Add("/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ1_emutu/outfile_12096.root");
//    ch ->Add("/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ1_emutu/outfile_12146.root");
//    MyAnalysis t1(ch);
//    t1.Loop("2017_LFVStScalarU_0_0.root","2017_LFVStScalarU", "mc" , "" , "2017" , "" , 0.102 , 41.53 , 500000);
delete t1;



}
