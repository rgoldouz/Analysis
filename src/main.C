#include "../include/MyAnalysis.h"
int main(){
    TChain* ch    = new TChain("IIHEAnalysis") ;
//    ch ->Add("/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_TuneCP5_PSweights_2016/200311_140933/0000/outfile_434.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runD/191203_070404/0000/outfile_85.root");
    MyAnalysis t1(ch);
    t1.Loop("2016_D_MuonEG_0_0.root", "data" , "MuonEG" , "2016" , "D" , 1 , 1 , 1);
//    t1.Loop("TTTo2L2Nu.root",  "mc" , "" , "2016" , "" , 87.31 , 35.92 , 67312164);
}
