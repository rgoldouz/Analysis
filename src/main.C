#include "MyAnalysis.h"

main(){
    TChain* ch    = new TChain("IIHEAnalysis") ;
//    ch   ->Add("/group/HEEP/SYNCH/MC2016_ttTolnu.root") ;
//    ch   ->Add("/group/HEEP/SYNCH/MuonEG_runH.root");    
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_781.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_736.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_655.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_574.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_529.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_493.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_448.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_367.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_286.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_944.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_863.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_818.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_782.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_737.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_656.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_575.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_494.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_449.root");
    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/outfile_368.root");
    MyAnalysis t1(ch);

    t1.Loop("/user/rgoldouz/NewAnalysis2020/Analysis/hists/test.root", "data","SingleElectron","2016", "C",1,1,1);
}




