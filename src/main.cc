#include "../include/MyAnalysis.h"
int main(){
    TChain* ch    = new TChain("Events");  //"IIHEAnalysis") ;
    //ch ->Add("/Users/jingyanli/eclipse-workspace/ROOT/outfile_2017_TT_vector_emutc.root");
    ch ->Add("root://cms-xrd-global.cern.ch///store/user/piedavid/topNanoAOD/v6-1-1/2017/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/TopNanoAODv6-1-1_2017/200613_065556/0000/tree_9.root");
    	//"/eos/cms/store/user/skinnari/TopLFV/outfile_2017_ST_vector_emutu/outfile_2718.root");
    MyAnalysis t1(ch);
    t1.Loop("testnano.root", "mc" , "" , "2017" , "" , 2.06 , 41.53 , 494000);
    return 0;
}
