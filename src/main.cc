#include "MyAnalysis.h"
int main(){
    TChain* ch    = new TChain("IIHEAnalysis") ;
    //ch ->Add("/Users/jingyanli/eclipse-workspace/ROOT/outfile_2017_TT_vector_emutc.root");
    ch ->Add("/eos/cms/store/user/skinnari/TopLFV/outfile_2017_ST_vector_emutu/outfile_2718.root");
    MyAnalysis t1(ch);
    t1.Loop("test.root", "mc" , "" , "2017" , "" , 2.06 , 41.53 , 494000);
    return 0;
}
