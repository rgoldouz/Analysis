import sys
import os
import subprocess
import readline
import string

data2016_samples = {}
mc2016_samples = {}
#data2016_samples ['DYM10to50'] = ['address', 'data/mc','dataset','year', 'run', 'cross section','lumi','Neventsraw']

#2016 MC
mc2016_samples['2016_DYM10to50'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_v3/191130_081246/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_v3_ext1/191130_081323/0000/'], 'mc','','2016', '','18610','35.92','78843820']
mc2016_samples['2016_DYM50'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50-amcFXFX/200310_182844/0000/'], 'mc','','2016', '','5765.4','35.92','80924255']
mc2016_samples['2016_TTTo2L2Nu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_TuneCP5_PSweights_2016/200311_140933/0000/'], 'mc','','2016', '','87.31','35.92','67312164']
mc2016_samples['2016_ST_tW'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_v3_ext1/191130_093520/0000/'], 'mc','','2016', '','19.47','35.92','3256650']
mc2016_samples['2016_ST_atW'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_antitop_v3_ext1/191130_092526/0000/'], 'mc','','2016', '','19.47','35.92','3111419']
mc2016_samples['2016_WWTo2L2Nu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WWTo2L2Nu_13TeV-powheg/crab_WWTo2L2Nu/191130_090952/0000/'], 'mc','','2016', '','12.178','35.92','1999000']
mc2016_samples['2016_ZZTo2L2Nu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ZZTo2L2Nu_13TeV_powheg_pythia8/crab_ZZTo2L2Nu/191130_091725/0000/'], 'mc','','2016', '','0.564','35.92','8931750']
mc2016_samples['2016_ZZTo4L'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ZZTo4L_13TeV_powheg_pythia8/crab_ZZTo4L/191130_092020/0000/'], 'mc','','2016', '','1.212','35.92','6669988']
mc2016_samples['2016_WZTo2L2Q'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo2L2Q/191130_091455/0000/'], 'mc','','2016', '','5.595','35.92','15879472']
mc2016_samples['2016_WZTo3LNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_v3_ext1_v1/191201_062001/0000/'], 'mc','','2016', '','4.43','35.92','18000000']
mc2016_samples['2016_WJetsToLNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amc_v3_v1/191130_144027/0000/','/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amc_v3_ext2_v1/191130_144155/0000/','/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amc_v3_ext2_v1/191130_144155/0001/'], 'mc','','2016', '','61526.7','35.92','178497178']
mc2016_samples['2016_TTWJetsToQQ'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToQQ/191130_110154/0000/'], 'mc','','2016', '','0.4062','35.92','430310']
mc2016_samples['2016_TTWJetsToLNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu_v3_ext1_v2/191201_070704/0000/','/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu/191130_110349/0000/'], 'mc','','2016', '','0.2043','35.92','2716249']
mc2016_samples['2016_TTZToLLNuNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuN/191130_110706/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_v3_ext1_v2/191201_071154/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_v3_ext3_v1/191201_071416/0000/'], 'mc','','2016', '','0.2529','35.92','6420825']
mc2016_samples['2016_TTZToQQ'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ/191130_110848/0000/'], 'mc','','2016', '','0.5297','35.92','351164']

mc2016_samples['2016_TTsys_hdampUP'] = [['/pnfs/iihe/cms/store/user/schenara/SYS_RunII_2016/TTTo2L2Nu_hdampUP_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_hdampUP_v3_v1/200311_144425/0000/', '/pnfs/iihe/cms/store/user/schenara/SYS_RunII_2016/TTTo2L2Nu_hdampUP_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_hdampUP_v3_ext1/200311_144331/0000/'], 'mc','','2016', '','87.31','35.92','14821776']
mc2016_samples['2016_TTsys_hdampDOWN'] = [['/pnfs/iihe/cms/store/user/schenara/SYS_RunII_2016/TTTo2L2Nu_hdampDOWN_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_hdampDOWN_v3_v1/200311_144843/0000/','/pnfs/iihe/cms/store/user/schenara/SYS_RunII_2016/TTTo2L2Nu_hdampDOWN_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_hdampDOWN_v3_ext1/200311_144755/0000/'], 'mc','','2016', '','87.31','35.92','14663032']
mc2016_samples['2016_TTsys_TuneCP5up'] = [['/pnfs/iihe/cms/store/user/schenara/SYS_RunII_2016/TTTo2L2Nu_TuneCP5up_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_TuneCP5up_v3_v1/200311_145034/0000/','/pnfs/iihe/cms/store/user/schenara/SYS_RunII_2016/TTTo2L2Nu_TuneCP5up_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_TuneCP5up_v3_ext1/200311_144936/0000/'], 'mc','','2016', '','87.31','35.92','14719314']
mc2016_samples['2016_TTsys_TuneCP5down'] = [['/pnfs/iihe/cms/store/user/schenara/SYS_RunII_2016/TTTo2L2Nu_TuneCP5down_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_TuneCP5down_v3_v1/200311_145235/0000/','/pnfs/iihe/cms/store/user/schenara/SYS_RunII_2016/TTTo2L2Nu_TuneCP5down_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_TuneCP5down_v3_ext1/200311_145122/0000/'], 'mc','','2016', '','87.31','35.92','14194688']
mc2016_samples['2016_TTsys_CR1QCDbased'] = [['/pnfs/iihe/cms/store/user/schenara/SYS_RunII_2016/TTTo2L2Nu_TuneCP5CR1_QCDbased_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_QCDbased/200311_145409/0000/'], 'mc','','2016', '','87.31','35.92','13849010']
mc2016_samples['2016_TTsys_CR2QCDbased'] = [['/pnfs/iihe/cms/store/user/schenara/SYS_RunII_2016/TTTo2L2Nu_TuneCP5CR2_GluonMove_PSweights_13TeV-powheg-pythia8/crab_TTTo2L2Nu_TuneCP5CR2_GluonMove_2016/200514_062950/0000/'], 'mc','','2016', '','87.31','35.92','14867026']
mc2016_samples['2016_TTsys_CRerdON'] = [['/pnfs/iihe/cms/store/user/schenara/SYS_RunII_2016/TTTo2L2Nu_TuneCP5_PSweights_erdON_13TeV-powheg-pythia8/crab_TTTo2L2Nu_erdON/200311_145319/0000/'], 'mc','','2016', '','87.31','35.92','9647222']

mc2016_samples['2016_LFVStVecC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_ST_vector_emutc/'], 'mc','','2016', '','0.0512','35.92','284000']
mc2016_samples['2016_LFVStVecU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_ST_vector_emutu/'], 'mc','','2016', '','0.515','35.92','492000']
mc2016_samples['2016_LFVTtVecC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_TT_vector_emutc/'], 'mc','','2016', '','0.032','35.92','478000']
mc2016_samples['2016_LFVTtVecU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_TT_vector_emutu/'], 'mc','','2016', '','0.032','35.92','488000']

mc2016_samples['2016_LFVStScalarC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ1_emutc/'], 'mc','','2016', '','0.008' ,'35.92','182000']
mc2016_samples['2016_LFVStScalarU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ1_emutu/'], 'mc','','2016', '','0.102' ,'35.92','302000']
mc2016_samples['2016_LFVTtScalarC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_TT_clequ1_emutc/'], 'mc','','2016', '','0.004' ,'35.92','466000']
mc2016_samples['2016_LFVTtScalarU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_TT_clequ1_emutu/'], 'mc','','2016', '','0.004' ,'35.92','322000']

mc2016_samples['2016_LFVStTensorC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ3_emutc/'], 'mc','','2016', '','0.187' ,'35.92','138000']
mc2016_samples['2016_LFVStTensorU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ3_emutu/'], 'mc','','2016', '','1.900' ,'35.92','462000']
mc2016_samples['2016_LFVTtTensorC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_TT_clequ3_emutc/'], 'mc','','2016', '','0.1876','35.92','290000']
mc2016_samples['2016_LFVTtTensorU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2016/IIHE_Ntuple/ntuple_SMEFTfr_TT_clequ3_emutu/'], 'mc','','2016', '','0.1876','35.92','178000']

#########################################
data2016_samples['2016_B_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runB/191203_070326/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runB/191203_070326/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runB/191203_070326/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runB/191203_070326/0003/'], 'data','MuonEG','2016', 'B','1','1','1']
data2016_samples['2016_C_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0001/'], 'data','MuonEG','2016', 'C','1','1','1']
data2016_samples['2016_D_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runD/191203_070404/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runD/191203_070404/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runD/191203_070404/0002/'], 'data','MuonEG','2016', 'D','1','1','1']
data2016_samples['2016_E_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runE/191203_071509/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runE/191203_071509/0001/'], 'data','MuonEG','2016', 'E','1','1','1']
data2016_samples['2016_F_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runF/191203_070443/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runF/191203_070443/0001/'], 'data','MuonEG','2016', 'F','1','1','1']
data2016_samples['2016_G_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runG/191203_070502/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runG/191203_070502/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runG/191203_070502/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runG/191203_070502/0003/'], 'data','MuonEG','2016', 'G','1','1','1']
data2016_samples['2016_H_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runH/191203_070522/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runH/191203_070522/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runH/191203_070522/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runH/191203_070522/0003/'], 'data','MuonEG','2016', 'H','1','1','1']

data2016_samples['2016_B_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runB/191203_065806/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runB/191203_065806/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runB/191203_065806/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runB/191203_065806/0003/'], 'data','SingleMuon','2016', 'B','1','1','1']
data2016_samples['2016_C_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runC/191203_065920/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runC/191203_065920/0001/'], 'data','SingleMuon','2016', 'C','1','1','1']
data2016_samples['2016_D_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runD/191203_065939/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runD/191203_065939/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runD/191203_065939/0002/'], 'data','SingleMuon','2016', 'D','1','1','1']
data2016_samples['2016_E_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runE/191203_065958/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runE/191203_065958/0001/'], 'data','SingleMuon','2016', 'E','1','1','1']
data2016_samples['2016_F_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runF/191203_070017/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runF/191203_070017/0001/'], 'data','SingleMuon','2016', 'F','1','1','1']
data2016_samples['2016_G_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runG/191203_070036/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runG/191203_070036/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runG/191203_070036/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runG/191203_070036/0003/'], 'data','SingleMuon','2016', 'G','1','1','1']
data2016_samples['2016_H_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runH/191203_070054/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runH/191203_070054/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runH/191203_070054/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runH/191203_070054/0003/'], 'data','SingleMuon','2016', 'H','1','1','1']

data2016_samples['2016_B_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runB/191203_070113/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runB/191203_070113/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runB/191203_070113/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runB/191203_070113/0003/'], 'data','SingleElectron','2016', 'B','1','1','1']
data2016_samples['2016_C_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runC/191203_070132/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runC/191203_070132/0001/'], 'data','SingleElectron','2016', 'C','1','1','1']
data2016_samples['2016_D_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runD/191203_070151/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runD/191203_070151/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runD/191203_070151/0002/'], 'data','SingleElectron','2016', 'D','1','1','1']
data2016_samples['2016_E_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runE/191203_070210/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runE/191203_070210/0001/'], 'data','SingleElectron','2016', 'E','1','1','1']
data2016_samples['2016_F_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runF/191203_070229/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runF/191203_070229/0001/'], 'data','SingleElectron','2016', 'F','1','1','1']
data2016_samples['2016_G_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runG/191203_070247/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runG/191203_070247/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runG/191203_070247/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runG/191203_070247/0003/'], 'data','SingleElectron','2016', 'G','1','1','1']
data2016_samples['2016_H_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runH/191203_070307/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runH/191203_070307/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runH/191203_070307/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runH/191203_070307/0003/'], 'data','SingleElectron','2016', 'H','1','1','1']


data2016_samples['2016_B_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runB/200306_134908/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runB/200306_134908/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runB/200306_134908/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runB/200306_134908/0003/'], 'data','DoubleEG','2016', 'B','1','1','1']
data2016_samples['2016_C_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runC/200306_134951/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runC/200306_134951/0001/'], 'data','DoubleEG','2016', 'C','1','1','1']
data2016_samples['2016_D_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runD/200306_135009/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runD/200306_135009/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runD/200306_135009/0002/'], 'data','DoubleEG','2016', 'D','1','1','1']
data2016_samples['2016_E_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runE/200306_135026/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runE/200306_135026/0001/'], 'data','DoubleEG','2016', 'E','1','1','1']
data2016_samples['2016_F_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runF/200306_135047/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runF/200306_135047/0001/'], 'data','DoubleEG','2016', 'F','1','1','1']
data2016_samples['2016_G_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runG/200306_135106/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runG/200306_135106/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runG/200306_135106/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runG/200306_135106/0003/'], 'data','DoubleEG','2016', 'G','1','1','1']
data2016_samples['2016_H_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runH/200306_135125/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runH/200306_135125/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runH/200306_135125/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runH/200306_135125/0003/'], 'data','DoubleEG','2016', 'H','1','1','1']

data2016_samples['2016_B_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runB/200306_135144/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runB/200306_135144/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runB/200306_135144/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runB/200306_135144/0003/'], 'data','DoubleMu','2016', 'B','1','1','1']
data2016_samples['2016_C_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runC/200306_135205/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runC/200306_135205/0001/'], 'data','DoubleMu','2016', 'C','1','1','1']
data2016_samples['2016_D_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runD/200306_135224/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runD/200306_135224/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runD/200306_135224/0002/'], 'data','DoubleMu','2016', 'D','1','1','1']
data2016_samples['2016_E_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runE/200306_135245/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runE/200306_135245/0001/'], 'data','DoubleMu','2016', 'E','1','1','1']
data2016_samples['2016_F_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runF/200306_135305/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runF/200306_135305/0001/'], 'data','DoubleMu','2016', 'F','1','1','1']
data2016_samples['2016_G_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runG/200306_135327/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runG/200306_135327/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runG/200306_135327/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runG/200306_135327/0003/'], 'data','DoubleMu','2016', 'G','1','1','1']
data2016_samples['2016_H_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runH/200306_135349/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runH/200306_135349/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runH/200306_135349/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runH/200306_135349/0003/'], 'data','DoubleMu','2016', 'H','1','1','1']
