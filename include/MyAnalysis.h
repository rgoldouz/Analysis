//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Feb 22 14:33:31 2020 by ROOT version 6.10/05
// from TTree IIHEAnalysis/IIHEAnalysis
// found on file: /pnfs/iihe/cms/store/user/schenara/SYS_2018/TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8/crab_TTTo2L2Nu_hdampUP_2018/200211_122927/0000/outfile_5_2016_71.root
//////////////////////////////////////////////////////////

#ifndef MyAnalysis_h
#define MyAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

using namespace std;
class MyAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       event;
   ULong64_t       run;
   ULong64_t       luminosityBlock;


   UInt_t          nElectron;
   vector<float>   *Electron_deltaEtaSC;
   vector<float>   *Electron_dr03EcalRecHitSumEt;
   vector<float>   *Electron_dr03HcalDepth1TowerSumEt;
   vector<float>   *Electron_dr03TkSumPt;
   vector<float>   *Electron_dr03TkSumPtHEEP;
   vector<float>   *Electron_dxy;
   vector<float>   *Electron_dxyErr;
   vector<float>   *Electron_dz;
   vector<float>   *Electron_dzErr;
   vector<float>   *Electron_eCorr;
   vector<float>   *Electron_eInvMinusPInv;
   vector<float>   *Electron_energyErr;
   vector<float>   *Electron_eta;
   vector<float>   *Electron_hoe;
   vector<float>   *Electron_ip3d;
   vector<float>   *Electron_jetPtRelv2;
   vector<float>   *Electron_jetRelIso;
   vector<float>   *Electron_mass;
   vector<float>   *Electron_miniPFRelIso_all;
   vector<float>   *Electron_miniPFRelIso_chg;
   vector<float>   *Electron_mvaFall17V1Iso;
   vector<float>   *Electron_mvaFall17V1noIso;
   vector<float>   *Electron_mvaFall17V2Iso;
   vector<float>   *Electron_mvaFall17V2noIso;
   vector<float>   *Electron_pfRelIso03_all;
   vector<float>   *Electron_pfRelIso03_chg;
   vector<float>   *Electron_phi;
   vector<float>   *Electron_pt;
   vector<float>   *Electron_r9;
   vector<float>   *Electron_sieie;
   vector<float>   *Electron_sip3d;
   vector<float>   *Electron_mvaTOP;
   vector<float>   *Electron_mvaTTH;


   vector<int>     *Electron_charge;
   vector<int>     *Electron_cutBased;
   vector<int>     *Electron_cutBased_Fall17_V1;
   vector<int>     *Electron_jetIdx;
   vector<int>     *Electron_jetNDauChargedMVASel;
   vector<int>     *Electron_pdgId;
   vector<int>     *Electron_photonIdx;
   vector<int>     *Electron_tightCharge;
   vector<int>     *Electron_vidNestedWPBitmap;
   vector<int>     *Electron_vidNestedWPBitmapHEEP;
   vector<int>     *Electron_convVeto;
   vector<int>     *Electron_cutBased_HEEP;
   vector<int>     *Electron_isPFcand;
   vector<int>     *Electron_lostHits;
   vector<int>     *Electron_mvaFall17V1Iso_WP80;
   vector<int>     *Electron_mvaFall17V1Iso_WP90;
   vector<int>     *Electron_mvaFall17V1Iso_WPL;
   vector<int>     *Electron_mvaFall17V1noIso_WP80;
   vector<int>     *Electron_mvaFall17V1noIso_WP90;
   vector<int>     *Electron_mvaFall17V1noIso_WPL;
   vector<int>     *Electron_mvaFall17V2Iso_WP80;
   vector<int>     *Electron_mvaFall17V2Iso_WP90;
   vector<int>     *Electron_mvaFall17V2Iso_WPL;
   vector<int>     *Electron_mvaFall17V2noIso_WP80;
   vector<int>     *Electron_mvaFall17V2noIso_WP90;
   vector<int>     *Electron_mvaFall17V2noIso_WPL;
   vector<int>     *Electron_seedGain;

   UInt_t          nMuon;

   vector<float>   *Muon_dxy;
   vector<float>   *Muon_dxyErr;
   vector<float>   *Muon_dz;
   vector<float>   *Muon_dzErr;
   vector<float>   *Muon_eta;
   vector<float>   *Muon_ip3d;
   vector<float>   *Muon_jetPtRelv2;
   vector<float>   *Muon_jetRelIso;
   vector<float>   *Muon_mass;
   vector<float>   *Muon_miniPFRelIso_all;
   vector<float>   *Muon_miniPFRelIso_chg;
   vector<float>   *Muon_pfRelIso03_all;
   vector<float>   *Muon_pfRelIso03_chg;
   vector<float>   *Muon_pfRelIso04_all;
   vector<float>   *Muon_phi;
   vector<float>   *Muon_pt;
   vector<float>   *Muon_ptErr;
   vector<float>   *Muon_segmentComp;
   vector<float>   *Muon_sip3d;
   vector<float>   *Muon_tkRelIso;
   vector<float>   *Muon_tunepRelPt;
   vector<float>   *Muon_mvaLowPt;
   vector<float>   *Muon_mvaTOP;
   vector<float>   *Muon_mvaTTH;
   vector<int>   *Muon_charge;
   vector<int>   *Muon_jetIdx;
   vector<int>   *Muon_jetNDauChargedMVASel;
   vector<int>   *Muon_nStations;
   vector<int>   *Muon_nTrackerLayers;
   vector<int>   *Muon_pdgId;
   vector<int>   *Muon_tightCharge;
   vector<int>   *Muon_fsrPhotonIdx;
   vector<int>   *Muon_highPtId;
   vector<int>   *Muon_inTimeMuon;
   vector<int>   *Muon_isGlobal;
   vector<int>   *Muon_isPFcand;
   vector<int>   *Muon_isTracker;
   vector<int>   *Muon_looseId;
   vector<int>   *Muon_mediumId;
   vector<int>   *Muon_mediumPromptId;
   vector<int>   *Muon_miniIsoId;
   vector<int>   *Muon_multiIsoId;
   vector<int>   *Muon_mvaId;
   vector<int>   *Muon_pfIsoId;
   vector<int>   *Muon_softId;
   vector<int>   *Muon_softMvaId;
   vector<int>   *Muon_tightId;
   vector<int>   *Muon_tkIsoId;
   vector<int>   *Muon_triggerIdLoose;

   vector<float>   *MET_phi;
   vector<float>   *MET_pt;

         UInt_t    nJet; 
   vector<float>   *Jet_area;
   vector<float>   *Jet_btagCMVA;
   vector<float>   *Jet_btagCSVV2; 
   vector<int>     *Jet_btagDeepB;
   vector<int>     *Jet_btagDeepC;
   vector<float>   *Jet_btagDeepFlavB;
   vector<float>   *Jet_btagDeepFlavC;
   vector<float>   *Jet_chEmEF;
   vector<float>   *Jet_chHEF;
   vector<float>   *Jet_eta;
   vector<float>   *Jet_jercCHF; 
   vector<float>   *Jet_jercCHPUF;
   vector<float>   *Jet_mass;
   vector<float>   *Jet_muEF;
   vector<float>   *Jet_muonSubtrFactor;
   vector<float>   *Jet_neEmEF;
   vector<float>   *Jet_neHEF;
   vector<float>   *Jet_phi;
   vector<float>   *Jet_pt;
   vector<float>   *Jet_qgl;
   vector<float>   *Jet_rawFactor;

   vector<float>   *Jet_bRegCorr;
   vector<float>   *Jet_bRegRes;
   vector<int>     *Jet_electronIdx1;
   vector<int>     *Jet_electronIdx2;
   vector<int>     *Jet_jetId;
   vector<int>     *Jet_muonIdx1;
   vector<int>     *Jet_muonIdx2;
   vector<int>     *Jet_nConstituents;
   vector<int>     *Jet_nElectrons;
   vector<int>     *Jet_nMuons;
   vector<int>     *Jet_puId;
   vector<int>     *Jet_partonFlavour;



         Int_t     HLT_Ele35_WPTight_Gsf;

         Int_t     HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;

         Int_t     HLT_IsoMu27;

         Int_t     HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;


         Int_t     HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ; 
         Int_t     HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;

         Int_t     HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;

         Int_t     PV_npvsGood;

         float     L1PreFiringWeight_Dn;
         float     L1PreFiringWeight_Nom;
         float     L1PreFiringWeight_Up;

   vector<float>   *GenPart_mass;
   vector<float>   *GenPart_phi;
   vector<float>   *GenPart_pt;
   vector<float>   *GenPart_eta;
   vector<int>     *GenPart_genPartIdxMother;
   vector<int>     *GenPart_pdgId;
   vector<int>     *GenPart_status;


   // List of branches
  TBranch        *b_event;
   TBranch        *b_run;
   TBranch        *b_luminosityBlock;

   TBranch        *b_nElectron;
   TBranch        *b_Electron_deltaEtaSC;
   TBranch        *b_Electron_dr03EcalRecHitSumEt;
   TBranch        *b_Electron_dr03HcalDepth1TowerSumEt;
   TBranch        *b_Electron_dr03TkSumPt;
   TBranch        *b_Electron_dr03TkSumPtHEEP;
   TBranch        *b_Electron_dxy;
   TBranch        *b_Electron_dxyErr;
   TBranch        *b_Electron_dz;
   TBranch        *b_Electron_dzErr;
   TBranch        *b_Electron_eCorr;
   TBranch        *b_Electron_eInvMinusPInv;
   TBranch        *b_Electron_energyErr;
   TBranch        *b_Electron_eta;
   TBranch        *b_Electron_hoe;
   TBranch        *b_Electron_ip3d;
   TBranch        *b_Electron_jetPtRelv2;
   TBranch        *b_Electron_jetRelIso;
   TBranch        *b_Electron_mass;
   TBranch        *b_Electron_miniPFRelIso_all;
   TBranch        *b_Electron_miniPFRelIso_chg;
   TBranch        *b_Electron_mvaFall17V1Iso;
   TBranch        *b_Electron_mvaFall17V1noIso;
   TBranch        *b_Electron_mvaFall17V2Iso;
   TBranch        *b_Electron_mvaFall17V2noIso;
   TBranch        *b_Electron_pfRelIso03_all;
   TBranch        *b_Electron_pfRelIso03_chg;
   TBranch        *b_Electron_phi;
   TBranch        *b_Electron_pt;
   TBranch        *b_Electron_r9;
   TBranch        *b_Electron_sieie;
   TBranch        *b_Electron_sip3d;
   TBranch        *b_Electron_mvaTOP;
   TBranch        *b_Electron_mvaTTH;


   TBranch        *b_Electron_charge;
   TBranch        *b_Electron_cutBased;
   TBranch        *b_Electron_cutBased_Fall17_V1;
   TBranch        *b_Electron_jetIdx;
   TBranch        *b_Electron_jetNDauChargedMVASel;
   TBranch        *b_Electron_pdgId;
   TBranch        *b_Electron_photonIdx;
   TBranch        *b_Electron_tightCharge;
   TBranch        *b_Electron_vidNestedWPBitmap;
   TBranch        *b_Electron_vidNestedWPBitmapHEEP;
   TBranch        *b_Electron_convVeto;
   TBranch        *b_Electron_cutBased_HEEP;
   TBranch        *b_Electron_isPFcand;
   TBranch        *b_Electron_lostHits;
   TBranch        *b_Electron_mvaFall17V1Iso_WP80;
   TBranch        *b_Electron_mvaFall17V1Iso_WP90;
   TBranch        *b_Electron_mvaFall17V1Iso_WPL;
   TBranch        *b_Electron_mvaFall17V1noIso_WP80;
   TBranch        *b_Electron_mvaFall17V1noIso_WP90;
   TBranch        *b_Electron_mvaFall17V1noIso_WPL;
   TBranch        *b_Electron_mvaFall17V2Iso_WP80;
   TBranch        *b_Electron_mvaFall17V2Iso_WP90;
   TBranch        *b_Electron_mvaFall17V2Iso_WPL;
   TBranch        *b_Electron_mvaFall17V2noIso_WP80;
   TBranch        *b_Electron_mvaFall17V2noIso_WP90;
   TBranch        *b_Electron_mvaFall17V2noIso_WPL;
   TBranch        *b_Electron_seedGain;

   TBranch        *b_nMuon;

   TBranch        *b_Muon_dxy;
   TBranch        *b_Muon_dxyErr;
   TBranch        *b_Muon_dz;
   TBranch        *b_Muon_dzErr;
   TBranch        *b_Muon_eta;
   TBranch        *b_Muon_ip3d;
   TBranch        *b_Muon_jetPtRelv2;
   TBranch        *b_Muon_jetRelIso;
   TBranch        *b_Muon_mass;
   TBranch        *b_Muon_miniPFRelIso_all;
   TBranch        *b_Muon_miniPFRelIso_chg;
   TBranch        *b_Muon_pfRelIso03_all;
   TBranch        *b_Muon_pfRelIso03_chg;
   TBranch        *b_Muon_pfRelIso04_all;
   TBranch        *b_Muon_phi;
   TBranch        *b_Muon_pt;
   TBranch        *b_Muon_ptErr;
   TBranch        *b_Muon_segmentComp;
   TBranch        *b_Muon_sip3d;
   TBranch        *b_Muon_tkRelIso;
   TBranch        *b_Muon_tunepRelPt;
   TBranch        *b_Muon_mvaLowPt;
   TBranch        *b_Muon_mvaTOP;
   TBranch        *b_Muon_mvaTTH;
   TBranch        *b_Muon_charge;
   TBranch        *b_Muon_jetIdx;
   TBranch        *b_Muon_jetNDauChargedMVASel;
   TBranch        *b_Muon_nStations;
   TBranch        *b_Muon_nTrackerLayers;
   TBranch        *b_Muon_pdgId;
   TBranch        *b_Muon_tightCharge;
   TBranch        *b_Muon_fsrPhotonIdx;
   TBranch        *b_Muon_highPtId;
   TBranch        *b_Muon_inTimeMuon;
   TBranch        *b_Muon_isGlobal;
   TBranch        *b_Muon_isPFcand;
   TBranch        *b_Muon_isTracker;
   TBranch        *b_Muon_looseId;
   TBranch        *b_Muon_mediumId;
   TBranch        *b_Muon_mediumPromptId;
   TBranch        *b_Muon_miniIsoId;
   TBranch        *b_Muon_multiIsoId;
   TBranch        *b_Muon_mvaId;
   TBranch        *b_Muon_pfIsoId;
   TBranch        *b_Muon_softId;
   TBranch        *b_Muon_softMvaId;
   TBranch        *b_Muon_tightId;
   TBranch        *b_Muon_tkIsoId;
   TBranch        *b_Muon_triggerIdLoose;

   TBranch        *b_MET_phi;
   TBranch        *b_MET_pt;


   TBranch        *b_nJet; 
   TBranch        *b_Jet_partonFlavour; 
   TBranch        *b_Jet_area;
   TBranch        *b_Jet_btagCMVA;
   TBranch        *b_Jet_btagCSVV2; 
   TBranch        *b_Jet_btagDeepB;
   TBranch        *b_Jet_btagDeepC;
   TBranch        *b_Jet_btagDeepFlavB;
   TBranch        *b_Jet_btagDeepFlavC;
   TBranch        *b_Jet_chEmEF;
   TBranch        *b_Jet_chHEF;
   TBranch        *b_Jet_eta;
   TBranch        *b_Jet_jercCHF; 
   TBranch        *b_Jet_jercCHPUF;
   TBranch        *b_Jet_mass;
   TBranch        *b_Jet_muEF;
   TBranch        *b_Jet_muonSubtrFactor;
   TBranch        *b_Jet_neEmEF;
   TBranch        *b_Jet_neHEF;
   TBranch        *b_Jet_phi;
   TBranch        *b_Jet_pt;
   TBranch        *b_Jet_qgl;
   TBranch        *b_Jet_rawFactor;
   TBranch        *b_Jet_bRegCorr;
   TBranch        *b_Jet_bRegRes;
   TBranch        *b_Jet_electronIdx1;
   TBranch        *b_Jet_electronIdx2;
   TBranch        *b_Jet_jetId;
   TBranch        *b_Jet_muonIdx1;
   TBranch        *b_Jet_muonIdx2;
   TBranch        *b_Jet_nConstituents;
   TBranch        *b_Jet_nElectrons;
   TBranch        *b_Jet_nMuons;
   TBranch        *b_Jet_puId;
         



   TBranch        *b_HLT_Ele35_WPTight_Gsf;

   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;

   TBranch        *b_HLT_IsoMu27;

   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;


   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ; 
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
    

   TBranch        *b_PV_npvsGood;

   TBranch        *b_L1PreFiringWeight_Dn;
   TBranch        *b_L1PreFiringWeight_Nom;
   TBranch        *b_L1PreFiringWeight_Up;

   TBranch        *b_GenPart_mass;
   TBranch        *b_GenPart_phi;
   TBranch        *b_GenPart_pt;
   TBranch        *b_GenPart_eta;
   TBranch        *b_GenPart_genPartIdxMother;
   TBranch        *b_GenPart_pdgId;
   TBranch        *b_GenPart_status;



   MyAnalysis(TTree *tree=0);
   virtual ~MyAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString, TString, TString, TString, TString, float,float,float);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyAnalysis_cxx
MyAnalysis::MyAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/pnfs/iihe/cms/store/user/schenara/SYS_2018/TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8/crab_TTTo2L2Nu_hdampUP_2018/200211_122927/0000/outfile_5_2016_71.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/pnfs/iihe/cms/store/user/schenara/SYS_2018/TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8/crab_TTTo2L2Nu_hdampUP_2018/200211_122927/0000/outfile_5_2016_71.root");
      }
      f->GetObject("IIHEAnalysis",tree);

   }
   Init(tree);
}

MyAnalysis::~MyAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);

   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_deltaEtaSC", &Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", &Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", &Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", &Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_dr03TkSumPt", &Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_dr03TkSumPtHEEP", &Electron_dr03TkSumPtHEEP, &b_Electron_dr03TkSumPtHEEP);
   fChain->SetBranchAddress("Electron_dxy", &Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", &Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_eta", &Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_mass", &Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_dz", &Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_eCorr", &Electron_eCorr, &b_Electron_eCorr);
   fChain->SetBranchAddress("Electron_phi", &Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", &Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_miniPFRelIso_all", &Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("Electron_miniPFRelIso_chg", &Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso", &Electron_mvaFall17V2Iso, &b_Electron_mvaFall17V2Iso);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso", &Electron_mvaFall17V2noIso, &b_Electron_mvaFall17V2noIso);
   fChain->SetBranchAddress("Electron_pfRelIso03_all", &Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_pfRelIso03_chg", &Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("Electron_mvaTOP", &Electron_mvaTOP, &b_Electron_mvaTOP);
   fChain->SetBranchAddress("Electron_mvaTTH", &Electron_mvaTTH, &b_Electron_mvaTTH);
   fChain->SetBranchAddress("Electron_sip3d", &Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_pdgId", &Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_charge", &Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_cutBased", &Electron_cutBased, &b_Electron_cutBased);


   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_dxy", &Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr", &Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dz", &Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr", &Muon_dzErr, &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d", &Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_mass", &Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_miniPFRelIso_all", &Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all)

   ;
   fChain->SetBranchAddress("Muon_pfRelIso03_all", &Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_mvaTOP", &Muon_mvaTOP, &b_Muon_mvaTOP);
   fChain->SetBranchAddress("Muon_mvaTTH", &Muon_mvaTTH, &b_Muon_mvaTTH);
   fChain->SetBranchAddress("Muon_mediumId", &Muon_mediumId, &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_mvaId", &Muon_mvaId, &b_Muon_mvaId);
   fChain->SetBranchAddress("Muon_miniIsoId", &Muon_miniIsoId, &b_Muon_miniIsoId);


   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);


   fChain->SetBranchAddress("nJet", &nJet, &b_nJet); 
   fChain->SetBranchAddress("Jet_partonFlavour", &Jet_partonFlavour, &b_Jet_partonFlavour); 

   fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_jetId", &Jet_jetId, &b_Jet_jetId);
   fChain->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB, &b_Jet_btagDeepB);


   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf, &b_HLT_Ele35_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);

   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);

   fChain->SetBranchAddress("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom, &b_L1PreFiringWeight_Nom);
   fChain->SetBranchAddress("L1PreFiringWeight_Dn", &L1PreFiringWeight_Dn, &b_L1PreFiringWeight_Dn);
   fChain->SetBranchAddress("L1PreFiringWeight_Up", &L1PreFiringWeight_Up, &b_L1PreFiringWeight_Up);


   fChain->SetBranchAddress("GenPart_mass", &GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi", &GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt", &GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_eta", &GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status", &GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);



   Notify();


}

Bool_t MyAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyAnalysis_cxx
