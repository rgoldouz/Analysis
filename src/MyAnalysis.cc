#define MyAnalysis_cxx
#include "MyAnalysis.h"
#include "PU_reWeighting.h"
#include "lepton_candidate.h"
#include "jet_candidate.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <vector>

void displayProgress(long current, long max){
  using std::cerr;
  if (max<2500) return;
  if (current%(max/2500)!=0 && current<max-1) return;

  int width = 52; // Hope the terminal is at least that wide.
  int barWidth = width - 2;
  cerr << "\x1B[2K"; // Clear line
  cerr << "\x1B[2000D"; // Cursor left
  cerr << '[';
  for(int i=0 ; i<barWidth ; ++i){ if(i<barWidth*current/max){ cerr << '=' ; }else{ cerr << ' ' ; } }
  cerr << ']';
  cerr << " " << Form("%8d/%8d (%5.2f%%)", (int)current, (int)max, 100.0*current/max) ;
  cerr.flush();
}

Double_t deltaPhi(Double_t phi1, Double_t phi2) {
  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
  return dPhi;
}


Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);
  return sqrt(dEta*dEta+dPhi*dPhi);
}

bool ComparePt(lepton_candidate *a, lepton_candidate *b) { return a->pt_ > b->pt_; }

float scale_factor( TH2F* h, float X, float Y , TString uncert){
  int NbinsX=h->GetXaxis()->GetNbins();
  int NbinsY=h->GetYaxis()->GetNbins();
  float x_min=h->GetXaxis()->GetBinLowEdge(1);
  float x_max=h->GetXaxis()->GetBinLowEdge(NbinsX)+h->GetXaxis()->GetBinWidth(NbinsX);
  float y_min=h->GetYaxis()->GetBinLowEdge(1);
  float y_max=h->GetYaxis()->GetBinLowEdge(NbinsY)+h->GetYaxis()->GetBinWidth(NbinsY);
  TAxis *Xaxis = h->GetXaxis();
  TAxis *Yaxis = h->GetYaxis();
  Int_t binx=1;
  Int_t biny=1;
  if(x_min < X && X < x_max) binx = Xaxis->FindBin(X);
  else binx= (X<=x_min) ? 1 : NbinsX ;
  if(y_min < Y && Y < y_max) biny = Yaxis->FindBin(Y);
  else biny= (Y<=y_min) ? 1 : NbinsY ;
  if(uncert=="up") return (h->GetBinContent(binx, biny)+h->GetBinError(binx, biny));
  else if(uncert=="down") return (h->GetBinContent(binx, biny)-h->GetBinError(binx, biny));
  else return  h->GetBinContent(binx, biny);
}


//float MuID_scale_factor( TH2F* h, float pt, float eta , TString uncert){
//  int ptbin  = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(pt)));
//  int etabin = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindBin(fabs(eta))));
//cout<<pt<<" "<<eta<<endl;
//cout<<ptbin<<" "<<etabin<<endl;
//  if (uncert=="up") return (h->GetBinContent(ptbin,etabin)+h->GetBinError(ptbin,etabin));
//  else if(uncert=="down" )return (h->GetBinContent(ptbin,etabin)-h->GetBinError(ptbin,etabin));
//  else return h->GetBinContent(ptbin,etabin);
//}
//
void MyAnalysis::Loop(TString fname, TString data, TString year, TString run, float xs, float lumi, float Nevent)
{

//Get scale factor and weight histograms
  TH2F  sf_Ele_Reco_H;
  TH2F  sf_Ele_ID_H;
  TH2F  sf_Mu_ID_H;
  TH2F  sf_Mu_ISO_H;
  TH2F  sf_trigger_H;
  PU wPU;

  if(data == "mc"){
    if(year == "2016"){
      TFile *f_Ele_Reco_Map = new TFile("input/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("input/2016LegacyReReco_ElectronTight_Fall17V2.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map_1 = new TFile("input/2016_RunBCDEF_SF_ID.root");
      TH2F *sf_Mu_ID_H_1 = (TH2F*)f_Mu_ID_Map_1->Get("NUM_TightID_DEN_genTracks_eta_pt");
      TFile *f_Mu_ID_Map_2 = new TFile("input/2016_RunGH_SF_ID.root");
      TH2F *sf_Mu_ID_H_2 = (TH2F*)f_Mu_ID_Map_2->Get("NUM_TightID_DEN_genTracks_eta_pt");
      sf_Mu_ID_H_1->Scale(0.55);
      sf_Mu_ID_H_2->Scale(0.45);
      sf_Mu_ID_H_1->Add(sf_Mu_ID_H_2);
      sf_Mu_ID_H = *sf_Mu_ID_H_1;

      TFile *f_Mu_ISO_Map_1 = new TFile("input/2016_RunBCDEF_SF_ISO.root");
      TH2F *sf_Mu_ISO_H_1 = (TH2F*)f_Mu_ISO_Map_1->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");
      TFile *f_Mu_ISO_Map_2 = new TFile("input/2016_RunGH_SF_ISO.root");
      TH2F *sf_Mu_ISO_H_2 = (TH2F*)f_Mu_ISO_Map_2->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");
      sf_Mu_ISO_H_1->Scale(0.55);
      sf_Mu_ISO_H_2->Scale(0.45);
      sf_Mu_ISO_H_1->Add(sf_Mu_ISO_H_2);
      sf_Mu_ISO_H = *sf_Mu_ISO_H_1;

      TFile *f_trigger = new TFile("input/TriggerSF_emu2016_pt.root");
      sf_trigger_H = *(TH2F*)f_trigger->Get("h_lep1Pt_lep2Pt_Step3");

      f_Ele_Reco_Map->Close();
      f_Ele_ID_Map->Close();
      f_Mu_ID_Map_1->Close();
      f_Mu_ID_Map_2->Close();
      f_Mu_ISO_Map_1->Close();
      f_Mu_ISO_Map_2->Close();
      f_trigger->Close();
    }
    if(year == "2017"){
      TFile *f_Ele_Reco_Map = new TFile("input/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("input/2017_ElectronTight.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map = new TFile("input/2017_RunBCDEF_SF_ID_syst.root");
      sf_Mu_ID_H = *(TH2F*)f_Mu_ID_Map->Get("NUM_TightID_DEN_genTracks_pt_abseta");

      TFile *f_Mu_ISO_Map = new TFile("input/2017_RunBCDEF_SF_ISO_syst.root");
      sf_Mu_ISO_H = *(TH2F*)f_Mu_ISO_Map->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

      TFile *f_trigger = new TFile("input/TriggerSF_emu2017_pt.root");
      sf_trigger_H = *(TH2F*)f_trigger->Get("h_lep1Pt_lep2Pt_Step3");

      f_Ele_Reco_Map->Close();
      f_Ele_ID_Map->Close();
      f_Mu_ID_Map->Close();
      f_Mu_ISO_Map->Close();
      f_trigger->Close();
    }
    if(year == "2018"){
      TFile *f_Ele_Reco_Map = new TFile("input/egammaEffi.txt_EGM2D_updatedAll.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("input/2018_ElectronTight.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map = new TFile("input/2018_RunABCD_SF_ID.root");
      sf_Mu_ID_H = *(TH2F*)f_Mu_ID_Map->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta");

      TFile *f_Mu_ISO_Map = new TFile("input/2018_RunABCD_SF_ISO.root");
      sf_Mu_ISO_H = *(TH2F*)f_Mu_ISO_Map->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

      TFile *f_trigger = new TFile("input/TriggerSF_emu2018_pt.root");
      sf_trigger_H = *(TH2F*)f_trigger->Get("h_lep1Pt_lep2Pt_Step3");

      f_Ele_Reco_Map->Close();
      f_Ele_ID_Map->Close();
      f_Mu_ID_Map->Close();
      f_Mu_ISO_Map->Close();
      f_trigger->Close();
    }
  }

  TFile file_out (fname,"RECREATE");
  TTree tree_out("analysis","main analysis") ;

  float leadingLepton_pt ;
  float leadingLepton_eta ;
  float leadingLepton_phi ;
  float leadingLepton_E ;
  float leadingLepton_charge ;
  float leadingLepton_type ;

  float subLeadingLepton_pt ;
  float subLeadingLepton_eta ;
  float subLeadingLepton_phi ;
  float subLeadingLepton_E ;
  float subLeadingLepton_type ;

  float met_pt;
  float met_phi;

  vector<float> jets_pt;
  vector<float> jets_eta;
  vector<float> jets_phi;
  vector<float> jets_E;
  vector<int>   jets_btag;

  float weight_lumi ;
  float weight_pu ;
  float sf_eleReco ;
  float sf_eleID ;
  float sf_muID ;
  float sf_muISO ;
  float sf_trigger;

  tree_out.Branch("leadingLepton_pt"    , &leadingLepton_pt    , "leadingLepton_pt/F"    ) ;
  tree_out.Branch("leadingLepton_eta"    , &leadingLepton_eta    , "leadingLepton_eta/F"    ) ;
  tree_out.Branch("leadingLepton_phi"    , &leadingLepton_phi    , "leadingLepton_phi/F"    ) ;
  tree_out.Branch("leadingLepton_E"    , &leadingLepton_E    , "leadingLepton_E/F"    ) ;
  tree_out.Branch("leadingLepton_charge"    , &leadingLepton_charge    , "leadingLepton_charge/F"    ) ;
  tree_out.Branch("leadingLepton_type"    , &leadingLepton_type    , "leadingLepton_type/F"    ) ;
  tree_out.Branch("subLeadingLepton_pt"    , &subLeadingLepton_pt    , "subLeadingLepton_pt/F"    ) ;
  tree_out.Branch("subLeadingLepton_eta"    , &subLeadingLepton_eta    , "subLeadingLepton_eta/F"    ) ;
  tree_out.Branch("subLeadingLepton_phi"    , &subLeadingLepton_phi    , "subLeadingLepton_phi/F"    ) ;
  tree_out.Branch("subLeadingLepton_E"    , &subLeadingLepton_E    , "subLeadingLepton_E/F"    ) ;
  tree_out.Branch("subLeadingLepton_type"    , &subLeadingLepton_type    , "subLeadingLepton_type/F"    ) ;
  tree_out.Branch("jets_pt"    , &jets_pt    ) ;
  tree_out.Branch("jets_eta"   , &jets_eta   ) ;
  tree_out.Branch("jets_phi"   , &jets_phi   ) ;
  tree_out.Branch("jets_E"     , &jets_E   ) ;
  tree_out.Branch("jets_btag"     , &jets_btag   ) ;
  tree_out.Branch("met_pt"    , &met_pt    , "met_pt/F"    ) ;
  tree_out.Branch("met_phi"    , &met_phi    , "met_phi/F"    ) ;
  tree_out.Branch("weight_lumi"    , &weight_lumi    , "weight_lumi/F"    ) ;
  tree_out.Branch("weight_pu"    , &weight_pu    , "weight_pu/F"    ) ;
  tree_out.Branch("sf_eleReco"    , &sf_eleReco    , "sf_eleReco/F"    ) ;
  tree_out.Branch("sf_eleID"    , &sf_eleID    , "sf_eleID/F"    ) ;
  tree_out.Branch("sf_muID"    , &sf_muID    , "sf_muID/F"    ) ;
  tree_out.Branch("sf_muISO"    , &sf_muISO    , "sf_muISO/F"    ) ;
  tree_out.Branch("sf_trigger"    , &sf_trigger    , "sf_trigger/F"    ) ;
//   tree_out.Branch(""    , &    , "/F"    ) ;


    cout<<"ev_event"<<"   "<<"sf_Ele_Reco"<<"   "<<"sf_Ele_ID"<<"      "<<"sf_Mu_ID"<<"   "<<"sf_Mu_ISO"<<"   "<<"sf_trigger"<<"   "<<"PU weight"<<endl;
  std::vector<lepton_candidate*> *selectedLeptons;
  std::vector<jet_candidate*> *selectedJets;
  lepton_candidate *leptest;
  jet_candidate *jettest;

  int nAccept=0;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  Long64_t ntr = fChain->GetEntries ();
//  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  for (Long64_t jentry=0; jentry<100;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    displayProgress(jentry, ntr) ;

    bool triggerPass = false;
    bool metFilterPass = false;
    float sf_Ele_Reco =1;
    float sf_Ele_ID =1;
    float sf_Mu_ID =1;
    float sf_Mu_ISO =1;
    float sf_Trigger =1;
    float weight_PU =1;
    float weight_Lumi =1;
//MET filters

    if(data == "mc"){
      if(year == "2016" || year == "2018" ){
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1)
        metFilterPass = true;
        }
        if(year == "2017"){
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1 && trig_Flag_ecalBadCalibReduced ==1)
        metFilterPass = true;
        }
    }

    if(data == "data"){
      if(year == "2016" || year == "2018"){
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1 &&  trig_Flag_eeBadScFilter_accept==1)
        metFilterPass = true;
        }
        if(year == "2017"){
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1 && trig_Flag_ecalBadCalibReduced ==1)
        metFilterPass = true;
        }
    }

//trigger

      if(data == "mc" && year == "2016"){
        if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele27_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept) {triggerPass =true;}
      }

      if(data == "mc" && year == "2017"){
        if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele35_WPTight_Gsf_accept || trig_HLT_IsoMu27_accept) {triggerPass =true;}
      }

      if(data == "mc" && year == "2018"){
        if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept) {triggerPass =true;}
      } 


    if(data == "data" && year == "2016" && run == "H"){
      if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Ele27_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept) {triggerPass =true;}
    }
    if(data == "data" && year == "2016" && run != "H"){
      if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele27_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept) {triggerPass =true;}    
    }
    if(data == "data" && year == "2017"){
      if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele35_WPTight_Gsf_accept || trig_HLT_IsoMu27_accept) {triggerPass =true;}
    }
    if(data == "data" && year == "2018"){
      if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept) {triggerPass =true;}
    }
    

int check_ev = 50712029; 
//if (ev_event==check_ev) {
//cout<<data<<year<<run<<endl;
//cout<<"triggerPass "<<triggerPass<<endl;
//cout<<"trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept "<<trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept<<endl;
//cout<<"trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept "<<trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept<<endl;
//cout<<"trig_HLT_Ele27_WPTight_Gsf_accept "<<trig_HLT_Ele27_WPTight_Gsf_accept<<endl;
//cout<<"trig_HLT_IsoMu24_accept "<<trig_HLT_IsoMu24_accept<<endl;
//cout<<"trig_HLT_IsoTkMu24_accept "<<trig_HLT_IsoTkMu24_accept<<endl;
//cout<<"metFilterPass "<<metFilterPass<<endl;
//cout<<"trig_Flag_goodVertices_accept "<<trig_Flag_goodVertices_accept<<endl;
//cout<<"trig_Flag_globalSuperTightHalo2016Filter_accept "<<trig_Flag_globalSuperTightHalo2016Filter_accept<<endl;
//cout<<"trig_Flag_HBHENoiseFilter_accept "<<trig_Flag_HBHENoiseFilter_accept<<endl;
//cout<<"trig_Flag_HBHENoiseIsoFilter_accept "<<trig_Flag_HBHENoiseIsoFilter_accept<<endl;
//cout<<"trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept "<<trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept<<endl;
//cout<<"trig_Flag_BadPFMuonFilter_accept "<<trig_Flag_BadPFMuonFilter_accept<<endl;
//cout<<"trig_Flag_eeBadScFilter_accept "<<trig_Flag_eeBadScFilter_accept<<endl;
//
//}
    if(!triggerPass) continue;
    if(!metFilterPass) continue;

// lepton selection
  selectedLeptons = new std::vector<lepton_candidate*>();
// electron
    for (int l=0;l<gsf_pt->size();l++){
      if((*gsf_pt)[l] <20 || abs((*gsf_eta)[l]) > 2.4 || (abs((*gsf_sc_eta)[l])> 1.4442 && (abs((*gsf_sc_eta)[l])< 1.566))) continue;
      if(!(*gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight)[l]) continue;
//      leptest = new lepton_candidate((*gsf_pt)[l],(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1);
      selectedLeptons->push_back(new lepton_candidate((*gsf_pt)[l],(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
      if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
      if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
    }
// Muon
    for (int l=0;l<mu_gt_pt->size();l++){
      if((*mu_gt_pt)[l] <20 || abs((*mu_gt_eta)[l]) > 2.4) continue;
      if(!(*mu_isTightMuon)[l]) continue;
      if((*mu_pfIsoDbCorrected04)[l] > 0.15) continue;
//      leptest = new lepton_candidate((*gsf_pt)[l],(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,10);
      selectedLeptons->push_back(new lepton_candidate((*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l],(*mu_gt_charge)[l],l,10));
      if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
      if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
      if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
      if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");

    }
    sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePt);
// dilepton selection

//if (ev_event==check_ev) {cout<<"lep size "<<(selectedLeptons->size()<2)<<" leading lep pt "<<((*selectedLeptons)[0].pt_ <25)<<" charge "<<((*selectedLeptons)[0].charge_ * (*selectedLeptons)[1].charge_ == 1)<<" emu "<< ((*selectedLeptons)[0].lep_ + (*selectedLeptons)[1].lep_ != 11) <<" mass "<<(((*selectedLeptons)[0].p4_ + (*selectedLeptons)[1].p4_).M()<20)<<endl;}
    if(selectedLeptons->size()<2 ||
      ((*selectedLeptons)[0]->pt_ <25) ||
      ((*selectedLeptons)[0]->charge_ * (*selectedLeptons)[1]->charge_ == 1) ||
      ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ != 11) ||
      ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()<20) {      
      for (int l=0;l<selectedLeptons->size();l++){
        delete (*selectedLeptons)[l];  
      }
      selectedLeptons->clear();
      selectedLeptons->shrink_to_fit();
      delete selectedLeptons;
      continue;
    }
nAccept++;
//jets
    selectedJets = new std::vector<jet_candidate*>();
    bool jetlepfail;
    for (int l=0;l<jet_pt->size();l++){
      if(data == "mc" && ((*jet_Smeared_pt)[l] <30 || abs((*jet_eta)[l]) > 2.4)) continue;
      if(data == "data" && ((*jet_pt)[l] <30 || abs((*jet_eta)[l])) > 2.4) continue;
      if(year == "2016" && !(*jet_isJetIDTightLepVeto_2016)[l]) continue;
      if(year == "2017" && !(*jet_isJetIDLepVeto_2017)[l]) continue;
      if(year == "2018" && !(*jet_isJetIDLepVeto_2018)[l]) continue;
      jetlepfail = false;
      for (int i=0;i<selectedLeptons->size();i++){
        if(deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,(*jet_eta)[l],(*jet_phi)[l]) < 0.4 ) jetlepfail=true;
      }
      if(jetlepfail) continue; 
      if(data == "mc"){
        jettest = new jet_candidate((*jet_Smeared_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l], year,l);
        selectedJets->push_back(new jet_candidate((*jet_Smeared_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l], year,l));
      }
      if(data == "data"){
        jettest = new jet_candidate((*jet_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l],year,l);
        selectedJets->push_back(new jet_candidate((*jet_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l],year,l));
      }
    }

    int nbjet=0;
    for (int l=0;l<selectedJets->size();l++){
      jets_pt.push_back((*selectedJets)[l]->pt_);
      jets_eta.push_back((*selectedJets)[l]->eta_);
      jets_phi.push_back((*selectedJets)[l]->phi_);
      jets_E.push_back((*selectedJets)[l]->p4_.E());
      jets_btag.push_back((*selectedJets)[l]->btag_);
      if((*selectedJets)[l]->btag_) nbjet++;
    }

    if (data == "mc") sf_Trigger = scale_factor(&sf_trigger_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"");
    if (data == "mc" && year == "2016") weight_PU = wPU.PU_2016(mc_trueNumInteractions,"nominal");
    if (data == "mc" && year == "2017") weight_PU = wPU.PU_2017(mc_trueNumInteractions,"nominal");
    if (data == "mc" && year == "2018") weight_PU = wPU.PU_2018(mc_trueNumInteractions,"nominal");
    if (data == "mc") weight_Lumi = (xs*lumi)/Nevent;

    leadingLepton_pt = (*selectedLeptons)[0]->pt_;
    leadingLepton_eta = (*selectedLeptons)[0]->eta_ ;
    leadingLepton_phi = (*selectedLeptons)[0]->phi_;
    leadingLepton_E = (*selectedLeptons)[0]->p4_.E();
    leadingLepton_charge = (*selectedLeptons)[0]->charge_;
    leadingLepton_type = (*selectedLeptons)[0]->lep_;

    subLeadingLepton_pt = (*selectedLeptons)[1]->pt_;
    subLeadingLepton_eta = (*selectedLeptons)[1]->eta_ ;
    subLeadingLepton_phi = (*selectedLeptons)[1]->phi_;
    subLeadingLepton_E = (*selectedLeptons)[1]->p4_.E();
    subLeadingLepton_type = (*selectedLeptons)[1]->lep_;

    met_pt = MET_FinalCollection_Pt;
    met_phi = MET_FinalCollection_phi;
    weight_lumi = weight_Lumi;
    weight_pu = weight_PU;
    sf_eleReco = sf_Ele_Reco ;
    sf_eleID = sf_Ele_ID ;
    sf_muID = sf_Mu_ID;
    sf_muISO = sf_Mu_ISO ;
    sf_trigger = sf_Trigger ;
    tree_out.Fill() ;

     cout<<ev_event<<"   "<<sf_Ele_Reco<<"   "<<sf_Ele_ID<<"      "<<sf_Mu_ID<<"   "<<sf_Mu_ISO<<"   "<<sf_Trigger<<"   "<<weight_PU<<endl;
    if(selectedJets->size()<3 || MET_FinalCollection_Pt>30 || nbjet !=1) continue;
//    nAccept++;

//if (ev_event==check_ev) 
//cout<<ev_event<<"   "<<selectedJets->size()<<"   "<<MET_FinalCollection_Pt<<"     "<<nbjet<<endl;

    for (int l=0;l<selectedLeptons->size();l++){
      delete (*selectedLeptons)[l];
    }
    for (int l=0;l<selectedJets->size();l++){
      delete (*selectedJets)[l];
    }
    selectedLeptons->clear();
    selectedLeptons->shrink_to_fit();
    delete selectedLeptons;
    selectedJets->clear();
    selectedJets->shrink_to_fit();
    delete selectedJets;
//  cout<<ev_event<<endl;
      // if (Cut(ientry) < 0) continue;
  } //end of event loop

cout<<"from "<<ntr<<" evnets, "<<nAccept<<" events are accepted"<<endl;
tree_out.Write() ;
file_out.Close() ;
}
