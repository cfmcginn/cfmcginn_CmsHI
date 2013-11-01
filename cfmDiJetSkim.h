//=============================================                                 
// Author: Chris McGinn                                                         
//                                                                              
// DiJet Analysis Class (MC)                                                    
//                                                                              
//=============================================  
#ifndef cfmDiJetSkim_h
#define cfmDiJetSkim_h

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"

enum sampleType {
  kHIDATA, //0
  kHIMC,   //1
  kPPDATA, //2
  kPPMC,   //3
  kPADATA, //4 
  kPAMC    //5
};

TString getSampleName ( sampleType colli) {
  if (colli == kHIDATA) return "pbpbDATA";
  if (colli == kHIMC) return "pbpbMC";
  if (colli == kPPDATA) return "ppDATA";
  if (colli == kPPMC) return "ppMC";
  if (colli == kPADATA) return "ppbDATA";
  if (colli == kPAMC) return "ppbMC";
  return "NULL";
}
TString getSampleName ( int colli) {
  if (colli == kHIDATA) return "pbpbDATA";
  if (colli == kHIMC) return "pbpbMC";
  if (colli == kPPDATA) return "ppDATA";
  if (colli == kPPMC) return "ppMC";
  if (colli == kPADATA) return "ppbDATA";
  if (colli == kPAMC) return "ppbMC";
  return "NULL";
}

TTree* trackTree_p;
TTree* jetTree_p;

//Track Tree Variables

Int_t nTrk_;
Float_t trkPt_;
Float_t trkPhi_;
Float_t trkEta_;

//Jet Tree Variables

const int MAXJETS = 504; //From SetupJetTree.h

Int_t run_;
Int_t evt_;
Int_t lumi_;
Int_t hiBin_;
Int_t nref_;

Float_t jtpt_[MAXJETS];
Float_t jtphi_[MAXJETS];
Float_t jteta_[MAXJETS];
Float_t refpt_[MAXJETS];
Float_t refphi_[MAXJETS];
Float_t refeta_[MAXJETS];

Float_t leadJtPt_;
Float_t leadJtPhi_;
Float_t leadJtEta_;
Float_t subLeadJtPt_;
Float_t subLeadJtPhi_;
Float_t subLeadJtEta_;

void SetBranches(bool montecarlo)
{
  //Track Tree Branches

  /*
  trackTree_p->Branch("nTrk", &nTrk_, "nTrk/I");
  trackTree_p->Branch("trkPt", &trkPt_, "trkPt/F");
  trackTree_p->Branch("trkPhi", &trkPhi_, "trkPhi/F");
  trackTree_p->Branch("trkEta", &trkEta_, "trkEta/F");
  */

  //Jet Tree Branches

  jetTree_p->Branch("run", &run_, "run/I");
  jetTree_p->Branch("evt", &evt_, "evt/I");
  jetTree_p->Branch("lumi", &lumi_, "lumi/I");
  jetTree_p->Branch("hiBin", &hiBin_, "hiBin/I");

  /*
  jetTree_p->Branch("nref", &nref_, "nref/I");
  jetTree_p->Branch("jtpt", &jtpt_, "jtpt/F");
  jetTree_p->Branch("jtphi", &jtphi_, "jtphi/F");
  jetTree_p->Branch("jteta", &jteta_, "jteta/F");
  jetTree_p->Branch("refpt", &refpt_, "refpt/F");
  jetTree_p->Branch("refphi", &refphi_, "refphi/F");
  jetTree_p->Branch("refeta", &refeta_, "refeta/F");
  */

  jetTree_p->Branch("leadJtPt", &leadJtPt_, "leadJtPt/I");
  jetTree_p->Branch("leadJtPhi", &leadJtPhi_, "leadJtPhi/I");
  jetTree_p->Branch("leadJtEta", &leadJtEta_, "leadJtEta/I");
  jetTree_p->Branch("subLeadJtPt", &subLeadJtPt_, "subLeadJtPt/I");
  jetTree_p->Branch("subLeadJtPhi", &subLeadJtPhi_, "subLeadJtPhi/I");
  jetTree_p->Branch("subLeadJtEta", &subLeadJtEta_, "subLeadJtEta/I");

  if(montecarlo){

  }
}


void GetBranches(bool montecarlo)
{
  //Track Tree Branches

  /*
  trackTree_p->SetBranchAddress("nTrk", &nTrk_);
  trackTree_p->SetBranchAddress("trkPt", &trkPt_);
  trackTree_p->SetBranchAddress("trkPhi", &trkPhi_);
  trackTree_p->SetBranchAddress("trkEta", &trkEta_);
  */

  //Jet Tree Branches

  jetTree_p->SetBranchAddress("run", &run_ );
  jetTree_p->SetBranchAddress("evt", &evt_ );
  jetTree_p->SetBranchAddress("lumi", &lumi_ );
  jetTree_p->SetBranchAddress("hiBin", &hiBin_ );

  /*
  jetTree_p->SetBranchAddress("nref", &nref_ );
  jetTree_p->SetBranchAddress("jtpt", &jtpt_ );
  jetTree_p->SetBranchAddress("jtphi", &jtphi_ );
  jetTree_p->SetBranchAddress("jteta", &jteta_ );
  jetTree_p->SetBranchAddress("refpt", &refpt_ );
  jetTree_p->SetBranchAddress("refphi", &refphi_ );
  jetTree_p->SetBranchAddress("refeta", &refeta_ );
  */


  jetTree_p->SetBranchAddress("leadJtPt", &leadJtPt_);
  jetTree_p->SetBranchAddress("leadJtPhi", &leadJtPhi_);
  jetTree_p->SetBranchAddress("leadJtEta", &leadJtEta_);
  jetTree_p->SetBranchAddress("subLeadJtPt", &subLeadJtPt_);
  jetTree_p->SetBranchAddress("subLeadJtPhi", &subLeadJtPhi_);
  jetTree_p->SetBranchAddress("subLeadJtEta", &subLeadJtEta_);

  if(montecarlo){

  }
}


void ReadDiJetSkim(TFile* inFile, bool montecarlo = false)
{
  //  trackTree_p = (TTree*)inFile->Get("trackTree");
  jetTree_p = (TTree*)inFile->Get("jetTree");

  GetBranches(montecarlo);
}


void InitDiJetSkim(bool montecarlo = false)
{
  //  trackTree_p = new TTree("trackTree", "trackTree");
  jetTree_p = new TTree("jetTree", "jetTree");

  SetBranches(montecarlo);
}

#endif
