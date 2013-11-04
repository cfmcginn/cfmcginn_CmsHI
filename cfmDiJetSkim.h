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
TTree* genTree_p;

//Track Tree Variables

const int MAXTRKS = 2302; //From SetupTrackTree.h
Int_t nTrk_;
Float_t trkPt_[MAXTRKS];
Float_t trkPhi_[MAXTRKS];
Float_t trkEta_[MAXTRKS];

//Jet Tree Variables

const int MAXJETS = 504; //From SetupJetTree.h

Int_t run_;
Int_t evt_;
Int_t lumi_;
Int_t hiBin_;
Int_t nref_;

Int_t nJt_;
Float_t jtpt_[MAXJETS];
Float_t jtphi_[MAXJETS];
Float_t jteta_[MAXJETS];
Float_t refpt_[MAXJETS];
Float_t refphi_[MAXJETS];
Float_t refeta_[MAXJETS];

Float_t gLeadJtPt_;
Float_t gLeadJtPhi_;
Float_t gLeadJtEta_;
Float_t gSubLeadJtPt_;
Float_t gSubLeadJtPhi_;
Float_t gSubLeadJtEta_;

Float_t rLeadJtPt_;
Float_t rLeadJtPhi_;
Float_t rLeadJtEta_;
Float_t rSubLeadJtPt_;
Float_t rSubLeadJtPhi_;
Float_t rSubLeadJtEta_;

//Gen Tree Variables

const int MAXGEN = 17893; //From SetupGenParticleTree.h

Int_t nGen_;
Float_t genPt_[MAXGEN];
Float_t genPhi_[MAXGEN];
Float_t genEta_[MAXGEN];

void SetBranches(bool montecarlo)
{
  //Track Tree Branches
  
  trackTree_p->Branch("nTrk", &nTrk_, "nTrk/I");
  trackTree_p->Branch("trkPt", &trkPt_, "trkPt[nTrk]/F");
  trackTree_p->Branch("trkPhi", &trkPhi_, "trkPhi[nTrk]/F");
  trackTree_p->Branch("trkEta", &trkEta_, "trkEta[nTrk]/F");
  
  
  //Jet Tree Branches

  jetTree_p->Branch("run", &run_, "run/I");
  jetTree_p->Branch("evt", &evt_, "evt/I");
  jetTree_p->Branch("lumi", &lumi_, "lumi/I");
  jetTree_p->Branch("hiBin", &hiBin_, "hiBin/I");

  /*  
  jetTree_p->Branch("nref", &nref_, "nref/I");
  jetTree_p->Branch("jtpt", &jtpt_, "jtpt[nJt]/F");
  jetTree_p->Branch("jtphi", &jtphi_, "jtphi[nJt]/F");
  jetTree_p->Branch("jteta", &jteta_, "jteta[nJt]/F");
  jetTree_p->Branch("refpt", &refpt_, "refpt[nJt]/F");
  jetTree_p->Branch("refphi", &refphi_, "refphi[nJt]/F");
  jetTree_p->Branch("refeta", &refeta_, "refeta[nJt]/F");
  */  

  jetTree_p->Branch("gLeadJtPt", &gLeadJtPt_, "gLeadJtPt/F");
  jetTree_p->Branch("gLeadJtPhi", &gLeadJtPhi_, "gLeadJtPhi/F");
  jetTree_p->Branch("gLeadJtEta", &gLeadJtEta_, "gLeadJtEta/F");
  jetTree_p->Branch("gSubLeadJtPt", &gSubLeadJtPt_, "gSubLeadJtPt/F");
  jetTree_p->Branch("gSubLeadJtPhi", &gSubLeadJtPhi_, "gSubLeadJtPhi/F");
  jetTree_p->Branch("gSubLeadJtEta", &gSubLeadJtEta_, "gSubLeadJtEta/F");

  jetTree_p->Branch("rLeadJtPt", &rLeadJtPt_, "rLeadJtPt/F");
  jetTree_p->Branch("rLeadJtPhi", &rLeadJtPhi_, "rLeadJtPhi/F");
  jetTree_p->Branch("rLeadJtEta", &rLeadJtEta_, "rLeadJtEta/F");
  jetTree_p->Branch("rSubLeadJtPt", &rSubLeadJtPt_, "rSubLeadJtPt/F");
  jetTree_p->Branch("rSubLeadJtPhi", &rSubLeadJtPhi_, "rSubLeadJtPhi/F");
  jetTree_p->Branch("rSubLeadJtEta", &rSubLeadJtEta_, "rSubLeadJtEta/F");


  if(montecarlo){
    //Gen Tree Branches

    genTree_p->Branch("nGen", &nGen_, "nGen/I");
    genTree_p->Branch("genPt", &genPt_, "genPt[nGen]/F");
    genTree_p->Branch("genPhi", &genPhi_, "genPhi[nGen]/F");
    genTree_p->Branch("genEta", &genEta_, "genEta[nGen]/F");
  }
}


void GetBranches(bool montecarlo)
{
  //Track Tree Branches

  trackTree_p->SetBranchAddress("nTrk", &nTrk_);
  trackTree_p->SetBranchAddress("trkPt", &trkPt_);
  trackTree_p->SetBranchAddress("trkPhi", &trkPhi_);
  trackTree_p->SetBranchAddress("trkEta", &trkEta_);

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

  jetTree_p->SetBranchAddress("gLeadJtPt", &gLeadJtPt_);
  jetTree_p->SetBranchAddress("gLeadJtPhi", &gLeadJtPhi_);
  jetTree_p->SetBranchAddress("gLeadJtEta", &gLeadJtEta_);
  jetTree_p->SetBranchAddress("gSubLeadJtPt", &gSubLeadJtPt_);
  jetTree_p->SetBranchAddress("gSubLeadJtPhi", &gSubLeadJtPhi_);
  jetTree_p->SetBranchAddress("gSubLeadJtEta", &gSubLeadJtEta_);

  jetTree_p->SetBranchAddress("rLeadJtPt", &rLeadJtPt_);
  jetTree_p->SetBranchAddress("rLeadJtPhi", &rLeadJtPhi_);
  jetTree_p->SetBranchAddress("rLeadJtEta", &rLeadJtEta_);
  jetTree_p->SetBranchAddress("rSubLeadJtPt", &rSubLeadJtPt_);
  jetTree_p->SetBranchAddress("rSubLeadJtPhi", &rSubLeadJtPhi_);
  jetTree_p->SetBranchAddress("rSubLeadJtEta", &rSubLeadJtEta_);

  if(montecarlo){
    //Gen Tree Branches

    genTree_p->SetBranchAddress("nGen", &nGen_);
    genTree_p->SetBranchAddress("genPt", &genPt_);
    genTree_p->SetBranchAddress("genPhi", &genPhi_);
    genTree_p->SetBranchAddress("genEta", &genEta_);
  }
}


void ReadDiJetSkim(TFile* inFile, bool montecarlo = false)
{
  trackTree_p = (TTree*)inFile->Get("trackTree");
  jetTree_p = (TTree*)inFile->Get("jetTree");

  if(montecarlo)
    genTree_p = (TTree*)inFile->Get("genTree");

  GetBranches(montecarlo);
}


void InitDiJetSkim(bool montecarlo = false)
{
  trackTree_p = new TTree("trackTree", "trackTree");
  jetTree_p = new TTree("jetTree", "jetTree");

  if(montecarlo)
    genTree_p = new TTree("genTree", "genTree");

  SetBranches(montecarlo);
}

#endif
