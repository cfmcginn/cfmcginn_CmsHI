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

Float_t trkCorrLow_[40];

//Track Tree Variables

const int MAXTRKS = 2302; //From SetupTrackTree.h
Int_t nTrk_;
Float_t trkPt_[MAXTRKS];
Float_t trkPtCorr_[MAXTRKS];
Float_t trkPtFact_[MAXTRKS];
Float_t trkPhi_[MAXTRKS];
Float_t trkEta_[MAXTRKS];
Float_t trkLeadDelPhi_[MAXTRKS];

Float_t rImbProjF_;
Float_t rImbProjH_;
Float_t rImbProjL_;
Float_t rImbPerpF_;
Float_t rImbPerpH_;
Float_t rImbPerpL_;

Float_t rImbProjFCorr_;
Float_t rImbProjHCorr_;
Float_t rImbProjLCorr_;
Float_t rImbPerpFCorr_;
Float_t rImbPerpHCorr_;
Float_t rImbPerpLCorr_;

Float_t gRImbPerpF_;
Float_t gRImbPerpH_;
Float_t gRImbPerpL_;
Float_t gRImbProjF_;
Float_t gRImbProjH_;
Float_t gRImbProjL_;

Float_t gRImbPerpFCorr_;
Float_t gRImbPerpHCorr_;
Float_t gRImbPerpLCorr_;
Float_t gRImbProjFCorr_;
Float_t gRImbProjHCorr_;
Float_t gRImbProjLCorr_;

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

Float_t gRLeadJtPt_;
Float_t gRLeadJtPhi_;
Float_t gRLeadJtEta_;
Float_t gRSubLeadJtPt_;
Float_t gRSubLeadJtPhi_;
Float_t gRSubLeadJtEta_;

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
Float_t genLeadDelPhi_[MAXGEN];

Float_t gImbProjF_;
Float_t gImbProjH_;
Float_t gImbProjL_;

Float_t gImbPerpF_;
Float_t gImbPerpH_;
Float_t gImbPerpL_;

Float_t rDivGPt_[10];



void SetBranches(bool montecarlo)
{
  //Track Tree Branches
  
  trackTree_p->Branch("nTrk", &nTrk_, "nTrk/I");
  trackTree_p->Branch("trkPt", &trkPt_, "trkPt[nTrk]/F");
  trackTree_p->Branch("trkPtCorr", &trkPtCorr_, "trkPtCorr[nTrk]/F");
  trackTree_p->Branch("trkPtFact", &trkPtFact_, "trkPtFact[nTrk]/F");
  trackTree_p->Branch("trkPhi", &trkPhi_, "trkPhi[nTrk]/F");
  trackTree_p->Branch("trkEta", &trkEta_, "trkEta[nTrk]/F");
  trackTree_p->Branch("trkLeadDelPhi", &trkLeadDelPhi_, "trkLeadDelPhi[nTrk]/F");
  
  trackTree_p->Branch("rImbProjF", &rImbProjF_, "rImbProjF/F");
  trackTree_p->Branch("rImbProjH", &rImbProjH_, "rImbProjH/F");
  trackTree_p->Branch("rImbProjL", &rImbProjL_, "rImbProjL/F");
  trackTree_p->Branch("rImbPerpF", &rImbPerpF_, "rImbPerpF/F");
  trackTree_p->Branch("rImbPerpH", &rImbPerpH_, "rImbPerpH/F");
  trackTree_p->Branch("rImbPerpL", &rImbPerpL_, "rImbPerpL/F");

  trackTree_p->Branch("rImbProjFCorr", &rImbProjFCorr_, "rImbProjFCorr/F");
  trackTree_p->Branch("rImbProjHCorr", &rImbProjHCorr_, "rImbProjHCorr/F");
  trackTree_p->Branch("rImbProjLCorr", &rImbProjLCorr_, "rImbProjLCorr/F");
  trackTree_p->Branch("rImbPerpFCorr", &rImbPerpFCorr_, "rImbPerpFCorr/F");
  trackTree_p->Branch("rImbPerpHCorr", &rImbPerpHCorr_, "rImbPerpHCorr/F");
  trackTree_p->Branch("rImbPerpLCorr", &rImbPerpLCorr_, "rImbPerpLCorr/F");

  if(montecarlo){
    //Track tree branches iff truth avail.
    trackTree_p->Branch("gRImbProjF", &gRImbProjF_, "gRImbProjF/F");
    trackTree_p->Branch("gRImbProjH", &gRImbProjH_, "gRImbProjH/F");
    trackTree_p->Branch("gRImbProjL", &gRImbProjL_, "gRImbProjL/F");    
    trackTree_p->Branch("gRImbPerpF", &gRImbPerpF_, "gRImbPerpF/F");
    trackTree_p->Branch("gRImbPerpH", &gRImbPerpH_, "gRImbPerpH/F");
    trackTree_p->Branch("gRImbPerpL", &gRImbPerpL_, "gRImbPerpL/F");

    trackTree_p->Branch("gRImbProjFCorr", &gRImbProjFCorr_, "gRImbProjFCorr/F");
    trackTree_p->Branch("gRImbProjHCorr", &gRImbProjHCorr_, "gRImbProjHCorr/F");
    trackTree_p->Branch("gRImbProjLCorr", &gRImbProjLCorr_, "gRImbProjLCorr/F");    
    trackTree_p->Branch("gRImbPerpFCorr", &gRImbPerpFCorr_, "gRImbPerpFCorr/F");
    trackTree_p->Branch("gRImbPerpHCorr", &gRImbPerpHCorr_, "gRImbPerpHCorr/F");
    trackTree_p->Branch("gRImbPerpLCorr", &gRImbPerpLCorr_, "gRImbPerpLCorr/F");
  }
  
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

  if(montecarlo){
    //Jet Tree branches iff truth avail.
    jetTree_p->Branch("gLeadJtPt", &gLeadJtPt_, "gLeadJtPt/F");
    jetTree_p->Branch("gLeadJtPhi", &gLeadJtPhi_, "gLeadJtPhi/F");
    jetTree_p->Branch("gLeadJtEta", &gLeadJtEta_, "gLeadJtEta/F");
    jetTree_p->Branch("gSubLeadJtPt", &gSubLeadJtPt_, "gSubLeadJtPt/F");
    jetTree_p->Branch("gSubLeadJtPhi", &gSubLeadJtPhi_, "gSubLeadJtPhi/F");
    jetTree_p->Branch("gSubLeadJtEta", &gSubLeadJtEta_, "gSubLeadJtEta/F");

    jetTree_p->Branch("gRLeadJtPt", &gRLeadJtPt_, "gRLeadJtPt/F");
    jetTree_p->Branch("gRLeadJtPhi", &gRLeadJtPhi_, "gRLeadJtPhi/F");
    jetTree_p->Branch("gRLeadJtEta", &gRLeadJtEta_, "gRLeadJtEta/F");
    jetTree_p->Branch("gRSubLeadJtPt", &gRSubLeadJtPt_, "gRSubLeadJtPt/F");
    jetTree_p->Branch("gRSubLeadJtPhi", &gRSubLeadJtPhi_, "gRSubLeadJtPhi/F");
    jetTree_p->Branch("gRSubLeadJtEta", &gRSubLeadJtEta_, "gRSubLeadJtEta/F");
  }

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
    genTree_p->Branch("genLeadDelPhi", &genLeadDelPhi_, "genLeadDelPhi[nGen]/F");

    genTree_p->Branch("gImbProjF", &gImbProjF_, "gImbProjF/F");
    genTree_p->Branch("gImbProjH", &gImbProjH_, "gImbProjH/F");
    genTree_p->Branch("gImbProjL", &gImbProjL_, "gImbProjL/F");

    genTree_p->Branch("gImbPerpF", &gImbPerpF_, "gImbPerpF/F");
    genTree_p->Branch("gImbPerpH", &gImbPerpH_, "gImbPerpH/F");
    genTree_p->Branch("gImbPerpL", &gImbPerpL_, "gImbPerpL/F");

    genTree_p->Branch("rDivGPt", &rDivGPt_, "rDivGPt[10]/F");
  }
}


void GetBranches(bool montecarlo)
{
  //Track Tree Branches

  trackTree_p->SetBranchAddress("nTrk", &nTrk_);
  trackTree_p->SetBranchAddress("trkPt", &trkPt_);
  trackTree_p->SetBranchAddress("trkPtCorr", &trkPtCorr_);
  trackTree_p->SetBranchAddress("trkPtFact", &trkPtFact_);
  trackTree_p->SetBranchAddress("trkPhi", &trkPhi_);
  trackTree_p->SetBranchAddress("trkEta", &trkEta_);
  trackTree_p->SetBranchAddress("trkLeadDelPhi", &trkLeadDelPhi_);

  trackTree_p->SetBranchAddress("rImbProjF", &rImbProjF_);
  trackTree_p->SetBranchAddress("rImbProjH", &rImbProjH_);
  trackTree_p->SetBranchAddress("rImbProjL", &rImbProjL_);
  trackTree_p->SetBranchAddress("rImbPerpF", &rImbPerpF_);
  trackTree_p->SetBranchAddress("rImbPerpH", &rImbPerpH_);
  trackTree_p->SetBranchAddress("rImbPerpL", &rImbPerpL_);

  trackTree_p->SetBranchAddress("rImbProjFCorr", &rImbProjFCorr_);
  trackTree_p->SetBranchAddress("rImbProjHCorr", &rImbProjHCorr_);
  trackTree_p->SetBranchAddress("rImbProjLCorr", &rImbProjLCorr_);
  trackTree_p->SetBranchAddress("rImbPerpFCorr", &rImbPerpFCorr_);
  trackTree_p->SetBranchAddress("rImbPerpHCorr", &rImbPerpHCorr_);
  trackTree_p->SetBranchAddress("rImbPerpLCorr", &rImbPerpLCorr_);

  if(montecarlo){
    //Track Tree Branches iff. Truth avail.
    trackTree_p->SetBranchAddress("gRImbProjF", &gRImbProjF_);
    trackTree_p->SetBranchAddress("gRImbProjH", &gRImbProjH_);
    trackTree_p->SetBranchAddress("gRImbProjL", &gRImbProjL_);    
    trackTree_p->SetBranchAddress("gRImbPerpF", &gRImbPerpF_);
    trackTree_p->SetBranchAddress("gRImbPerpH", &gRImbPerpH_);
    trackTree_p->SetBranchAddress("gRImbPerpL", &gRImbPerpL_);

    trackTree_p->SetBranchAddress("gRImbProjFCorr", &gRImbProjFCorr_);
    trackTree_p->SetBranchAddress("gRImbProjHCorr", &gRImbProjHCorr_);
    trackTree_p->SetBranchAddress("gRImbProjLCorr", &gRImbProjLCorr_);    
    trackTree_p->SetBranchAddress("gRImbPerpFCorr", &gRImbPerpFCorr_);
    trackTree_p->SetBranchAddress("gRImbPerpHCorr", &gRImbPerpHCorr_);
    trackTree_p->SetBranchAddress("gRImbPerpLCorr", &gRImbPerpLCorr_);
  }

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

  if(montecarlo){
    //Jet Tree Branches iff truth avail.
    jetTree_p->SetBranchAddress("gLeadJtPt", &gLeadJtPt_);
    jetTree_p->SetBranchAddress("gLeadJtPhi", &gLeadJtPhi_);
    jetTree_p->SetBranchAddress("gLeadJtEta", &gLeadJtEta_);
    jetTree_p->SetBranchAddress("gSubLeadJtPt", &gSubLeadJtPt_);
    jetTree_p->SetBranchAddress("gSubLeadJtPhi", &gSubLeadJtPhi_);
    jetTree_p->SetBranchAddress("gSubLeadJtEta", &gSubLeadJtEta_);
    
    jetTree_p->SetBranchAddress("gRLeadJtPt", &gRLeadJtPt_);
    jetTree_p->SetBranchAddress("gRLeadJtPhi", &gRLeadJtPhi_);
    jetTree_p->SetBranchAddress("gRLeadJtEta", &gRLeadJtEta_);
    jetTree_p->SetBranchAddress("gRSubLeadJtPt", &gRSubLeadJtPt_);
    jetTree_p->SetBranchAddress("gRSubLeadJtPhi", &gRSubLeadJtPhi_);
    jetTree_p->SetBranchAddress("gRSubLeadJtEta", &gRSubLeadJtEta_);
  }

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
    genTree_p->SetBranchAddress("genLeadDelPhi", &genLeadDelPhi_);

    genTree_p->SetBranchAddress("gImbProjF", &gImbProjF_);
    genTree_p->SetBranchAddress("gImbProjH", &gImbProjH_);
    genTree_p->SetBranchAddress("gImbProjL", &gImbProjL_);

    genTree_p->SetBranchAddress("gImbPerpF", &gImbPerpF_);
    genTree_p->SetBranchAddress("gImbPerpH", &gImbPerpH_);
    genTree_p->SetBranchAddress("gImbPerpL", &gImbPerpL_);

    genTree_p->SetBranchAddress("rDivGPt", &rDivGPt_);
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


void InitJetVar(bool montecarlo = false)
{
  if(montecarlo){
    gLeadJtPt_ = -10;
    gSubLeadJtPt_ = -10;
    gLeadJtPhi_ = -10;
    gSubLeadJtPhi_ = -10;
    gLeadJtEta_ = -10;
    gSubLeadJtEta_ = -10;
    
    gRLeadJtPt_ = -10;
    gRSubLeadJtPt_ = -10;
    gRLeadJtPhi_ = -10;
    gRSubLeadJtPhi_ = -10;
    gRLeadJtEta_ = -10;
    gRSubLeadJtEta_ = -10;
  }

  rLeadJtPt_ = -10;
  rSubLeadJtPt_ = -10;
  rLeadJtPhi_ = -10;
  rSubLeadJtPhi_ = -10;
  rLeadJtEta_ = -10;
  rSubLeadJtEta_ = -10;
}

void InitProjPerp(bool montecarlo = false)
{
  rImbProjF_ = 0;
  rImbProjH_ = 0;
  rImbProjL_ = 0;
  rImbPerpF_ = 0;
  rImbPerpH_ = 0;
  rImbPerpL_ = 0;

  rImbProjFCorr_ = 0;
  rImbProjHCorr_ = 0;
  rImbProjLCorr_ = 0;
  rImbPerpFCorr_ = 0;
  rImbPerpHCorr_ = 0;
  rImbPerpLCorr_ = 0;

  if(montecarlo){
    gImbProjF_ = 0;
    gImbProjH_ = 0;
    gImbProjL_ = 0;    
    gImbPerpF_ = 0;
    gImbPerpH_ = 0;
    gImbPerpL_ = 0;

    gRImbProjF_ = 0;
    gRImbProjH_ = 0;
    gRImbProjL_ = 0;
    gRImbPerpF_ = 0;
    gRImbPerpH_ = 0;
    gRImbPerpL_ = 0;

    gRImbProjFCorr_ = 0;
    gRImbProjHCorr_ = 0;
    gRImbProjLCorr_ = 0;
    gRImbPerpFCorr_ = 0;
    gRImbPerpHCorr_ = 0;
    gRImbPerpLCorr_ = 0;
  }
}

//Low track correction defined naively from TGraphAsymmErrors Print from .9 to 4.9, .1 bins.
void defTrkCorr(){
  trkCorrLow_[0] = .389839;
  trkCorrLow_[1] = .498488;
  trkCorrLow_[2] = .518702;
  trkCorrLow_[3] = .527675;
  trkCorrLow_[4] = .536142;
  trkCorrLow_[5] = .542951;
  trkCorrLow_[6] = .549862;
  trkCorrLow_[7] = .555522;
  trkCorrLow_[8] = .560996;
  trkCorrLow_[9] = .566915;
  
  trkCorrLow_[10] = .572628;
  trkCorrLow_[11] = .577247;
  trkCorrLow_[12] = .581227;
  trkCorrLow_[13] = .588351;
  trkCorrLow_[14] = .590406;
  trkCorrLow_[15] = .593747;
  trkCorrLow_[16] = .597831;
  trkCorrLow_[17] = .605439;
  trkCorrLow_[18] = .609709;
  trkCorrLow_[19] = .6124;
  
  trkCorrLow_[20] = .616573;
  trkCorrLow_[21] = .625854;
  trkCorrLow_[22] = .622986;
  trkCorrLow_[23] = .631909;
  trkCorrLow_[24] = .635628;
  trkCorrLow_[25] = .630028;
  trkCorrLow_[26] = .647423;
  trkCorrLow_[27] = .636658;
  trkCorrLow_[28] = .645713;
  trkCorrLow_[29] = .651468;
  
  trkCorrLow_[30] = .649657;
  trkCorrLow_[31] = .645197;
  trkCorrLow_[32] = .65497;
  trkCorrLow_[33] = .658405;
  trkCorrLow_[34] = .66165;
  trkCorrLow_[35] = .665341;
  trkCorrLow_[36] = .656652;
  trkCorrLow_[37] = .663836;
  trkCorrLow_[38] = .672488;
  trkCorrLow_[39] = .664461;
}

#endif
