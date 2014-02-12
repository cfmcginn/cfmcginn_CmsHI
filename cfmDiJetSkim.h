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
Float_t trkPtRPF_[MAXTRKS];
Float_t trkPtGRPF_[MAXTRKS];
Float_t trkPtCorr_[MAXTRKS];
Float_t trkPtFact_[MAXTRKS];
Float_t trkPhi_[MAXTRKS];
Float_t trkEta_[MAXTRKS];
Float_t trkRMin_[MAXTRKS];
Float_t trkPFLeadDelPhi_[MAXTRKS];
Float_t trkCaloLeadDelPhi_[MAXTRKS];

Float_t rPFImbProjF_;
Float_t rPFImbProjH_;
Float_t rPFImbProjL_;
Float_t rPFImbPerpF_;
Float_t rPFImbPerpH_;
Float_t rPFImbPerpL_;

Float_t rPFImbProj5_1_;
Float_t rPFImbProj1_2_;
Float_t rPFImbProj2_4_;
Float_t rPFImbProj4_8_;
Float_t rPFImbProj8_100_;

Float_t rPFImbProjFCorr_;
Float_t rPFImbProjHCorr_;
Float_t rPFImbProjLCorr_;
Float_t rPFImbPerpFCorr_;
Float_t rPFImbPerpHCorr_;
Float_t rPFImbPerpLCorr_;

Float_t rPFImbProjCorr5_1_;
Float_t rPFImbProjCorr1_2_;
Float_t rPFImbProjCorr2_4_;
Float_t rPFImbProjCorr4_8_;
Float_t rPFImbProjCorr8_100_;

Float_t rCaloImbProjF_;
Float_t rCaloImbProjH_;
Float_t rCaloImbProjL_;
Float_t rCaloImbPerpF_;
Float_t rCaloImbPerpH_;
Float_t rCaloImbPerpL_;

Float_t rCaloImbProjFCorr_;
Float_t rCaloImbProjHCorr_;
Float_t rCaloImbProjLCorr_;
Float_t rCaloImbPerpFCorr_;
Float_t rCaloImbPerpHCorr_;
Float_t rCaloImbPerpLCorr_;

Float_t gRPFImbPerpF_;
Float_t gRPFImbPerpH_;
Float_t gRPFImbPerpL_;
Float_t gRPFImbProjF_;
Float_t gRPFImbProjH_;
Float_t gRPFImbProjL_;

Float_t gRPFImbProj5_1_;
Float_t gRPFImbProj1_2_;
Float_t gRPFImbProj2_4_;
Float_t gRPFImbProj4_8_;
Float_t gRPFImbProj8_100_;

Float_t gRPFImbPerpFCorr_;
Float_t gRPFImbPerpHCorr_;
Float_t gRPFImbPerpLCorr_;
Float_t gRPFImbProjFCorr_;
Float_t gRPFImbProjHCorr_;
Float_t gRPFImbProjLCorr_;

Float_t gRPFImbProjCorr5_1_;
Float_t gRPFImbProjCorr1_2_;
Float_t gRPFImbProjCorr2_4_;
Float_t gRPFImbProjCorr4_8_;
Float_t gRPFImbProjCorr8_100_;

//Jet Tree Variables

const int MAXJETS = 504; //From SetupJetTree.h

Int_t run_;
Int_t evt_;
Int_t lumi_;
Int_t hiBin_;
Int_t nref_;

Int_t nJtPF_;
Float_t jtptPF_[MAXJETS];
Float_t jtphiPF_[MAXJETS];
Float_t jtetaPF_[MAXJETS];
Float_t jtptCalo_[MAXJETS];
Float_t jtphiCalo_[MAXJETS];
Float_t jtetaCalo_[MAXJETS];
Float_t refpt_[MAXJETS];
Float_t refphi_[MAXJETS];
Float_t refeta_[MAXJETS];

Bool_t truthSet_;
Bool_t recoPFSet_;
Bool_t recoCaloSet_;

Float_t gLeadJtPt_;
Float_t gLeadJtPhi_;
Float_t gLeadJtEta_;
Float_t gSubLeadJtPt_;
Float_t gSubLeadJtPhi_;
Float_t gSubLeadJtEta_;

Float_t gRPFLeadJtPt_;
Float_t gRPFLeadJtPhi_;
Float_t gRPFLeadJtEta_;
Float_t gRPFSubLeadJtPt_;
Float_t gRPFSubLeadJtPhi_;
Float_t gRPFSubLeadJtEta_;

Float_t gRCaloLeadJtPt_;
Float_t gRCaloLeadJtPhi_;
Float_t gRCaloLeadJtEta_;
Float_t gRCaloSubLeadJtPt_;
Float_t gRCaloSubLeadJtPhi_;
Float_t gRCaloSubLeadJtEta_;

Float_t rPFLeadJtPt_;
Float_t rPFLeadJtPhi_;
Float_t rPFLeadJtEta_;
Float_t rPFSubLeadJtPt_;
Float_t rPFSubLeadJtPhi_;
Float_t rPFSubLeadJtEta_;

Float_t rCaloLeadJtPt_;
Float_t rCaloLeadJtPhi_;
Float_t rCaloLeadJtEta_;
Float_t rCaloSubLeadJtPt_;
Float_t rCaloSubLeadJtPhi_;
Float_t rCaloSubLeadJtEta_;

//Gen Tree Variables

const int MAXGEN = 50000; //From SetupGenParticleTree.h

Int_t nGen_;
Float_t genPt_[MAXGEN];
Float_t genPtJT_[MAXGEN];
Float_t genPhi_[MAXGEN];
Float_t genEta_[MAXGEN];
Float_t genLeadDelPhi_[MAXGEN];

Float_t gImbProjF_;
Float_t gImbProjH_;
Float_t gImbProjL_;

Float_t gImbPerpF_;
Float_t gImbPerpH_;
Float_t gImbPerpL_;

Float_t gImbProj5_1_;
Float_t gImbProj1_2_;
Float_t gImbProj2_4_;
Float_t gImbProj4_8_;
Float_t gImbProj8_100_;

Float_t rDivGPt_[10];



void SetBranches(bool montecarlo)
{
  //Track Tree Branches
  
  trackTree_p->Branch("nTrk", &nTrk_, "nTrk/I");
  trackTree_p->Branch("trkPt", &trkPt_, "trkPt[nTrk]/F");
  trackTree_p->Branch("trkPtRPF", &trkPtRPF_, "trkPtRPF[nTrk]/F");
  trackTree_p->Branch("trkPtGRPF", &trkPtGRPF_, "trkPtGRPF[nTrk]/F");
  trackTree_p->Branch("trkPtCorr", &trkPtCorr_, "trkPtCorr[nTrk]/F");
  trackTree_p->Branch("trkPtFact", &trkPtFact_, "trkPtFact[nTrk]/F");
  trackTree_p->Branch("trkPhi", &trkPhi_, "trkPhi[nTrk]/F");
  trackTree_p->Branch("trkEta", &trkEta_, "trkEta[nTrk]/F");
  trackTree_p->Branch("trkRMin", &trkRMin_, "trkRMin[nTrk]/F");
  trackTree_p->Branch("trkPFLeadDelPhi", &trkPFLeadDelPhi_, "trkPFLeadDelPhi[nTrk]/F");
  trackTree_p->Branch("trkCaloLeadDelPhi", &trkCaloLeadDelPhi_, "trkCaloLeadDelPhi[nTrk]/F");
  
  trackTree_p->Branch("rPFImbProjF", &rPFImbProjF_, "rPFImbProjF/F");
  trackTree_p->Branch("rPFImbProjH", &rPFImbProjH_, "rPFImbProjH/F");
  trackTree_p->Branch("rPFImbProjL", &rPFImbProjL_, "rPFImbProjL/F");
  trackTree_p->Branch("rPFImbPerpF", &rPFImbPerpF_, "rPFImbPerpF/F");
  trackTree_p->Branch("rPFImbPerpH", &rPFImbPerpH_, "rPFImbPerpH/F");
  trackTree_p->Branch("rPFImbPerpL", &rPFImbPerpL_, "rPFImbPerpL/F");

  trackTree_p->Branch("rPFImbProj5_1", &rPFImbProj5_1_, "rPFImbProj5_1/F");
  trackTree_p->Branch("rPFImbProj1_2", &rPFImbProj1_2_, "rPFImbProj1_2/F");
  trackTree_p->Branch("rPFImbProj2_4", &rPFImbProj2_4_, "rPFImbProj2_4/F");
  trackTree_p->Branch("rPFImbProj4_8", &rPFImbProj4_8_, "rPFImbProj4_8/F");
  trackTree_p->Branch("rPFImbProj8_100", &rPFImbProj8_100_, "rPFImbProj8_100/F");

  trackTree_p->Branch("rPFImbProjFCorr", &rPFImbProjFCorr_, "rPFImbProjFCorr/F");
  trackTree_p->Branch("rPFImbProjHCorr", &rPFImbProjHCorr_, "rPFImbProjHCorr/F");
  trackTree_p->Branch("rPFImbProjLCorr", &rPFImbProjLCorr_, "rPFImbProjLCorr/F");
  trackTree_p->Branch("rPFImbPerpFCorr", &rPFImbPerpFCorr_, "rPFImbPerpFCorr/F");
  trackTree_p->Branch("rPFImbPerpHCorr", &rPFImbPerpHCorr_, "rPFImbPerpHCorr/F");
  trackTree_p->Branch("rPFImbPerpLCorr", &rPFImbPerpLCorr_, "rPFImbPerpLCorr/F");

  trackTree_p->Branch("rPFImbProjCorr5_1", &rPFImbProjCorr5_1_, "rPFImbProjCorr5_1/F");
  trackTree_p->Branch("rPFImbProjCorr1_2", &rPFImbProjCorr1_2_, "rPFImbProjCorr1_2/F");
  trackTree_p->Branch("rPFImbProjCorr2_4", &rPFImbProjCorr2_4_, "rPFImbProjCorr2_4/F");
  trackTree_p->Branch("rPFImbProjCorr4_8", &rPFImbProjCorr4_8_, "rPFImbProjCorr4_8/F");
  trackTree_p->Branch("rPFImbProjCorr8_100", &rPFImbProjCorr8_100_, "rPFImbProjCorr8_100/F");

  trackTree_p->Branch("rCaloImbProjF", &rCaloImbProjF_, "rCaloImbProjF/F");
  trackTree_p->Branch("rCaloImbProjH", &rCaloImbProjH_, "rCaloImbProjH/F");
  trackTree_p->Branch("rCaloImbProjL", &rCaloImbProjL_, "rCaloImbProjL/F");
  trackTree_p->Branch("rCaloImbPerpF", &rCaloImbPerpF_, "rCaloImbPerpF/F");
  trackTree_p->Branch("rCaloImbPerpH", &rCaloImbPerpH_, "rCaloImbPerpH/F");
  trackTree_p->Branch("rCaloImbPerpL", &rCaloImbPerpL_, "rCaloImbPerpL/F");

  trackTree_p->Branch("rCaloImbProjFCorr", &rCaloImbProjFCorr_, "rCaloImbProjFCorr/F");
  trackTree_p->Branch("rCaloImbProjHCorr", &rCaloImbProjHCorr_, "rCaloImbProjHCorr/F");
  trackTree_p->Branch("rCaloImbProjLCorr", &rCaloImbProjLCorr_, "rCaloImbProjLCorr/F");
  trackTree_p->Branch("rCaloImbPerpFCorr", &rCaloImbPerpFCorr_, "rCaloImbPerpFCorr/F");
  trackTree_p->Branch("rCaloImbPerpHCorr", &rCaloImbPerpHCorr_, "rCaloImbPerpHCorr/F");
  trackTree_p->Branch("rCaloImbPerpLCorr", &rCaloImbPerpLCorr_, "rCaloImbPerpLCorr/F");

  if(montecarlo){
    //Track tree branches iff truth avail.
    trackTree_p->Branch("gRPFImbProjF", &gRPFImbProjF_, "gRPFImbProjF/F");
    trackTree_p->Branch("gRPFImbProjH", &gRPFImbProjH_, "gRPFImbProjH/F");
    trackTree_p->Branch("gRPFImbProjL", &gRPFImbProjL_, "gRPFImbProjL/F");    
    trackTree_p->Branch("gRPFImbPerpF", &gRPFImbPerpF_, "gRPFImbPerpF/F");
    trackTree_p->Branch("gRPFImbPerpH", &gRPFImbPerpH_, "gRPFImbPerpH/F");
    trackTree_p->Branch("gRPFImbPerpL", &gRPFImbPerpL_, "gRPFImbPerpL/F");

    trackTree_p->Branch("gRPFImbProj5_1", &gRPFImbProj5_1_, "gRPFImbProj5_1/F");
    trackTree_p->Branch("gRPFImbProj1_2", &gRPFImbProj1_2_, "gRPFImbProj1_2/F");
    trackTree_p->Branch("gRPFImbProj2_4", &gRPFImbProj2_4_, "gRPFImbProj2_4/F");
    trackTree_p->Branch("gRPFImbProj4_8", &gRPFImbProj4_8_, "gRPFImbProj4_8/F");
    trackTree_p->Branch("gRPFImbProj8_100", &gRPFImbProj8_100_, "gRPFImbProj8_100/F");

    trackTree_p->Branch("gRPFImbProjFCorr", &gRPFImbProjFCorr_, "gRPFImbProjFCorr/F");
    trackTree_p->Branch("gRPFImbProjHCorr", &gRPFImbProjHCorr_, "gRPFImbProjHCorr/F");
    trackTree_p->Branch("gRPFImbProjLCorr", &gRPFImbProjLCorr_, "gRPFImbProjLCorr/F");    
    trackTree_p->Branch("gRPFImbPerpFCorr", &gRPFImbPerpFCorr_, "gRPFImbPerpFCorr/F");
    trackTree_p->Branch("gRPFImbPerpHCorr", &gRPFImbPerpHCorr_, "gRPFImbPerpHCorr/F");
    trackTree_p->Branch("gRPFImbPerpLCorr", &gRPFImbPerpLCorr_, "gRPFImbPerpLCorr/F");

    trackTree_p->Branch("gRPFImbProjCorr5_1", &gRPFImbProjCorr5_1_, "gRPFImbProjCorr5_1/F");
    trackTree_p->Branch("gRPFImbProjCorr1_2", &gRPFImbProjCorr1_2_, "gRPFImbProjCorr1_2/F");
    trackTree_p->Branch("gRPFImbProjCorr2_4", &gRPFImbProjCorr2_4_, "gRPFImbProjCorr2_4/F");
    trackTree_p->Branch("gRPFImbProjCorr4_8", &gRPFImbProjCorr4_8_, "gRPFImbProjCorr4_8/F");
    trackTree_p->Branch("gRPFImbProjCorr8_100", &gRPFImbProjCorr8_100_, "gRPFImbProjCorr8_100/F");
  }
  
  //Jet Tree Branches

  jetTree_p->Branch("run", &run_, "run/I");
  jetTree_p->Branch("evt", &evt_, "evt/I");
  jetTree_p->Branch("lumi", &lumi_, "lumi/I");
  jetTree_p->Branch("hiBin", &hiBin_, "hiBin/I");

  /*  
  jetTree_p->Branch("nref", &nref_, "nref/I");
  jetTree_p->Branch("jtptPF", &jtptPF_, "jtptPF[nJt]/F");
  jetTree_p->Branch("jtphiPF", &jtphiPF_, "jtphiPF[nJt]/F");
  jetTree_p->Branch("jtetaPF", &jtetaPF_, "jtetaPF[nJt]/F");
  jetTree_p->Branch("jtptCalo", &jtptCalo_, "jtptCalo[nJt]/F");
  jetTree_p->Branch("jtphiCalo", &jtphiCalo_, "jtphiCalo[nJt]/F");
  jetTree_p->Branch("jtetaCalo", &jtetaCalo_, "jtetaCalo[nJt]/F");
  jetTree_p->Branch("refpt", &refpt_, "refpt[nJt]/F");
  jetTree_p->Branch("refphi", &refphi_, "refphi[nJt]/F");
  jetTree_p->Branch("refeta", &refeta_, "refeta[nJt]/F");
  */  

  jetTree_p->Branch("recoPFSet", &recoPFSet_, "recoPFSet/O");
  jetTree_p->Branch("recoCaloSet", &recoCaloSet_, "recoCaloSet/O");

  if(montecarlo){
    //Jet Tree branches iff truth avail.
    jetTree_p->Branch("truthSet", &truthSet_, "truthSet/O");

    jetTree_p->Branch("gLeadJtPt", &gLeadJtPt_, "gLeadJtPt/F");
    jetTree_p->Branch("gLeadJtPhi", &gLeadJtPhi_, "gLeadJtPhi/F");
    jetTree_p->Branch("gLeadJtEta", &gLeadJtEta_, "gLeadJtEta/F");
    jetTree_p->Branch("gSubLeadJtPt", &gSubLeadJtPt_, "gSubLeadJtPt/F");
    jetTree_p->Branch("gSubLeadJtPhi", &gSubLeadJtPhi_, "gSubLeadJtPhi/F");
    jetTree_p->Branch("gSubLeadJtEta", &gSubLeadJtEta_, "gSubLeadJtEta/F");

    jetTree_p->Branch("gRPFLeadJtPt", &gRPFLeadJtPt_, "gRPFLeadJtPt/F");
    jetTree_p->Branch("gRPFLeadJtPhi", &gRPFLeadJtPhi_, "gRPFLeadJtPhi/F");
    jetTree_p->Branch("gRPFLeadJtEta", &gRPFLeadJtEta_, "gRPFLeadJtEta/F");
    jetTree_p->Branch("gRPFSubLeadJtPt", &gRPFSubLeadJtPt_, "gRPFSubLeadJtPt/F");
    jetTree_p->Branch("gRPFSubLeadJtPhi", &gRPFSubLeadJtPhi_, "gRPFSubLeadJtPhi/F");
    jetTree_p->Branch("gRPFSubLeadJtEta", &gRPFSubLeadJtEta_, "gRPFSubLeadJtEta/F");

    jetTree_p->Branch("gRCaloLeadJtPt", &gRCaloLeadJtPt_, "gRCaloLeadJtPt/F");
    jetTree_p->Branch("gRCaloLeadJtPhi", &gRCaloLeadJtPhi_, "gRCaloLeadJtPhi/F");
    jetTree_p->Branch("gRCaloLeadJtEta", &gRCaloLeadJtEta_, "gRCaloLeadJtEta/F");
    jetTree_p->Branch("gRCaloSubLeadJtPt", &gRCaloSubLeadJtPt_, "gRCaloSubLeadJtPt/F");
    jetTree_p->Branch("gRCaloSubLeadJtPhi", &gRCaloSubLeadJtPhi_, "gRCaloSubLeadJtPhi/F");
    jetTree_p->Branch("gRCaloSubLeadJtEta", &gRCaloSubLeadJtEta_, "gRCaloSubLeadJtEta/F");
  }


  jetTree_p->Branch("rPFLeadJtPt", &rPFLeadJtPt_, "rPFLeadJtPt/F");
  jetTree_p->Branch("rPFLeadJtPhi", &rPFLeadJtPhi_, "rPFLeadJtPhi/F");
  jetTree_p->Branch("rPFLeadJtEta", &rPFLeadJtEta_, "rPFLeadJtEta/F");
  jetTree_p->Branch("rPFSubLeadJtPt", &rPFSubLeadJtPt_, "rPFSubLeadJtPt/F");
  jetTree_p->Branch("rPFSubLeadJtPhi", &rPFSubLeadJtPhi_, "rPFSubLeadJtPhi/F");
  jetTree_p->Branch("rPFSubLeadJtEta", &rPFSubLeadJtEta_, "rPFSubLeadJtEta/F");

  jetTree_p->Branch("rCaloLeadJtPt", &rCaloLeadJtPt_, "rCaloLeadJtPt/F");
  jetTree_p->Branch("rCaloLeadJtPhi", &rCaloLeadJtPhi_, "rCaloLeadJtPhi/F");
  jetTree_p->Branch("rCaloLeadJtEta", &rCaloLeadJtEta_, "rCaloLeadJtEta/F");
  jetTree_p->Branch("rCaloSubLeadJtPt", &rCaloSubLeadJtPt_, "rCaloSubLeadJtPt/F");
  jetTree_p->Branch("rCaloSubLeadJtPhi", &rCaloSubLeadJtPhi_, "rCaloSubLeadJtPhi/F");
  jetTree_p->Branch("rCaloSubLeadJtEta", &rCaloSubLeadJtEta_, "rCaloSubLeadJtEta/F");


  if(montecarlo){
    //Gen Tree Branches

    genTree_p->Branch("nGen", &nGen_, "nGen/I");
    genTree_p->Branch("genPt", &genPt_, "genPt[nGen]/F");
    genTree_p->Branch("genPtJT", &genPtJT_, "genPtJT[nGen]/F");
    genTree_p->Branch("genPhi", &genPhi_, "genPhi[nGen]/F");
    genTree_p->Branch("genEta", &genEta_, "genEta[nGen]/F");
    genTree_p->Branch("genLeadDelPhi", &genLeadDelPhi_, "genLeadDelPhi[nGen]/F");

    genTree_p->Branch("gImbProjF", &gImbProjF_, "gImbProjF/F");
    genTree_p->Branch("gImbProjH", &gImbProjH_, "gImbProjH/F");
    genTree_p->Branch("gImbProjL", &gImbProjL_, "gImbProjL/F");

    genTree_p->Branch("gImbPerpF", &gImbPerpF_, "gImbPerpF/F");
    genTree_p->Branch("gImbPerpH", &gImbPerpH_, "gImbPerpH/F");
    genTree_p->Branch("gImbPerpL", &gImbPerpL_, "gImbPerpL/F");

    genTree_p->Branch("gImbProj5_1", &gImbProj5_1_, "gImbProj5_1/F");
    genTree_p->Branch("gImbProj1_2", &gImbProj1_2_, "gImbProj1_2/F");
    genTree_p->Branch("gImbProj2_4", &gImbProj2_4_, "gImbProj2_4/F");
    genTree_p->Branch("gImbProj4_8", &gImbProj4_8_, "gImbProj4_8/F");
    genTree_p->Branch("gImbProj8_100", &gImbProj8_100_, "gImbProj8_100/F");

    genTree_p->Branch("rDivGPt", &rDivGPt_, "rDivGPt[10]/F");
  }
}


void GetBranches(bool montecarlo)
{
  //Track Tree Branches

  trackTree_p->SetBranchAddress("nTrk", &nTrk_);
  trackTree_p->SetBranchAddress("trkPt", &trkPt_);
  trackTree_p->SetBranchAddress("trkPtRPF", &trkPtRPF_);
  trackTree_p->SetBranchAddress("trkPtGRPF", &trkPtGRPF_);
  trackTree_p->SetBranchAddress("trkPtCorr", &trkPtCorr_);
  trackTree_p->SetBranchAddress("trkPtFact", &trkPtFact_);
  trackTree_p->SetBranchAddress("trkPhi", &trkPhi_);
  trackTree_p->SetBranchAddress("trkEta", &trkEta_);
  trackTree_p->SetBranchAddress("trkRMin", &trkRMin_);
  trackTree_p->SetBranchAddress("trkPFLeadDelPhi", &trkPFLeadDelPhi_);
  trackTree_p->SetBranchAddress("trkCaloLeadDelPhi", &trkCaloLeadDelPhi_);

  trackTree_p->SetBranchAddress("rPFImbProjF", &rPFImbProjF_);
  trackTree_p->SetBranchAddress("rPFImbProjH", &rPFImbProjH_);
  trackTree_p->SetBranchAddress("rPFImbProjL", &rPFImbProjL_);
  trackTree_p->SetBranchAddress("rPFImbPerpF", &rPFImbPerpF_);
  trackTree_p->SetBranchAddress("rPFImbPerpH", &rPFImbPerpH_);
  trackTree_p->SetBranchAddress("rPFImbPerpL", &rPFImbPerpL_);

  trackTree_p->SetBranchAddress("rPFImbProj5_1", &rPFImbProj5_1_);
  trackTree_p->SetBranchAddress("rPFImbProj1_2", &rPFImbProj1_2_);
  trackTree_p->SetBranchAddress("rPFImbProj2_4", &rPFImbProj2_4_);
  trackTree_p->SetBranchAddress("rPFImbProj4_8", &rPFImbProj4_8_);
  trackTree_p->SetBranchAddress("rPFImbProj8_100", &rPFImbProj8_100_);

  trackTree_p->SetBranchAddress("rPFImbProjFCorr", &rPFImbProjFCorr_);
  trackTree_p->SetBranchAddress("rPFImbProjHCorr", &rPFImbProjHCorr_);
  trackTree_p->SetBranchAddress("rPFImbProjLCorr", &rPFImbProjLCorr_);
  trackTree_p->SetBranchAddress("rPFImbPerpFCorr", &rPFImbPerpFCorr_);
  trackTree_p->SetBranchAddress("rPFImbPerpHCorr", &rPFImbPerpHCorr_);
  trackTree_p->SetBranchAddress("rPFImbPerpLCorr", &rPFImbPerpLCorr_);

  trackTree_p->SetBranchAddress("rPFImbProjCorr5_1", &rPFImbProjCorr5_1_);
  trackTree_p->SetBranchAddress("rPFImbProjCorr1_2", &rPFImbProjCorr1_2_);
  trackTree_p->SetBranchAddress("rPFImbProjCorr2_4", &rPFImbProjCorr2_4_);
  trackTree_p->SetBranchAddress("rPFImbProjCorr4_8", &rPFImbProjCorr4_8_);
  trackTree_p->SetBranchAddress("rPFImbProjCorr8_100", &rPFImbProjCorr8_100_);

  trackTree_p->SetBranchAddress("rCaloImbProjF", &rCaloImbProjF_);
  trackTree_p->SetBranchAddress("rCaloImbProjH", &rCaloImbProjH_);
  trackTree_p->SetBranchAddress("rCaloImbProjL", &rCaloImbProjL_);
  trackTree_p->SetBranchAddress("rCaloImbPerpF", &rCaloImbPerpF_);
  trackTree_p->SetBranchAddress("rCaloImbPerpH", &rCaloImbPerpH_);
  trackTree_p->SetBranchAddress("rCaloImbPerpL", &rCaloImbPerpL_);

  trackTree_p->SetBranchAddress("rCaloImbProjFCorr", &rCaloImbProjFCorr_);
  trackTree_p->SetBranchAddress("rCaloImbProjHCorr", &rCaloImbProjHCorr_);
  trackTree_p->SetBranchAddress("rCaloImbProjLCorr", &rCaloImbProjLCorr_);
  trackTree_p->SetBranchAddress("rCaloImbPerpFCorr", &rCaloImbPerpFCorr_);
  trackTree_p->SetBranchAddress("rCaloImbPerpHCorr", &rCaloImbPerpHCorr_);
  trackTree_p->SetBranchAddress("rCaloImbPerpLCorr", &rCaloImbPerpLCorr_);

  if(montecarlo){
    //Track Tree Branches iff. Truth avail.
    trackTree_p->SetBranchAddress("gRPFImbProjF", &gRPFImbProjF_);
    trackTree_p->SetBranchAddress("gRPFImbProjH", &gRPFImbProjH_);
    trackTree_p->SetBranchAddress("gRPFImbProjL", &gRPFImbProjL_);    
    trackTree_p->SetBranchAddress("gRPFImbPerpF", &gRPFImbPerpF_);
    trackTree_p->SetBranchAddress("gRPFImbPerpH", &gRPFImbPerpH_);
    trackTree_p->SetBranchAddress("gRPFImbPerpL", &gRPFImbPerpL_);

    trackTree_p->SetBranchAddress("gRPFImbProj5_1", &gRPFImbProj5_1_);
    trackTree_p->SetBranchAddress("gRPFImbProj1_2", &gRPFImbProj1_2_);
    trackTree_p->SetBranchAddress("gRPFImbProj2_4", &gRPFImbProj2_4_);
    trackTree_p->SetBranchAddress("gRPFImbProj4_8", &gRPFImbProj4_8_);
    trackTree_p->SetBranchAddress("gRPFImbProj8_100", &gRPFImbProj8_100_);

    trackTree_p->SetBranchAddress("gRPFImbProjFCorr", &gRPFImbProjFCorr_);
    trackTree_p->SetBranchAddress("gRPFImbProjHCorr", &gRPFImbProjHCorr_);
    trackTree_p->SetBranchAddress("gRPFImbProjLCorr", &gRPFImbProjLCorr_);    
    trackTree_p->SetBranchAddress("gRPFImbPerpFCorr", &gRPFImbPerpFCorr_);
    trackTree_p->SetBranchAddress("gRPFImbPerpHCorr", &gRPFImbPerpHCorr_);
    trackTree_p->SetBranchAddress("gRPFImbPerpLCorr", &gRPFImbPerpLCorr_);

    trackTree_p->SetBranchAddress("gRPFImbProjCorr5_1", &gRPFImbProjCorr5_1_);
    trackTree_p->SetBranchAddress("gRPFImbProjCorr1_2", &gRPFImbProjCorr1_2_);
    trackTree_p->SetBranchAddress("gRPFImbProjCorr2_4", &gRPFImbProjCorr2_4_);
    trackTree_p->SetBranchAddress("gRPFImbProjCorr4_8", &gRPFImbProjCorr4_8_);
    trackTree_p->SetBranchAddress("gRPFImbProjCorr8_100", &gRPFImbProjCorr8_100_);
  }

  //Jet Tree Branches

  jetTree_p->SetBranchAddress("run", &run_ );
  jetTree_p->SetBranchAddress("evt", &evt_ );
  jetTree_p->SetBranchAddress("lumi", &lumi_ );
  jetTree_p->SetBranchAddress("hiBin", &hiBin_ );

  /*
  jetTree_p->SetBranchAddress("nref", &nref_ );
  jetTree_p->SetBranchAddress("jtptPF", &jtptPF_ );
  jetTree_p->SetBranchAddress("jtphiPF", &jtphiPF_ );
  jetTree_p->SetBranchAddress("jtetaPF", &jtetaPF_ );
  jetTree_p->SetBranchAddress("jtptCalo", &jtptCalo_ );
  jetTree_p->SetBranchAddress("jtphiCalo", &jtphiCalo_ );
  jetTree_p->SetBranchAddress("jtetaCalo", &jtetaCalo_ );
  jetTree_p->SetBranchAddress("refpt", &refpt_ );
  jetTree_p->SetBranchAddress("refphi", &refphi_ );
  jetTree_p->SetBranchAddress("refeta", &refeta_ );
  */  

  jetTree_p->SetBranchAddress("recoPFSet", &recoPFSet_);
  jetTree_p->SetBranchAddress("recoCaloSet", &recoCaloSet_);

  if(montecarlo){
    //Jet Tree Branches iff truth avail.
    jetTree_p->SetBranchAddress("truthSet", &truthSet_);

    jetTree_p->SetBranchAddress("gLeadJtPt", &gLeadJtPt_);
    jetTree_p->SetBranchAddress("gLeadJtPhi", &gLeadJtPhi_);
    jetTree_p->SetBranchAddress("gLeadJtEta", &gLeadJtEta_);
    jetTree_p->SetBranchAddress("gSubLeadJtPt", &gSubLeadJtPt_);
    jetTree_p->SetBranchAddress("gSubLeadJtPhi", &gSubLeadJtPhi_);
    jetTree_p->SetBranchAddress("gSubLeadJtEta", &gSubLeadJtEta_);
    
    jetTree_p->SetBranchAddress("gRPFLeadJtPt", &gRPFLeadJtPt_);
    jetTree_p->SetBranchAddress("gRPFLeadJtPhi", &gRPFLeadJtPhi_);
    jetTree_p->SetBranchAddress("gRPFLeadJtEta", &gRPFLeadJtEta_);
    jetTree_p->SetBranchAddress("gRPFSubLeadJtPt", &gRPFSubLeadJtPt_);
    jetTree_p->SetBranchAddress("gRPFSubLeadJtPhi", &gRPFSubLeadJtPhi_);
    jetTree_p->SetBranchAddress("gRPFSubLeadJtEta", &gRPFSubLeadJtEta_);

    jetTree_p->SetBranchAddress("gRCaloLeadJtPt", &gRCaloLeadJtPt_);
    jetTree_p->SetBranchAddress("gRCaloLeadJtPhi", &gRCaloLeadJtPhi_);
    jetTree_p->SetBranchAddress("gRCaloLeadJtEta", &gRCaloLeadJtEta_);
    jetTree_p->SetBranchAddress("gRCaloSubLeadJtPt", &gRCaloSubLeadJtPt_);
    jetTree_p->SetBranchAddress("gRCaloSubLeadJtPhi", &gRCaloSubLeadJtPhi_);
    jetTree_p->SetBranchAddress("gRCaloSubLeadJtEta", &gRCaloSubLeadJtEta_);
  }

  jetTree_p->SetBranchAddress("rPFLeadJtPt", &rPFLeadJtPt_);
  jetTree_p->SetBranchAddress("rPFLeadJtPhi", &rPFLeadJtPhi_);
  jetTree_p->SetBranchAddress("rPFLeadJtEta", &rPFLeadJtEta_);
  jetTree_p->SetBranchAddress("rPFSubLeadJtPt", &rPFSubLeadJtPt_);
  jetTree_p->SetBranchAddress("rPFSubLeadJtPhi", &rPFSubLeadJtPhi_);
  jetTree_p->SetBranchAddress("rPFSubLeadJtEta", &rPFSubLeadJtEta_);

  jetTree_p->SetBranchAddress("rCaloLeadJtPt", &rCaloLeadJtPt_);
  jetTree_p->SetBranchAddress("rCaloLeadJtPhi", &rCaloLeadJtPhi_);
  jetTree_p->SetBranchAddress("rCaloLeadJtEta", &rCaloLeadJtEta_);
  jetTree_p->SetBranchAddress("rCaloSubLeadJtPt", &rCaloSubLeadJtPt_);
  jetTree_p->SetBranchAddress("rCaloSubLeadJtPhi", &rCaloSubLeadJtPhi_);
  jetTree_p->SetBranchAddress("rCaloSubLeadJtEta", &rCaloSubLeadJtEta_);

  if(montecarlo){
    //Gen Tree Branches

    genTree_p->SetBranchAddress("nGen", &nGen_);
    genTree_p->SetBranchAddress("genPt", &genPt_);
    genTree_p->SetBranchAddress("genPtJT", &genPtJT_);
    genTree_p->SetBranchAddress("genPhi", &genPhi_);
    genTree_p->SetBranchAddress("genEta", &genEta_);
    genTree_p->SetBranchAddress("genLeadDelPhi", &genLeadDelPhi_);

    genTree_p->SetBranchAddress("gImbProjF", &gImbProjF_);
    genTree_p->SetBranchAddress("gImbProjH", &gImbProjH_);
    genTree_p->SetBranchAddress("gImbProjL", &gImbProjL_);

    genTree_p->SetBranchAddress("gImbPerpF", &gImbPerpF_);
    genTree_p->SetBranchAddress("gImbPerpH", &gImbPerpH_);
    genTree_p->SetBranchAddress("gImbPerpL", &gImbPerpL_);

    genTree_p->SetBranchAddress("gImbProj5_1", &gImbProj5_1_);
    genTree_p->SetBranchAddress("gImbProj1_2", &gImbProj1_2_);
    genTree_p->SetBranchAddress("gImbProj2_4", &gImbProj2_4_);
    genTree_p->SetBranchAddress("gImbProj4_8", &gImbProj4_8_);
    genTree_p->SetBranchAddress("gImbProj8_100", &gImbProj8_100_);

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
    truthSet_ = false;

    gLeadJtPt_ = -10;
    gSubLeadJtPt_ = -10;
    gLeadJtPhi_ = -10;
    gSubLeadJtPhi_ = -10;
    gLeadJtEta_ = -10;
    gSubLeadJtEta_ = -10;
    
    gRPFLeadJtPt_ = -10;
    gRPFSubLeadJtPt_ = -10;
    gRPFLeadJtPhi_ = -10;
    gRPFSubLeadJtPhi_ = -10;
    gRPFLeadJtEta_ = -10;
    gRPFSubLeadJtEta_ = -10;

    gRCaloLeadJtPt_ = -10;
    gRCaloSubLeadJtPt_ = -10;
    gRCaloLeadJtPhi_ = -10;
    gRCaloSubLeadJtPhi_ = -10;
    gRCaloLeadJtEta_ = -10;
    gRCaloSubLeadJtEta_ = -10;
  }

  recoPFSet_ = false;
  recoCaloSet_ = false;

  rPFLeadJtPt_ = -10;
  rPFSubLeadJtPt_ = -10;
  rPFLeadJtPhi_ = -10;
  rPFSubLeadJtPhi_ = -10;
  rPFLeadJtEta_ = -10;
  rPFSubLeadJtEta_ = -10;

  rCaloLeadJtPt_ = -10;
  rCaloSubLeadJtPt_ = -10;
  rCaloLeadJtPhi_ = -10;
  rCaloSubLeadJtPhi_ = -10;
  rCaloLeadJtEta_ = -10;
  rCaloSubLeadJtEta_ = -10;
}

void InitProjPerp(bool montecarlo = false)
{
  rPFImbProjF_ = 0;
  rPFImbProjH_ = 0;
  rPFImbProjL_ = 0;
  rPFImbPerpF_ = 0;
  rPFImbPerpH_ = 0;
  rPFImbPerpL_ = 0;

  rPFImbProj5_1_ = 0;
  rPFImbProj1_2_ = 0;
  rPFImbProj2_4_ = 0;
  rPFImbProj4_8_ = 0;
  rPFImbProj8_100_ = 0;

  rPFImbProjFCorr_ = 0;
  rPFImbProjHCorr_ = 0;
  rPFImbProjLCorr_ = 0;
  rPFImbPerpFCorr_ = 0;
  rPFImbPerpHCorr_ = 0;
  rPFImbPerpLCorr_ = 0;

  rPFImbProjCorr5_1_ = 0;
  rPFImbProjCorr1_2_ = 0;
  rPFImbProjCorr2_4_ = 0;
  rPFImbProjCorr4_8_ = 0;
  rPFImbProjCorr8_100_ = 0;

  rCaloImbProjF_ = 0;
  rCaloImbProjH_ = 0;
  rCaloImbProjL_ = 0;
  rCaloImbPerpF_ = 0;
  rCaloImbPerpH_ = 0;
  rCaloImbPerpL_ = 0;

  rCaloImbProjFCorr_ = 0;
  rCaloImbProjHCorr_ = 0;
  rCaloImbProjLCorr_ = 0;
  rCaloImbPerpFCorr_ = 0;
  rCaloImbPerpHCorr_ = 0;
  rCaloImbPerpLCorr_ = 0;

  if(montecarlo){
    gImbProjF_ = 0;
    gImbProjH_ = 0;
    gImbProjL_ = 0;    
    gImbPerpF_ = 0;
    gImbPerpH_ = 0;
    gImbPerpL_ = 0;

    gImbProj5_1_ = 0;
    gImbProj1_2_ = 0;
    gImbProj2_4_ = 0;
    gImbProj4_8_ = 0;
    gImbProj8_100_ = 0;

    gRPFImbProjF_ = 0;
    gRPFImbProjH_ = 0;
    gRPFImbProjL_ = 0;
    gRPFImbPerpF_ = 0;
    gRPFImbPerpH_ = 0;
    gRPFImbPerpL_ = 0;

    gRPFImbProj5_1_ = 0;
    gRPFImbProj1_2_ = 0;
    gRPFImbProj2_4_ = 0;
    gRPFImbProj4_8_ = 0;
    gRPFImbProj8_100_ = 0;

    gRPFImbProjFCorr_ = 0;
    gRPFImbProjHCorr_ = 0;
    gRPFImbProjLCorr_ = 0;
    gRPFImbPerpFCorr_ = 0;
    gRPFImbPerpHCorr_ = 0;
    gRPFImbPerpLCorr_ = 0;

    gRPFImbProjCorr5_1_ = 0;
    gRPFImbProjCorr1_2_ = 0;
    gRPFImbProjCorr2_4_ = 0;
    gRPFImbProjCorr4_8_ = 0;
    gRPFImbProjCorr8_100_ = 0;
  }
}

#endif
