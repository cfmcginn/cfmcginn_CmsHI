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
Float_t trkPtPF_[MAXTRKS];
Float_t trkPtCalo_[MAXTRKS];
Float_t trkPtT_[MAXTRKS];
Float_t trkPhi_[MAXTRKS];
Float_t trkEta_[MAXTRKS];
Float_t trkPFLeadDelPhi_[MAXTRKS];
Float_t trkCaloLeadDelPhi_[MAXTRKS];

Float_t trkRMinPF_[MAXTRKS];
Float_t trkPtCorrPF_[MAXTRKS];
Float_t trkPtFactPF_[MAXTRKS];
Float_t trkRMinCalo_[MAXTRKS];
Float_t trkPtCorrCalo_[MAXTRKS];
Float_t trkPtFactCalo_[MAXTRKS];
Float_t trkRMinT_[MAXTRKS];
Float_t trkPtCorrT_[MAXTRKS];
Float_t trkPtFactT_[MAXTRKS];


Float_t rPFImbProjF_;
Float_t rPFImbProjH_;
Float_t rPFImbProjL_;
Float_t rPFImbPerpF_;
Float_t rPFImbPerpH_;
Float_t rPFImbPerpL_;

Float_t rPFImbProj0_1_;
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

Float_t rPFImbProj0_1Corr_;
Float_t rPFImbProj1_2Corr_;
Float_t rPFImbProj2_4Corr_;
Float_t rPFImbProj4_8Corr_;
Float_t rPFImbProj8_100Corr_;

Float_t rCaloImbProjF_;
Float_t rCaloImbProjH_;
Float_t rCaloImbProjL_;
Float_t rCaloImbPerpF_;
Float_t rCaloImbPerpH_;
Float_t rCaloImbPerpL_;

Float_t rCaloImbProj0_1_;
Float_t rCaloImbProj1_2_;
Float_t rCaloImbProj2_4_;
Float_t rCaloImbProj4_8_;
Float_t rCaloImbProj8_100_;

Float_t rCaloImbProjFCorr_;
Float_t rCaloImbProjHCorr_;
Float_t rCaloImbProjLCorr_;
Float_t rCaloImbPerpFCorr_;
Float_t rCaloImbPerpHCorr_;
Float_t rCaloImbPerpLCorr_;

Float_t rCaloImbProj0_1Corr_;
Float_t rCaloImbProj1_2Corr_;
Float_t rCaloImbProj2_4Corr_;
Float_t rCaloImbProj4_8Corr_;
Float_t rCaloImbProj8_100Corr_;

Float_t rTImbPerpF_;
Float_t rTImbPerpH_;
Float_t rTImbPerpL_;
Float_t rTImbProjF_;
Float_t rTImbProjH_;
Float_t rTImbProjL_;

Float_t rTImbProj0_1_;
Float_t rTImbProj1_2_;
Float_t rTImbProj2_4_;
Float_t rTImbProj4_8_;
Float_t rTImbProj8_100_;

Float_t rTImbPerpFCorr_;
Float_t rTImbPerpHCorr_;
Float_t rTImbPerpLCorr_;
Float_t rTImbProjFCorr_;
Float_t rTImbProjHCorr_;
Float_t rTImbProjLCorr_;

Float_t rTImbProj0_1Corr_;
Float_t rTImbProj1_2Corr_;
Float_t rTImbProj2_4Corr_;
Float_t rTImbProj4_8Corr_;
Float_t rTImbProj8_100Corr_;

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

Float_t TLeadJtPt_;
Float_t TLeadJtPhi_;
Float_t TLeadJtEta_;
Float_t TSubLeadJtPt_;
Float_t TSubLeadJtPhi_;
Float_t TSubLeadJtEta_;
Float_t TJtAsymm_;

Float_t PFLeadJtPt_;
Float_t PFLeadJtPhi_;
Float_t PFLeadJtEta_;
Float_t PFSubLeadJtPt_;
Float_t PFSubLeadJtPhi_;
Float_t PFSubLeadJtEta_;
Float_t PFJtAsymm_;

Float_t CaloLeadJtPt_;
Float_t CaloLeadJtPhi_;
Float_t CaloLeadJtEta_;
Float_t CaloSubLeadJtPt_;
Float_t CaloSubLeadJtPhi_;
Float_t CaloSubLeadJtEta_;
Float_t CaloJtAsymm_;

//Gen Tree Variables

const int MAXGEN = 50000; //From SetupGenParticleTree.h

Int_t nGen_;
Float_t genPt_[MAXGEN];
Float_t genPtPF_[MAXGEN];
Float_t genPtCalo_[MAXGEN];
Float_t genPtT_[MAXGEN];
Float_t genPhi_[MAXGEN];
Float_t genEta_[MAXGEN];
Float_t genLeadDelPhi_[MAXGEN];

Float_t gTImbProjF_;
Float_t gTImbProjH_;
Float_t gTImbProjL_;
Float_t gTImbPerpF_;
Float_t gTImbPerpH_;
Float_t gTImbPerpL_;

Float_t gTImbProj0_1_;
Float_t gTImbProj1_2_;
Float_t gTImbProj2_4_;
Float_t gTImbProj4_8_;
Float_t gTImbProj8_100_;

Float_t gPFImbProjF_;
Float_t gPFImbProjH_;
Float_t gPFImbProjL_;
Float_t gPFImbPerpF_;
Float_t gPFImbPerpH_;
Float_t gPFImbPerpL_;

Float_t gPFImbProj0_1_;
Float_t gPFImbProj1_2_;
Float_t gPFImbProj2_4_;
Float_t gPFImbProj4_8_;
Float_t gPFImbProj8_100_;

Float_t gCaloImbProjF_;
Float_t gCaloImbProjH_;
Float_t gCaloImbProjL_;
Float_t gCaloImbPerpF_;
Float_t gCaloImbPerpH_;
Float_t gCaloImbPerpL_;

Float_t gCaloImbProj0_1_;
Float_t gCaloImbProj1_2_;
Float_t gCaloImbProj2_4_;
Float_t gCaloImbProj4_8_;
Float_t gCaloImbProj8_100_;



void SetBranches(bool montecarlo)
{
  //Track Tree Branches
  
  trackTree_p->Branch("nTrk", &nTrk_, "nTrk/I");
  trackTree_p->Branch("trkPt", &trkPt_, "trkPt[nTrk]/F");
  trackTree_p->Branch("trkPtPF", &trkPtPF_, "trkPtPF[nTrk]/F");
  trackTree_p->Branch("trkPtCalo", &trkPtCalo_, "trkPtCalo[nTrk]/F");
  trackTree_p->Branch("trkPtT", &trkPtT_, "trkPtT[nTrk]/F");
  trackTree_p->Branch("trkPhi", &trkPhi_, "trkPhi[nTrk]/F");
  trackTree_p->Branch("trkEta", &trkEta_, "trkEta[nTrk]/F");
  trackTree_p->Branch("trkPFLeadDelPhi", &trkPFLeadDelPhi_, "trkPFLeadDelPhi[nTrk]/F");
  trackTree_p->Branch("trkCaloLeadDelPhi", &trkCaloLeadDelPhi_, "trkCaloLeadDelPhi[nTrk]/F");

  trackTree_p->Branch("trkRMinPF", &trkRMinPF_, "trkRMinPF[nTrk]/F");
  trackTree_p->Branch("trkPtCorrPF", &trkPtCorrPF_, "trkPtCorrPF[nTrk]/F");
  trackTree_p->Branch("trkPtFactPF", &trkPtFactPF_, "trkPtFactPF[nTrk]/F");
  trackTree_p->Branch("trkRMinCalo", &trkRMinCalo_, "trkRMinCalo[nTrk]/F");
  trackTree_p->Branch("trkPtCorrCalo", &trkPtCorrCalo_, "trkPtCorrCalo[nTrk]/F");
  trackTree_p->Branch("trkPtFactCalo", &trkPtFactCalo_, "trkPtFactCalo[nTrk]/F");
  trackTree_p->Branch("trkRMinT", &trkRMinT_, "trkRMinT[nTrk]/F");
  trackTree_p->Branch("trkPtCorrT", &trkPtCorrT_, "trkPtCorrT[nTrk]/F");
  trackTree_p->Branch("trkPtFactT", &trkPtFactT_, "trkPtFactT[nTrk]/F");
  
  trackTree_p->Branch("rPFImbProjF", &rPFImbProjF_, "rPFImbProjF/F");
  trackTree_p->Branch("rPFImbProjH", &rPFImbProjH_, "rPFImbProjH/F");
  trackTree_p->Branch("rPFImbProjL", &rPFImbProjL_, "rPFImbProjL/F");
  trackTree_p->Branch("rPFImbPerpF", &rPFImbPerpF_, "rPFImbPerpF/F");
  trackTree_p->Branch("rPFImbPerpH", &rPFImbPerpH_, "rPFImbPerpH/F");
  trackTree_p->Branch("rPFImbPerpL", &rPFImbPerpL_, "rPFImbPerpL/F");

  trackTree_p->Branch("rPFImbProj0_1", &rPFImbProj0_1_, "rPFImbProj0_1/F");
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

  trackTree_p->Branch("rPFImbProj0_1Corr", &rPFImbProj0_1Corr_, "rPFImbProj0_1Corr/F");
  trackTree_p->Branch("rPFImbProj1_2Corr", &rPFImbProj1_2Corr_, "rPFImbProj1_2Corr/F");
  trackTree_p->Branch("rPFImbProj2_4Corr", &rPFImbProj2_4Corr_, "rPFImbProj2_4Corr/F");
  trackTree_p->Branch("rPFImbProj4_8Corr", &rPFImbProj4_8Corr_, "rPFImbProj4_8Corr/F");
  trackTree_p->Branch("rPFImbProj8_100Corr", &rPFImbProj8_100Corr_, "rPFImbProj8_100Corr/F");

  trackTree_p->Branch("rCaloImbProjF", &rCaloImbProjF_, "rCaloImbProjF/F");
  trackTree_p->Branch("rCaloImbProjH", &rCaloImbProjH_, "rCaloImbProjH/F");
  trackTree_p->Branch("rCaloImbProjL", &rCaloImbProjL_, "rCaloImbProjL/F");
  trackTree_p->Branch("rCaloImbPerpF", &rCaloImbPerpF_, "rCaloImbPerpF/F");
  trackTree_p->Branch("rCaloImbPerpH", &rCaloImbPerpH_, "rCaloImbPerpH/F");
  trackTree_p->Branch("rCaloImbPerpL", &rCaloImbPerpL_, "rCaloImbPerpL/F");

  trackTree_p->Branch("rCaloImbProj0_1", &rCaloImbProj0_1_, "rCaloImbProj0_1/F");
  trackTree_p->Branch("rCaloImbProj1_2", &rCaloImbProj1_2_, "rCaloImbProj1_2/F");
  trackTree_p->Branch("rCaloImbProj2_4", &rCaloImbProj2_4_, "rCaloImbProj2_4/F");
  trackTree_p->Branch("rCaloImbProj4_8", &rCaloImbProj4_8_, "rCaloImbProj4_8/F");
  trackTree_p->Branch("rCaloImbProj8_100", &rCaloImbProj8_100_, "rCaloImbProj8_100/F");

  trackTree_p->Branch("rCaloImbProjFCorr", &rCaloImbProjFCorr_, "rCaloImbProjFCorr/F");
  trackTree_p->Branch("rCaloImbProjHCorr", &rCaloImbProjHCorr_, "rCaloImbProjHCorr/F");
  trackTree_p->Branch("rCaloImbProjLCorr", &rCaloImbProjLCorr_, "rCaloImbProjLCorr/F");
  trackTree_p->Branch("rCaloImbPerpFCorr", &rCaloImbPerpFCorr_, "rCaloImbPerpFCorr/F");
  trackTree_p->Branch("rCaloImbPerpHCorr", &rCaloImbPerpHCorr_, "rCaloImbPerpHCorr/F");
  trackTree_p->Branch("rCaloImbPerpLCorr", &rCaloImbPerpLCorr_, "rCaloImbPerpLCorr/F");

  trackTree_p->Branch("rCaloImbProj0_1Corr", &rCaloImbProj0_1Corr_, "rCaloImbProj0_1Corr/F");
  trackTree_p->Branch("rCaloImbProj1_2Corr", &rCaloImbProj1_2Corr_, "rCaloImbProj1_2Corr/F");
  trackTree_p->Branch("rCaloImbProj2_4Corr", &rCaloImbProj2_4Corr_, "rCaloImbProj2_4Corr/F");
  trackTree_p->Branch("rCaloImbProj4_8Corr", &rCaloImbProj4_8Corr_, "rCaloImbProj4_8Corr/F");
  trackTree_p->Branch("rCaloImbProj8_100Corr", &rCaloImbProj8_100Corr_, "rCaloImbProj8_100Corr/F");

  if(montecarlo){
    //Track tree branches iff truth avail.
    trackTree_p->Branch("rTImbProjF", &rTImbProjF_, "rTImbProjF/F");
    trackTree_p->Branch("rTImbProjH", &rTImbProjH_, "rTImbProjH/F");
    trackTree_p->Branch("rTImbProjL", &rTImbProjL_, "rTImbProjL/F");    
    trackTree_p->Branch("rTImbPerpF", &rTImbPerpF_, "rTImbPerpF/F");
    trackTree_p->Branch("rTImbPerpH", &rTImbPerpH_, "rTImbPerpH/F");
    trackTree_p->Branch("rTImbPerpL", &rTImbPerpL_, "rTImbPerpL/F");

    trackTree_p->Branch("rTImbProj0_1", &rTImbProj0_1_, "rTImbProj0_1/F");
    trackTree_p->Branch("rTImbProj1_2", &rTImbProj1_2_, "rTImbProj1_2/F");
    trackTree_p->Branch("rTImbProj2_4", &rTImbProj2_4_, "rTImbProj2_4/F");
    trackTree_p->Branch("rTImbProj4_8", &rTImbProj4_8_, "rTImbProj4_8/F");
    trackTree_p->Branch("rTImbProj8_100", &rTImbProj8_100_, "rTImbProj8_100/F");

    trackTree_p->Branch("rTImbProjFCorr", &rTImbProjFCorr_, "rTImbProjFCorr/F");
    trackTree_p->Branch("rTImbProjHCorr", &rTImbProjHCorr_, "rTImbProjHCorr/F");
    trackTree_p->Branch("rTImbProjLCorr", &rTImbProjLCorr_, "rTImbProjLCorr/F");    
    trackTree_p->Branch("rTImbPerpFCorr", &rTImbPerpFCorr_, "rTImbPerpFCorr/F");
    trackTree_p->Branch("rTImbPerpHCorr", &rTImbPerpHCorr_, "rTImbPerpHCorr/F");
    trackTree_p->Branch("rTImbPerpLCorr", &rTImbPerpLCorr_, "rTImbPerpLCorr/F");

    trackTree_p->Branch("rTImbProj0_1Corr", &rTImbProj0_1Corr_, "rTImbProj0_1Corr/F");
    trackTree_p->Branch("rTImbProj1_2Corr", &rTImbProj1_2Corr_, "rTImbProj1_2Corr/F");
    trackTree_p->Branch("rTImbProj2_4Corr", &rTImbProj2_4Corr_, "rTImbProj2_4Corr/F");
    trackTree_p->Branch("rTImbProj4_8Corr", &rTImbProj4_8Corr_, "rTImbProj4_8Corr/F");
    trackTree_p->Branch("rTImbProj8_100Corr", &rTImbProj8_100Corr_, "rTImbProj8_100Corr/F");
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

    jetTree_p->Branch("TLeadJtPt", &TLeadJtPt_, "TLeadJtPt/F");
    jetTree_p->Branch("TLeadJtPhi", &TLeadJtPhi_, "TLeadJtPhi/F");
    jetTree_p->Branch("TLeadJtEta", &TLeadJtEta_, "TLeadJtEta/F");
    jetTree_p->Branch("TSubLeadJtPt", &TSubLeadJtPt_, "TSubLeadJtPt/F");
    jetTree_p->Branch("TSubLeadJtPhi", &TSubLeadJtPhi_, "TSubLeadJtPhi/F");
    jetTree_p->Branch("TSubLeadJtEta", &TSubLeadJtEta_, "TSubLeadJtEta/F");
    jetTree_p->Branch("TJtAsymm", &TJtAsymm_, "TJtAsymm/F");
  }

  jetTree_p->Branch("PFLeadJtPt", &PFLeadJtPt_, "PFLeadJtPt/F");
  jetTree_p->Branch("PFLeadJtPhi", &PFLeadJtPhi_, "PFLeadJtPhi/F");
  jetTree_p->Branch("PFLeadJtEta", &PFLeadJtEta_, "PFLeadJtEta/F");
  jetTree_p->Branch("PFSubLeadJtPt", &PFSubLeadJtPt_, "PFSubLeadJtPt/F");
  jetTree_p->Branch("PFSubLeadJtPhi", &PFSubLeadJtPhi_, "PFSubLeadJtPhi/F");
  jetTree_p->Branch("PFSubLeadJtEta", &PFSubLeadJtEta_, "PFSubLeadJtEta/F");
  jetTree_p->Branch("PFJtAsymm", &PFJtAsymm_, "PFJtAsymm/F");
    
  jetTree_p->Branch("CaloLeadJtPt", &CaloLeadJtPt_, "CaloLeadJtPt/F");
  jetTree_p->Branch("CaloLeadJtPhi", &CaloLeadJtPhi_, "CaloLeadJtPhi/F");
  jetTree_p->Branch("CaloLeadJtEta", &CaloLeadJtEta_, "CaloLeadJtEta/F");
  jetTree_p->Branch("CaloSubLeadJtPt", &CaloSubLeadJtPt_, "CaloSubLeadJtPt/F");
  jetTree_p->Branch("CaloSubLeadJtPhi", &CaloSubLeadJtPhi_, "CaloSubLeadJtPhi/F");
  jetTree_p->Branch("CaloSubLeadJtEta", &CaloSubLeadJtEta_, "CaloSubLeadJtEta/F");
  jetTree_p->Branch("CaloJtAsymm", &CaloJtAsymm_, "CaloJtAsymm/F");

  if(montecarlo){
    //Gen Tree Branches

    genTree_p->Branch("nGen", &nGen_, "nGen/I");
    genTree_p->Branch("genPt", &genPt_, "genPt[nGen]/F");
    genTree_p->Branch("genPtPF", &genPtPF_, "genPtPF[nGen]/F");
    genTree_p->Branch("genPtCalo", &genPtCalo_, "genPtCalo[nGen]/F");
    genTree_p->Branch("genPtT", &genPtT_, "genPtT[nGen]/F");
    genTree_p->Branch("genPhi", &genPhi_, "genPhi[nGen]/F");
    genTree_p->Branch("genEta", &genEta_, "genEta[nGen]/F");
    genTree_p->Branch("genLeadDelPhi", &genLeadDelPhi_, "genLeadDelPhi[nGen]/F");

    genTree_p->Branch("gTImbProjF", &gTImbProjF_, "gTImbProjF/F");
    genTree_p->Branch("gTImbProjH", &gTImbProjH_, "gTImbProjH/F");
    genTree_p->Branch("gTImbProjL", &gTImbProjL_, "gTImbProjL/F");
    genTree_p->Branch("gTImbPerpF", &gTImbPerpF_, "gTImbPerpF/F");
    genTree_p->Branch("gTImbPerpH", &gTImbPerpH_, "gTImbPerpH/F");
    genTree_p->Branch("gTImbPerpL", &gTImbPerpL_, "gTImbPerpL/F");

    genTree_p->Branch("gTImbProj0_1", &gTImbProj0_1_, "gTImbProj0_1/F");
    genTree_p->Branch("gTImbProj1_2", &gTImbProj1_2_, "gTImbProj1_2/F");
    genTree_p->Branch("gTImbProj2_4", &gTImbProj2_4_, "gTImbProj2_4/F");
    genTree_p->Branch("gTImbProj4_8", &gTImbProj4_8_, "gTImbProj4_8/F");
    genTree_p->Branch("gTImbProj8_100", &gTImbProj8_100_, "gTImbProj8_100/F");

    genTree_p->Branch("gPFImbProjF", &gPFImbProjF_, "gPFImbProjF/F");
    genTree_p->Branch("gPFImbProjH", &gPFImbProjH_, "gPFImbProjH/F");
    genTree_p->Branch("gPFImbProjL", &gPFImbProjL_, "gPFImbProjL/F");
    genTree_p->Branch("gPFImbPerpF", &gPFImbPerpF_, "gPFImbPerpF/F");
    genTree_p->Branch("gPFImbPerpH", &gPFImbPerpH_, "gPFImbPerpH/F");
    genTree_p->Branch("gPFImbPerpL", &gPFImbPerpL_, "gPFImbPerpL/F");

    genTree_p->Branch("gPFImbProj0_1", &gPFImbProj0_1_, "gPFImbProj0_1/F");
    genTree_p->Branch("gPFImbProj1_2", &gPFImbProj1_2_, "gPFImbProj1_2/F");
    genTree_p->Branch("gPFImbProj2_4", &gPFImbProj2_4_, "gPFImbProj2_4/F");
    genTree_p->Branch("gPFImbProj4_8", &gPFImbProj4_8_, "gPFImbProj4_8/F");
    genTree_p->Branch("gPFImbProj8_100", &gPFImbProj8_100_, "gPFImbProj8_100/F");

    genTree_p->Branch("gCaloImbProjF", &gCaloImbProjF_, "gCaloImbProjF/F");
    genTree_p->Branch("gCaloImbProjH", &gCaloImbProjH_, "gCaloImbProjH/F");
    genTree_p->Branch("gCaloImbProjL", &gCaloImbProjL_, "gCaloImbProjL/F");
    genTree_p->Branch("gCaloImbPerpF", &gCaloImbPerpF_, "gCaloImbPerpF/F");
    genTree_p->Branch("gCaloImbPerpH", &gCaloImbPerpH_, "gCaloImbPerpH/F");
    genTree_p->Branch("gCaloImbPerpL", &gCaloImbPerpL_, "gCaloImbPerpL/F");

    genTree_p->Branch("gCaloImbProj0_1", &gCaloImbProj0_1_, "gCaloImbProj0_1/F");
    genTree_p->Branch("gCaloImbProj1_2", &gCaloImbProj1_2_, "gCaloImbProj1_2/F");
    genTree_p->Branch("gCaloImbProj2_4", &gCaloImbProj2_4_, "gCaloImbProj2_4/F");
    genTree_p->Branch("gCaloImbProj4_8", &gCaloImbProj4_8_, "gCaloImbProj4_8/F");
    genTree_p->Branch("gCaloImbProj8_100", &gCaloImbProj8_100_, "gCaloImbProj8_100/F");
  }
}


void GetBranches(bool montecarlo)
{
  //Track Tree Branches

  trackTree_p->SetBranchAddress("nTrk", &nTrk_);
  trackTree_p->SetBranchAddress("trkPt", &trkPt_);
  trackTree_p->SetBranchAddress("trkPtPF", &trkPtPF_);
  trackTree_p->SetBranchAddress("trkPtCalo", &trkPtCalo_);
  trackTree_p->SetBranchAddress("trkPtT", &trkPtT_);
  trackTree_p->SetBranchAddress("trkPhi", &trkPhi_);
  trackTree_p->SetBranchAddress("trkEta", &trkEta_);
  trackTree_p->SetBranchAddress("trkPFLeadDelPhi", &trkPFLeadDelPhi_);
  trackTree_p->SetBranchAddress("trkCaloLeadDelPhi", &trkCaloLeadDelPhi_);

  trackTree_p->SetBranchAddress("trkRMinPF", &trkRMinPF_);
  trackTree_p->SetBranchAddress("trkPtCorrPF", &trkPtCorrPF_);
  trackTree_p->SetBranchAddress("trkPtFactPF", &trkPtFactPF_);
  trackTree_p->SetBranchAddress("trkRMinCalo", &trkRMinCalo_);
  trackTree_p->SetBranchAddress("trkPtCorrCalo", &trkPtCorrCalo_);
  trackTree_p->SetBranchAddress("trkPtFactCalo", &trkPtFactCalo_);
  trackTree_p->SetBranchAddress("trkRMinT", &trkRMinT_);
  trackTree_p->SetBranchAddress("trkPtCorrT", &trkPtCorrT_);
  trackTree_p->SetBranchAddress("trkPtFactT", &trkPtFactT_);

  trackTree_p->SetBranchAddress("rPFImbProjF", &rPFImbProjF_);
  trackTree_p->SetBranchAddress("rPFImbProjH", &rPFImbProjH_);
  trackTree_p->SetBranchAddress("rPFImbProjL", &rPFImbProjL_);
  trackTree_p->SetBranchAddress("rPFImbPerpF", &rPFImbPerpF_);
  trackTree_p->SetBranchAddress("rPFImbPerpH", &rPFImbPerpH_);
  trackTree_p->SetBranchAddress("rPFImbPerpL", &rPFImbPerpL_);

  trackTree_p->SetBranchAddress("rPFImbProj0_1", &rPFImbProj0_1_);
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

  trackTree_p->SetBranchAddress("rPFImbProj0_1Corr", &rPFImbProj0_1Corr_);
  trackTree_p->SetBranchAddress("rPFImbProj1_2Corr", &rPFImbProj1_2Corr_);
  trackTree_p->SetBranchAddress("rPFImbProj2_4Corr", &rPFImbProj2_4Corr_);
  trackTree_p->SetBranchAddress("rPFImbProj4_8Corr", &rPFImbProj4_8Corr_);
  trackTree_p->SetBranchAddress("rPFImbProj8_100Corr", &rPFImbProj8_100Corr_);

  trackTree_p->SetBranchAddress("rCaloImbProjF", &rCaloImbProjF_);
  trackTree_p->SetBranchAddress("rCaloImbProjH", &rCaloImbProjH_);
  trackTree_p->SetBranchAddress("rCaloImbProjL", &rCaloImbProjL_);
  trackTree_p->SetBranchAddress("rCaloImbPerpF", &rCaloImbPerpF_);
  trackTree_p->SetBranchAddress("rCaloImbPerpH", &rCaloImbPerpH_);
  trackTree_p->SetBranchAddress("rCaloImbPerpL", &rCaloImbPerpL_);

  trackTree_p->SetBranchAddress("rCaloImbProj0_1", &rCaloImbProj0_1_);
  trackTree_p->SetBranchAddress("rCaloImbProj1_2", &rCaloImbProj1_2_);
  trackTree_p->SetBranchAddress("rCaloImbProj2_4", &rCaloImbProj2_4_);
  trackTree_p->SetBranchAddress("rCaloImbProj4_8", &rCaloImbProj4_8_);
  trackTree_p->SetBranchAddress("rCaloImbProj8_100", &rCaloImbProj8_100_);

  trackTree_p->SetBranchAddress("rCaloImbProjFCorr", &rCaloImbProjFCorr_);
  trackTree_p->SetBranchAddress("rCaloImbProjHCorr", &rCaloImbProjHCorr_);
  trackTree_p->SetBranchAddress("rCaloImbProjLCorr", &rCaloImbProjLCorr_);
  trackTree_p->SetBranchAddress("rCaloImbPerpFCorr", &rCaloImbPerpFCorr_);
  trackTree_p->SetBranchAddress("rCaloImbPerpHCorr", &rCaloImbPerpHCorr_);
  trackTree_p->SetBranchAddress("rCaloImbPerpLCorr", &rCaloImbPerpLCorr_);

  trackTree_p->SetBranchAddress("rCaloImbProj0_1Corr", &rCaloImbProj0_1Corr_);
  trackTree_p->SetBranchAddress("rCaloImbProj1_2Corr", &rCaloImbProj1_2Corr_);
  trackTree_p->SetBranchAddress("rCaloImbProj2_4Corr", &rCaloImbProj2_4Corr_);
  trackTree_p->SetBranchAddress("rCaloImbProj4_8Corr", &rCaloImbProj4_8Corr_);
  trackTree_p->SetBranchAddress("rCaloImbProj8_100Corr", &rCaloImbProj8_100Corr_);

  if(montecarlo){
    //Track Tree Branches iff. Truth avail.
    trackTree_p->SetBranchAddress("rTImbProjF", &rTImbProjF_);
    trackTree_p->SetBranchAddress("rTImbProjH", &rTImbProjH_);
    trackTree_p->SetBranchAddress("rTImbProjL", &rTImbProjL_);    
    trackTree_p->SetBranchAddress("rTImbPerpF", &rTImbPerpF_);
    trackTree_p->SetBranchAddress("rTImbPerpH", &rTImbPerpH_);
    trackTree_p->SetBranchAddress("rTImbPerpL", &rTImbPerpL_);

    trackTree_p->SetBranchAddress("rTImbProj0_1", &rTImbProj0_1_);
    trackTree_p->SetBranchAddress("rTImbProj1_2", &rTImbProj1_2_);
    trackTree_p->SetBranchAddress("rTImbProj2_4", &rTImbProj2_4_);
    trackTree_p->SetBranchAddress("rTImbProj4_8", &rTImbProj4_8_);
    trackTree_p->SetBranchAddress("rTImbProj8_100", &rTImbProj8_100_);

    trackTree_p->SetBranchAddress("rTImbProjFCorr", &rTImbProjFCorr_);
    trackTree_p->SetBranchAddress("rTImbProjHCorr", &rTImbProjHCorr_);
    trackTree_p->SetBranchAddress("rTImbProjLCorr", &rTImbProjLCorr_);    
    trackTree_p->SetBranchAddress("rTImbPerpFCorr", &rTImbPerpFCorr_);
    trackTree_p->SetBranchAddress("rTImbPerpHCorr", &rTImbPerpHCorr_);
    trackTree_p->SetBranchAddress("rTImbPerpLCorr", &rTImbPerpLCorr_);

    trackTree_p->SetBranchAddress("rTImbProj0_1Corr", &rTImbProj0_1Corr_);
    trackTree_p->SetBranchAddress("rTImbProj1_2Corr", &rTImbProj1_2Corr_);
    trackTree_p->SetBranchAddress("rTImbProj2_4Corr", &rTImbProj2_4Corr_);
    trackTree_p->SetBranchAddress("rTImbProj4_8Corr", &rTImbProj4_8Corr_);
    trackTree_p->SetBranchAddress("rTImbProj8_100Corr", &rTImbProj8_100Corr_);
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

    jetTree_p->SetBranchAddress("TLeadJtPt", &TLeadJtPt_);
    jetTree_p->SetBranchAddress("TLeadJtPhi", &TLeadJtPhi_);
    jetTree_p->SetBranchAddress("TLeadJtEta", &TLeadJtEta_);
    jetTree_p->SetBranchAddress("TSubLeadJtPt", &TSubLeadJtPt_);
    jetTree_p->SetBranchAddress("TSubLeadJtPhi", &TSubLeadJtPhi_);
    jetTree_p->SetBranchAddress("TSubLeadJtEta", &TSubLeadJtEta_);
    jetTree_p->SetBranchAddress("TJtAsymm", &TJtAsymm_);
  }

  jetTree_p->SetBranchAddress("PFLeadJtPt", &PFLeadJtPt_);
  jetTree_p->SetBranchAddress("PFLeadJtPhi", &PFLeadJtPhi_);
  jetTree_p->SetBranchAddress("PFLeadJtEta", &PFLeadJtEta_);
  jetTree_p->SetBranchAddress("PFSubLeadJtPt", &PFSubLeadJtPt_);
  jetTree_p->SetBranchAddress("PFSubLeadJtPhi", &PFSubLeadJtPhi_);
  jetTree_p->SetBranchAddress("PFSubLeadJtEta", &PFSubLeadJtEta_);
  jetTree_p->SetBranchAddress("PFJtAsymm", &PFJtAsymm_);

  jetTree_p->SetBranchAddress("CaloLeadJtPt", &CaloLeadJtPt_);
  jetTree_p->SetBranchAddress("CaloLeadJtPhi", &CaloLeadJtPhi_);
  jetTree_p->SetBranchAddress("CaloLeadJtEta", &CaloLeadJtEta_);
  jetTree_p->SetBranchAddress("CaloSubLeadJtPt", &CaloSubLeadJtPt_);
  jetTree_p->SetBranchAddress("CaloSubLeadJtPhi", &CaloSubLeadJtPhi_);
  jetTree_p->SetBranchAddress("CaloSubLeadJtEta", &CaloSubLeadJtEta_);
  jetTree_p->SetBranchAddress("CaloJtAsymm", &CaloJtAsymm_);

  if(montecarlo){
    //Gen Tree Branches

    genTree_p->SetBranchAddress("nGen", &nGen_);
    genTree_p->SetBranchAddress("genPt", &genPt_);
    genTree_p->SetBranchAddress("genPtPF", &genPtPF_);
    genTree_p->SetBranchAddress("genPtCalo", &genPtCalo_);
    genTree_p->SetBranchAddress("genPtT", &genPtT_);
    genTree_p->SetBranchAddress("genPhi", &genPhi_);
    genTree_p->SetBranchAddress("genEta", &genEta_);
    genTree_p->SetBranchAddress("genLeadDelPhi", &genLeadDelPhi_);

    genTree_p->SetBranchAddress("gTImbProjF", &gTImbProjF_);
    genTree_p->SetBranchAddress("gTImbProjH", &gTImbProjH_);
    genTree_p->SetBranchAddress("gTImbProjL", &gTImbProjL_);
    genTree_p->SetBranchAddress("gTImbPerpF", &gTImbPerpF_);
    genTree_p->SetBranchAddress("gTImbPerpH", &gTImbPerpH_);
    genTree_p->SetBranchAddress("gTImbPerpL", &gTImbPerpL_);

    genTree_p->SetBranchAddress("gTImbProj0_1", &gTImbProj0_1_);
    genTree_p->SetBranchAddress("gTImbProj1_2", &gTImbProj1_2_);
    genTree_p->SetBranchAddress("gTImbProj2_4", &gTImbProj2_4_);
    genTree_p->SetBranchAddress("gTImbProj4_8", &gTImbProj4_8_);
    genTree_p->SetBranchAddress("gTImbProj8_100", &gTImbProj8_100_);

    genTree_p->SetBranchAddress("gPFImbProjF", &gPFImbProjF_);
    genTree_p->SetBranchAddress("gPFImbProjH", &gPFImbProjH_);
    genTree_p->SetBranchAddress("gPFImbProjL", &gPFImbProjL_);
    genTree_p->SetBranchAddress("gPFImbPerpF", &gPFImbPerpF_);
    genTree_p->SetBranchAddress("gPFImbPerpH", &gPFImbPerpH_);
    genTree_p->SetBranchAddress("gPFImbPerpL", &gPFImbPerpL_);

    genTree_p->SetBranchAddress("gPFImbProj0_1", &gPFImbProj0_1_);
    genTree_p->SetBranchAddress("gPFImbProj1_2", &gPFImbProj1_2_);
    genTree_p->SetBranchAddress("gPFImbProj2_4", &gPFImbProj2_4_);
    genTree_p->SetBranchAddress("gPFImbProj4_8", &gPFImbProj4_8_);
    genTree_p->SetBranchAddress("gPFImbProj8_100", &gPFImbProj8_100_);

    genTree_p->SetBranchAddress("gCaloImbProjF", &gCaloImbProjF_);
    genTree_p->SetBranchAddress("gCaloImbProjH", &gCaloImbProjH_);
    genTree_p->SetBranchAddress("gCaloImbProjL", &gCaloImbProjL_);
    genTree_p->SetBranchAddress("gCaloImbPerpF", &gCaloImbPerpF_);
    genTree_p->SetBranchAddress("gCaloImbPerpH", &gCaloImbPerpH_);
    genTree_p->SetBranchAddress("gCaloImbPerpL", &gCaloImbPerpL_);

    genTree_p->SetBranchAddress("gCaloImbProj0_1", &gCaloImbProj0_1_);
    genTree_p->SetBranchAddress("gCaloImbProj1_2", &gCaloImbProj1_2_);
    genTree_p->SetBranchAddress("gCaloImbProj2_4", &gCaloImbProj2_4_);
    genTree_p->SetBranchAddress("gCaloImbProj4_8", &gCaloImbProj4_8_);
    genTree_p->SetBranchAddress("gCaloImbProj8_100", &gCaloImbProj8_100_);
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

    TLeadJtPt_ = -10;
    TSubLeadJtPt_ = -10;
    TLeadJtPhi_ = -10;
    TSubLeadJtPhi_ = -10;
    TLeadJtEta_ = -10;
    TSubLeadJtEta_ = -10;
    TJtAsymm_ = -10;
  }

  recoPFSet_ = false;
  recoCaloSet_ = false;

  PFLeadJtPt_ = -10;
  PFSubLeadJtPt_ = -10;
  PFLeadJtPhi_ = -10;
  PFSubLeadJtPhi_ = -10;
  PFLeadJtEta_ = -10;
  PFSubLeadJtEta_ = -10;
  PFJtAsymm_ = -10;

  CaloLeadJtPt_ = -10;
  CaloSubLeadJtPt_ = -10;
  CaloLeadJtPhi_ = -10;
  CaloSubLeadJtPhi_ = -10;
  CaloLeadJtEta_ = -10;
  CaloSubLeadJtEta_ = -10;
  CaloJtAsymm_ = -10;
}

void InitProjPerp(bool montecarlo = false)
{
  rPFImbProjF_ = 0;
  rPFImbProjH_ = 0;
  rPFImbProjL_ = 0;
  rPFImbPerpF_ = 0;
  rPFImbPerpH_ = 0;
  rPFImbPerpL_ = 0;

  rPFImbProj0_1_ = 0;
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

  rPFImbProj0_1Corr_ = 0;
  rPFImbProj1_2Corr_ = 0;
  rPFImbProj2_4Corr_ = 0;
  rPFImbProj4_8Corr_ = 0;
  rPFImbProj8_100Corr_ = 0;

  rCaloImbProjF_ = 0;
  rCaloImbProjH_ = 0;
  rCaloImbProjL_ = 0;
  rCaloImbPerpF_ = 0;
  rCaloImbPerpH_ = 0;
  rCaloImbPerpL_ = 0;

  rCaloImbProj0_1_ = 0;
  rCaloImbProj1_2_ = 0;
  rCaloImbProj2_4_ = 0;
  rCaloImbProj4_8_ = 0;
  rCaloImbProj8_100_ = 0;

  rCaloImbProjFCorr_ = 0;
  rCaloImbProjHCorr_ = 0;
  rCaloImbProjLCorr_ = 0;
  rCaloImbPerpFCorr_ = 0;
  rCaloImbPerpHCorr_ = 0;
  rCaloImbPerpLCorr_ = 0;

  rCaloImbProj0_1Corr_ = 0;
  rCaloImbProj1_2Corr_ = 0;
  rCaloImbProj2_4Corr_ = 0;
  rCaloImbProj4_8Corr_ = 0;
  rCaloImbProj8_100Corr_ = 0;

  if(montecarlo){
    gTImbProjF_ = 0;
    gTImbProjH_ = 0;
    gTImbProjL_ = 0;    
    gTImbPerpF_ = 0;
    gTImbPerpH_ = 0;
    gTImbPerpL_ = 0;

    gTImbProj0_1_ = 0;
    gTImbProj1_2_ = 0;
    gTImbProj2_4_ = 0;
    gTImbProj4_8_ = 0;
    gTImbProj8_100_ = 0;

    gPFImbProjF_ = 0;
    gPFImbProjH_ = 0;
    gPFImbProjL_ = 0;    
    gPFImbPerpF_ = 0;
    gPFImbPerpH_ = 0;
    gPFImbPerpL_ = 0;

    gPFImbProj0_1_ = 0;
    gPFImbProj1_2_ = 0;
    gPFImbProj2_4_ = 0;
    gPFImbProj4_8_ = 0;
    gPFImbProj8_100_ = 0;

    gCaloImbProjF_ = 0;
    gCaloImbProjH_ = 0;
    gCaloImbProjL_ = 0;    
    gCaloImbPerpF_ = 0;
    gCaloImbPerpH_ = 0;
    gCaloImbPerpL_ = 0;

    gCaloImbProj0_1_ = 0;
    gCaloImbProj1_2_ = 0;
    gCaloImbProj2_4_ = 0;
    gCaloImbProj4_8_ = 0;
    gCaloImbProj8_100_ = 0;

    rTImbProjF_ = 0;
    rTImbProjH_ = 0;
    rTImbProjL_ = 0;
    rTImbPerpF_ = 0;
    rTImbPerpH_ = 0;
    rTImbPerpL_ = 0;

    rTImbProj0_1_ = 0;
    rTImbProj1_2_ = 0;
    rTImbProj2_4_ = 0;
    rTImbProj4_8_ = 0;
    rTImbProj8_100_ = 0;

    rTImbProjFCorr_ = 0;
    rTImbProjHCorr_ = 0;
    rTImbProjLCorr_ = 0;
    rTImbPerpFCorr_ = 0;
    rTImbPerpHCorr_ = 0;
    rTImbPerpLCorr_ = 0;

    rTImbProj0_1Corr_ = 0;
    rTImbProj1_2Corr_ = 0;
    rTImbProj2_4Corr_ = 0;
    rTImbProj4_8Corr_ = 0;
    rTImbProj8_100Corr_ = 0;
  }
}

#endif
