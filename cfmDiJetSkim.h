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

const int MAXTRKS = 12000; //From SetupTrackTree.h
Int_t nTrk_;
Float_t trkPt_[MAXTRKS];
Float_t trkPtPF_[MAXTRKS];
Float_t trkPtCalo_[MAXTRKS];
Float_t trkPtT_[MAXTRKS];
Float_t trkPtVsPF_[MAXTRKS];
Float_t trkPtVsCalo_[MAXTRKS];
Float_t trkPhi_[MAXTRKS];
Float_t trkEta_[MAXTRKS];
Float_t trkPFLeadDelPhi_[MAXTRKS];
Float_t trkCaloLeadDelPhi_[MAXTRKS];

Float_t trkRLeadPF_[MAXTRKS];
Float_t trkRMinPF_[MAXTRKS];
Float_t trkPtCorrPF_[MAXTRKS];
Float_t trkPtFactPF_[MAXTRKS];

Float_t trkRMinCalo_[MAXTRKS];
Float_t trkPtCorrCalo_[MAXTRKS];
Float_t trkPtFactCalo_[MAXTRKS];

Float_t trkRMinT_[MAXTRKS];
Float_t trkPtCorrT_[MAXTRKS];
Float_t trkPtFactT_[MAXTRKS];

Float_t trkRMinVsPF_[MAXTRKS];
Float_t trkPtCorrVsPF_[MAXTRKS];
Float_t trkPtFactVsPF_[MAXTRKS];

Float_t trkRMinVsCalo_[MAXTRKS];
Float_t trkPtCorrVsCalo_[MAXTRKS];
Float_t trkPtFactVsCalo_[MAXTRKS];

//Tracks proj. onto PF, All, Cone, and NotCone

Float_t rPFImbProjF_;
Float_t rPFImbPerpF_;
Float_t rPFImbProj0_1_;
Float_t rPFImbProj1_2_;
Float_t rPFImbProj2_4_;
Float_t rPFImbProj4_8_;
Float_t rPFImbProj8_100_;
Float_t rPFImbProjCF_;
Float_t rPFImbPerpCF_;
Float_t rPFImbProjC0_1_;
Float_t rPFImbProjC1_2_;
Float_t rPFImbProjC2_4_;
Float_t rPFImbProjC4_8_;
Float_t rPFImbProjC8_100_;
Float_t rPFImbProjNCF_;
Float_t rPFImbPerpNCF_;
Float_t rPFImbProjNC0_1_;
Float_t rPFImbProjNC1_2_;
Float_t rPFImbProjNC2_4_;
Float_t rPFImbProjNC4_8_;
Float_t rPFImbProjNC8_100_;


//Corr. Tracks proj. onto PF, All, Cone, and NotCone

Float_t rPFImbProjFCorr_;
Float_t rPFImbPerpFCorr_;
Float_t rPFImbProj0_1Corr_;
Float_t rPFImbProj1_2Corr_;
Float_t rPFImbProj2_4Corr_;
Float_t rPFImbProj4_8Corr_;
Float_t rPFImbProj8_100Corr_;
Float_t rPFImbProjCFCorr_;
Float_t rPFImbPerpCFCorr_;
Float_t rPFImbProjC0_1Corr_;
Float_t rPFImbProjC1_2Corr_;
Float_t rPFImbProjC2_4Corr_;
Float_t rPFImbProjC4_8Corr_;
Float_t rPFImbProjC8_100Corr_;
Float_t rPFImbProjNCFCorr_;
Float_t rPFImbPerpNCFCorr_;
Float_t rPFImbProjNC0_1Corr_;
Float_t rPFImbProjNC1_2Corr_;
Float_t rPFImbProjNC2_4Corr_;
Float_t rPFImbProjNC4_8Corr_;
Float_t rPFImbProjNC8_100Corr_;

//Tracks proj. onto Calo

Float_t rCaloImbProjF_;
Float_t rCaloImbPerpF_;
Float_t rCaloImbProj0_1_;
Float_t rCaloImbProj1_2_;
Float_t rCaloImbProj2_4_;
Float_t rCaloImbProj4_8_;
Float_t rCaloImbProj8_100_;
Float_t rCaloImbProjCF_;
Float_t rCaloImbPerpCF_;
Float_t rCaloImbProjC0_1_;
Float_t rCaloImbProjC1_2_;
Float_t rCaloImbProjC2_4_;
Float_t rCaloImbProjC4_8_;
Float_t rCaloImbProjC8_100_;
Float_t rCaloImbProjNCF_;
Float_t rCaloImbPerpNCF_;
Float_t rCaloImbProjNC0_1_;
Float_t rCaloImbProjNC1_2_;
Float_t rCaloImbProjNC2_4_;
Float_t rCaloImbProjNC4_8_;
Float_t rCaloImbProjNC8_100_;

//Corr. Tracks proj. onto Calo

Float_t rCaloImbProjFCorr_;
Float_t rCaloImbPerpFCorr_;
Float_t rCaloImbProj0_1Corr_;
Float_t rCaloImbProj1_2Corr_;
Float_t rCaloImbProj2_4Corr_;
Float_t rCaloImbProj4_8Corr_;
Float_t rCaloImbProj8_100Corr_;
Float_t rCaloImbProjCFCorr_;
Float_t rCaloImbPerpCFCorr_;
Float_t rCaloImbProjC0_1Corr_;
Float_t rCaloImbProjC1_2Corr_;
Float_t rCaloImbProjC2_4Corr_;
Float_t rCaloImbProjC4_8Corr_;
Float_t rCaloImbProjC8_100Corr_;
Float_t rCaloImbProjNCFCorr_;
Float_t rCaloImbPerpNCFCorr_;
Float_t rCaloImbProjNC0_1Corr_;
Float_t rCaloImbProjNC1_2Corr_;
Float_t rCaloImbProjNC2_4Corr_;
Float_t rCaloImbProjNC4_8Corr_;
Float_t rCaloImbProjNC8_100Corr_;


//Tracks proj. onto Vs PF

Float_t rVsPFImbProjF_;
Float_t rVsPFImbPerpF_;
Float_t rVsPFImbProj0_1_;
Float_t rVsPFImbProj1_2_;
Float_t rVsPFImbProj2_4_;
Float_t rVsPFImbProj4_8_;
Float_t rVsPFImbProj8_100_;
Float_t rVsPFImbProjCF_;
Float_t rVsPFImbPerpCF_;
Float_t rVsPFImbProjC0_1_;
Float_t rVsPFImbProjC1_2_;
Float_t rVsPFImbProjC2_4_;
Float_t rVsPFImbProjC4_8_;
Float_t rVsPFImbProjC8_100_;
Float_t rVsPFImbProjNCF_;
Float_t rVsPFImbPerpNCF_;
Float_t rVsPFImbProjNC0_1_;
Float_t rVsPFImbProjNC1_2_;
Float_t rVsPFImbProjNC2_4_;
Float_t rVsPFImbProjNC4_8_;
Float_t rVsPFImbProjNC8_100_;

//Corr. Tracks proj. onto Vs PF

Float_t rVsPFImbProjFCorr_;
Float_t rVsPFImbPerpFCorr_;
Float_t rVsPFImbProj0_1Corr_;
Float_t rVsPFImbProj1_2Corr_;
Float_t rVsPFImbProj2_4Corr_;
Float_t rVsPFImbProj4_8Corr_;
Float_t rVsPFImbProj8_100Corr_;
Float_t rVsPFImbProjCFCorr_;
Float_t rVsPFImbPerpCFCorr_;
Float_t rVsPFImbProjC0_1Corr_;
Float_t rVsPFImbProjC1_2Corr_;
Float_t rVsPFImbProjC2_4Corr_;
Float_t rVsPFImbProjC4_8Corr_;
Float_t rVsPFImbProjC8_100Corr_;
Float_t rVsPFImbProjNCFCorr_;
Float_t rVsPFImbPerpNCFCorr_;
Float_t rVsPFImbProjNC0_1Corr_;
Float_t rVsPFImbProjNC1_2Corr_;
Float_t rVsPFImbProjNC2_4Corr_;
Float_t rVsPFImbProjNC4_8Corr_;
Float_t rVsPFImbProjNC8_100Corr_;

//Tracks proj. onto Vs Calo

Float_t rVsCaloImbProjF_;
Float_t rVsCaloImbPerpF_;
Float_t rVsCaloImbProj0_1_;
Float_t rVsCaloImbProj1_2_;
Float_t rVsCaloImbProj2_4_;
Float_t rVsCaloImbProj4_8_;
Float_t rVsCaloImbProj8_100_;
Float_t rVsCaloImbProjCF_;
Float_t rVsCaloImbPerpCF_;
Float_t rVsCaloImbProjC0_1_;
Float_t rVsCaloImbProjC1_2_;
Float_t rVsCaloImbProjC2_4_;
Float_t rVsCaloImbProjC4_8_;
Float_t rVsCaloImbProjC8_100_;
Float_t rVsCaloImbProjNCF_;
Float_t rVsCaloImbPerpNCF_;
Float_t rVsCaloImbProjNC0_1_;
Float_t rVsCaloImbProjNC1_2_;
Float_t rVsCaloImbProjNC2_4_;
Float_t rVsCaloImbProjNC4_8_;
Float_t rVsCaloImbProjNC8_100_;

//Corr. Tracks proj. onto Vs Calo

Float_t rVsCaloImbProjFCorr_;
Float_t rVsCaloImbPerpFCorr_;
Float_t rVsCaloImbProj0_1Corr_;
Float_t rVsCaloImbProj1_2Corr_;
Float_t rVsCaloImbProj2_4Corr_;
Float_t rVsCaloImbProj4_8Corr_;
Float_t rVsCaloImbProj8_100Corr_;
Float_t rVsCaloImbProjCFCorr_;
Float_t rVsCaloImbPerpCFCorr_;
Float_t rVsCaloImbProjC0_1Corr_;
Float_t rVsCaloImbProjC1_2Corr_;
Float_t rVsCaloImbProjC2_4Corr_;
Float_t rVsCaloImbProjC4_8Corr_;
Float_t rVsCaloImbProjC8_100Corr_;
Float_t rVsCaloImbProjNCFCorr_;
Float_t rVsCaloImbPerpNCFCorr_;
Float_t rVsCaloImbProjNC0_1Corr_;
Float_t rVsCaloImbProjNC1_2Corr_;
Float_t rVsCaloImbProjNC2_4Corr_;
Float_t rVsCaloImbProjNC4_8Corr_;
Float_t rVsCaloImbProjNC8_100Corr_;

//Tracks proj. onto Truth

Float_t rTImbPerpF_;
Float_t rTImbProjF_;
Float_t rTImbProj0_1_;
Float_t rTImbProj1_2_;
Float_t rTImbProj2_4_;
Float_t rTImbProj4_8_;
Float_t rTImbProj8_100_;
Float_t rTImbProjCF_;
Float_t rTImbPerpCF_;
Float_t rTImbProjC0_1_;
Float_t rTImbProjC1_2_;
Float_t rTImbProjC2_4_;
Float_t rTImbProjC4_8_;
Float_t rTImbProjC8_100_;
Float_t rTImbProjNCF_;
Float_t rTImbPerpNCF_;
Float_t rTImbProjNC0_1_;
Float_t rTImbProjNC1_2_;
Float_t rTImbProjNC2_4_;
Float_t rTImbProjNC4_8_;
Float_t rTImbProjNC8_100_;

//Corr. Tracks proj. onto Truth

Float_t rTImbPerpFCorr_;
Float_t rTImbProjFCorr_;
Float_t rTImbProj0_1Corr_;
Float_t rTImbProj1_2Corr_;
Float_t rTImbProj2_4Corr_;
Float_t rTImbProj4_8Corr_;
Float_t rTImbProj8_100Corr_;
Float_t rTImbProjCFCorr_;
Float_t rTImbPerpCFCorr_;
Float_t rTImbProjC0_1Corr_;
Float_t rTImbProjC1_2Corr_;
Float_t rTImbProjC2_4Corr_;
Float_t rTImbProjC4_8Corr_;
Float_t rTImbProjC8_100Corr_;
Float_t rTImbProjNCFCorr_;
Float_t rTImbPerpNCFCorr_;
Float_t rTImbProjNC0_1Corr_;
Float_t rTImbProjNC1_2Corr_;
Float_t rTImbProjNC2_4Corr_;
Float_t rTImbProjNC4_8Corr_;
Float_t rTImbProjNC8_100Corr_;



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

//Set Bool

Bool_t truthSet_;
Bool_t recoPFSet_;
Bool_t recoVsPFSet_;
Bool_t recoCaloSet_;
Bool_t recoVsCaloSet_;

//Truth Jt

Float_t TLeadJtPt_;
Float_t TLeadJtPhi_;
Float_t TLeadJtEta_;
Float_t TSubLeadJtPt_;
Float_t TSubLeadJtPhi_;
Float_t TSubLeadJtEta_;
Float_t TJtAsymm_;

//PF Jt

Float_t PFLeadJtPt_;
Float_t PFLeadJtPhi_;
Float_t PFLeadJtEta_;
Float_t PFSubLeadJtPt_;
Float_t PFSubLeadJtPhi_;
Float_t PFSubLeadJtEta_;
Float_t PFJtAsymm_;

//Calo Jt

Float_t CaloLeadJtPt_;
Float_t CaloLeadJtPhi_;
Float_t CaloLeadJtEta_;
Float_t CaloSubLeadJtPt_;
Float_t CaloSubLeadJtPhi_;
Float_t CaloSubLeadJtEta_;
Float_t CaloJtAsymm_;

//Vs PF Jt

Float_t VsPFLeadJtPt_;
Float_t VsPFLeadJtPhi_;
Float_t VsPFLeadJtEta_;
Float_t VsPFSubLeadJtPt_;
Float_t VsPFSubLeadJtPhi_;
Float_t VsPFSubLeadJtEta_;
Float_t VsPFJtAsymm_;

//Vs Calo Jt

Float_t VsCaloLeadJtPt_;
Float_t VsCaloLeadJtPhi_;
Float_t VsCaloLeadJtEta_;
Float_t VsCaloSubLeadJtPt_;
Float_t VsCaloSubLeadJtPhi_;
Float_t VsCaloSubLeadJtEta_;
Float_t VsCaloJtAsymm_;

//Gen Tree Variables

const int MAXGEN = 50000; //From SetupGenParticleTree.h

Int_t nGen_;
Float_t genPt_[MAXGEN];
Float_t genPtPF_[MAXGEN];
Float_t genPtCalo_[MAXGEN];
Float_t genPtT_[MAXGEN];
Float_t genPtVsPF_[MAXGEN];
Float_t genPtVsCalo_[MAXGEN];
Float_t genPhi_[MAXGEN];
Float_t genEta_[MAXGEN];
Float_t genLeadDelPhi_[MAXGEN];

//Gen. proj. onto Truth

Float_t gTImbProjF_;
Float_t gTImbPerpF_;

Float_t gTImbProj0_1_;
Float_t gTImbProj1_2_;
Float_t gTImbProj2_4_;
Float_t gTImbProj4_8_;
Float_t gTImbProj8_100_;

//Gen. proj. onto PF

Float_t gPFImbProjF_;
Float_t gPFImbPerpF_;

Float_t gPFImbProj0_1_;
Float_t gPFImbProj1_2_;
Float_t gPFImbProj2_4_;
Float_t gPFImbProj4_8_;
Float_t gPFImbProj8_100_;

//Gen. proj. onto Calo

Float_t gCaloImbProjF_;
Float_t gCaloImbPerpF_;

Float_t gCaloImbProj0_1_;
Float_t gCaloImbProj1_2_;
Float_t gCaloImbProj2_4_;
Float_t gCaloImbProj4_8_;
Float_t gCaloImbProj8_100_;


//Gen. proj. onto Vs PF

Float_t gVsPFImbProjF_;
Float_t gVsPFImbPerpF_;

Float_t gVsPFImbProj0_1_;
Float_t gVsPFImbProj1_2_;
Float_t gVsPFImbProj2_4_;
Float_t gVsPFImbProj4_8_;
Float_t gVsPFImbProj8_100_;

//Gen. proj. onto Vs Calo

Float_t gVsCaloImbProjF_;
Float_t gVsCaloImbPerpF_;

Float_t gVsCaloImbProj0_1_;
Float_t gVsCaloImbProj1_2_;
Float_t gVsCaloImbProj2_4_;
Float_t gVsCaloImbProj4_8_;
Float_t gVsCaloImbProj8_100_;


void SetBranches(bool montecarlo)
{
  //Track Tree Branches

  std::cout << "Branches Set" << std::endl;
  
  trackTree_p->Branch("nTrk", &nTrk_, "nTrk/I");
  trackTree_p->Branch("trkPt", &trkPt_, "trkPt[nTrk]/F");
  trackTree_p->Branch("trkPtPF", &trkPtPF_, "trkPtPF[nTrk]/F");
  trackTree_p->Branch("trkPtCalo", &trkPtCalo_, "trkPtCalo[nTrk]/F");
  trackTree_p->Branch("trkPtVsPF", &trkPtVsPF_, "trkPtVsPF[nTrk]/F");
  trackTree_p->Branch("trkPtVsCalo", &trkPtVsCalo_, "trkPtVsCalo[nTrk]/F");

  if(montecarlo)
    trackTree_p->Branch("trkPtT", &trkPtT_, "trkPtT[nTrk]/F");

  trackTree_p->Branch("trkPhi", &trkPhi_, "trkPhi[nTrk]/F");
  trackTree_p->Branch("trkEta", &trkEta_, "trkEta[nTrk]/F");
  trackTree_p->Branch("trkPFLeadDelPhi", &trkPFLeadDelPhi_, "trkPFLeadDelPhi[nTrk]/F");
  trackTree_p->Branch("trkCaloLeadDelPhi", &trkCaloLeadDelPhi_, "trkCaloLeadDelPhi[nTrk]/F");

  trackTree_p->Branch("trkRLeadPF", &trkRLeadPF_, "trkRLeadPF[nTrk]/F");
  trackTree_p->Branch("trkRMinPF", &trkRMinPF_, "trkRMinPF[nTrk]/F");
  trackTree_p->Branch("trkPtCorrPF", &trkPtCorrPF_, "trkPtCorrPF[nTrk]/F");
  trackTree_p->Branch("trkPtFactPF", &trkPtFactPF_, "trkPtFactPF[nTrk]/F");
  trackTree_p->Branch("trkRMinCalo", &trkRMinCalo_, "trkRMinCalo[nTrk]/F");
  trackTree_p->Branch("trkPtCorrCalo", &trkPtCorrCalo_, "trkPtCorrCalo[nTrk]/F");
  trackTree_p->Branch("trkPtFactCalo", &trkPtFactCalo_, "trkPtFactCalo[nTrk]/F");

  if(montecarlo){
    trackTree_p->Branch("trkRMinT", &trkRMinT_, "trkRMinT[nTrk]/F");
    trackTree_p->Branch("trkPtCorrT", &trkPtCorrT_, "trkPtCorrT[nTrk]/F");
    trackTree_p->Branch("trkPtFactT", &trkPtFactT_, "trkPtFactT[nTrk]/F");
  }

  trackTree_p->Branch("trkRMinVsPF", &trkRMinVsPF_, "trkRMinVsPF[nTrk]/F");
  trackTree_p->Branch("trkPtCorrVsPF", &trkPtCorrVsPF_, "trkPtCorrVsPF[nTrk]/F");
  trackTree_p->Branch("trkPtFactVsPF", &trkPtFactVsPF_, "trkPtFactVsPF[nTrk]/F");
  trackTree_p->Branch("trkRMinVsCalo", &trkRMinVsCalo_, "trkRMinVsCalo[nTrk]/F");
  trackTree_p->Branch("trkPtCorrVsCalo", &trkPtCorrVsCalo_, "trkPtCorrVsCalo[nTrk]/F");
  trackTree_p->Branch("trkPtFactVsCalo", &trkPtFactVsCalo_, "trkPtFactVsCalo[nTrk]/F");  

  //Tracks proj. onto PF, All, Cone, and NotCone

  trackTree_p->Branch("rPFImbProjF", &rPFImbProjF_, "rPFImbProjF/F");
  trackTree_p->Branch("rPFImbPerpF", &rPFImbPerpF_, "rPFImbPerpF/F");
  trackTree_p->Branch("rPFImbProj0_1", &rPFImbProj0_1_, "rPFImbProj0_1/F");
  trackTree_p->Branch("rPFImbProj1_2", &rPFImbProj1_2_, "rPFImbProj1_2/F");
  trackTree_p->Branch("rPFImbProj2_4", &rPFImbProj2_4_, "rPFImbProj2_4/F");
  trackTree_p->Branch("rPFImbProj4_8", &rPFImbProj4_8_, "rPFImbProj4_8/F");
  trackTree_p->Branch("rPFImbProj8_100", &rPFImbProj8_100_, "rPFImbProj8_100/F");

  trackTree_p->Branch("rPFImbProjCF", &rPFImbProjCF_, "rPFImbProjCF/F");
  trackTree_p->Branch("rPFImbPerpCF", &rPFImbPerpCF_, "rPFImbPerpCF/F");
  trackTree_p->Branch("rPFImbProjC0_1", &rPFImbProjC0_1_, "rPFImbProjC0_1/F");
  trackTree_p->Branch("rPFImbProjC1_2", &rPFImbProjC1_2_, "rPFImbProjC1_2/F");
  trackTree_p->Branch("rPFImbProjC2_4", &rPFImbProjC2_4_, "rPFImbProjC2_4/F");
  trackTree_p->Branch("rPFImbProjC4_8", &rPFImbProjC4_8_, "rPFImbProjC4_8/F");
  trackTree_p->Branch("rPFImbProjC8_100", &rPFImbProjC8_100_, "rPFImbProjC8_100/F");
  trackTree_p->Branch("rPFImbProjNCF", &rPFImbProjNCF_, "rPFImbProjNCF/F");
  trackTree_p->Branch("rPFImbPerpNCF", &rPFImbPerpNCF_, "rPFImbPerpNCF/F");
  trackTree_p->Branch("rPFImbProjNC0_1", &rPFImbProjNC0_1_, "rPFImbProjNC0_1/F");
  trackTree_p->Branch("rPFImbProjNC1_2", &rPFImbProjNC1_2_, "rPFImbProjNC1_2/F");
  trackTree_p->Branch("rPFImbProjNC2_4", &rPFImbProjNC2_4_, "rPFImbProjNC2_4/F");
  trackTree_p->Branch("rPFImbProjNC4_8", &rPFImbProjNC4_8_, "rPFImbProjNC4_8/F");
  trackTree_p->Branch("rPFImbProjNC8_100", &rPFImbProjNC8_100_, "rPFImbProjNC8_100/F");


  //Corr. Tracks proj. onto PF, All, Cone, and NotCone

  trackTree_p->Branch("rPFImbProjFCorr", &rPFImbProjFCorr_, "rPFImbProjFCorr/F");
  trackTree_p->Branch("rPFImbPerpFCorr", &rPFImbPerpFCorr_, "rPFImbPerpFCorr/F");
  trackTree_p->Branch("rPFImbProj0_1Corr", &rPFImbProj0_1Corr_, "rPFImbProj0_1Corr/F");
  trackTree_p->Branch("rPFImbProj1_2Corr", &rPFImbProj1_2Corr_, "rPFImbProj1_2Corr/F");
  trackTree_p->Branch("rPFImbProj2_4Corr", &rPFImbProj2_4Corr_, "rPFImbProj2_4Corr/F");
  trackTree_p->Branch("rPFImbProj4_8Corr", &rPFImbProj4_8Corr_, "rPFImbProj4_8Corr/F");
  trackTree_p->Branch("rPFImbProj8_100Corr", &rPFImbProj8_100Corr_, "rPFImbProj8_100Corr/F");
  trackTree_p->Branch("rPFImbProjCFCorr", &rPFImbProjCFCorr_, "rPFImbProjCFCorr/F");
  trackTree_p->Branch("rPFImbPerpCFCorr", &rPFImbPerpCFCorr_, "rPFImbPerpCFCorr/F");
  trackTree_p->Branch("rPFImbProjC0_1Corr", &rPFImbProjC0_1Corr_, "rPFImbProjC0_1Corr/F");
  trackTree_p->Branch("rPFImbProjC1_2Corr", &rPFImbProjC1_2Corr_, "rPFImbProjC1_2Corr/F");
  trackTree_p->Branch("rPFImbProjC2_4Corr", &rPFImbProjC2_4Corr_, "rPFImbProjC2_4Corr/F");
  trackTree_p->Branch("rPFImbProjC4_8Corr", &rPFImbProjC4_8Corr_, "rPFImbProjC4_8Corr/F");
  trackTree_p->Branch("rPFImbProjC8_100Corr", &rPFImbProjC8_100Corr_, "rPFImbProjC8_100Corr/F");
  trackTree_p->Branch("rPFImbProjNCFCorr", &rPFImbProjNCFCorr_, "rPFImbProjNCFCorr/F");
  trackTree_p->Branch("rPFImbPerpNCFCorr", &rPFImbPerpNCFCorr_, "rPFImbPerpNCFCorr/F");
  trackTree_p->Branch("rPFImbProjNC0_1Corr", &rPFImbProjNC0_1Corr_, "rPFImbProjNC0_1Corr/F");
  trackTree_p->Branch("rPFImbProjNC1_2Corr", &rPFImbProjNC1_2Corr_, "rPFImbProjNC1_2Corr/F");
  trackTree_p->Branch("rPFImbProjNC2_4Corr", &rPFImbProjNC2_4Corr_, "rPFImbProjNC2_4Corr/F");
  trackTree_p->Branch("rPFImbProjNC4_8Corr", &rPFImbProjNC4_8Corr_, "rPFImbProjNC4_8Corr/F");
  trackTree_p->Branch("rPFImbProjNC8_100Corr", &rPFImbProjNC8_100Corr_, "rPFImbProjNC8_100Corr/F");

  //Tracks proj. onto Calo

  trackTree_p->Branch("rCaloImbProjF", &rCaloImbProjF_, "rCaloImbProjF/F");
  trackTree_p->Branch("rCaloImbPerpF", &rCaloImbPerpF_, "rCaloImbPerpF/F");

  trackTree_p->Branch("rCaloImbProj0_1", &rCaloImbProj0_1_, "rCaloImbProj0_1/F");
  trackTree_p->Branch("rCaloImbProj1_2", &rCaloImbProj1_2_, "rCaloImbProj1_2/F");
  trackTree_p->Branch("rCaloImbProj2_4", &rCaloImbProj2_4_, "rCaloImbProj2_4/F");
  trackTree_p->Branch("rCaloImbProj4_8", &rCaloImbProj4_8_, "rCaloImbProj4_8/F");
  trackTree_p->Branch("rCaloImbProj8_100", &rCaloImbProj8_100_, "rCaloImbProj8_100/F");

  //Corr. Tracks proj. onto Calo

  trackTree_p->Branch("rCaloImbProjFCorr", &rCaloImbProjFCorr_, "rCaloImbProjFCorr/F");
  trackTree_p->Branch("rCaloImbPerpFCorr", &rCaloImbPerpFCorr_, "rCaloImbPerpFCorr/F");

  trackTree_p->Branch("rCaloImbProj0_1Corr", &rCaloImbProj0_1Corr_, "rCaloImbProj0_1Corr/F");
  trackTree_p->Branch("rCaloImbProj1_2Corr", &rCaloImbProj1_2Corr_, "rCaloImbProj1_2Corr/F");
  trackTree_p->Branch("rCaloImbProj2_4Corr", &rCaloImbProj2_4Corr_, "rCaloImbProj2_4Corr/F");
  trackTree_p->Branch("rCaloImbProj4_8Corr", &rCaloImbProj4_8Corr_, "rCaloImbProj4_8Corr/F");
  trackTree_p->Branch("rCaloImbProj8_100Corr", &rCaloImbProj8_100Corr_, "rCaloImbProj8_100Corr/F");

  //Tracks proj. onto Vs PF

  trackTree_p->Branch("rVsPFImbProjF", &rVsPFImbProjF_, "rVsPFImbProjF/F");
  trackTree_p->Branch("rVsPFImbPerpF", &rVsPFImbPerpF_, "rVsPFImbPerpF/F");

  trackTree_p->Branch("rVsPFImbProj0_1", &rVsPFImbProj0_1_, "rVsPFImbProj0_1/F");
  trackTree_p->Branch("rVsPFImbProj1_2", &rVsPFImbProj1_2_, "rVsPFImbProj1_2/F");
  trackTree_p->Branch("rVsPFImbProj2_4", &rVsPFImbProj2_4_, "rVsPFImbProj2_4/F");
  trackTree_p->Branch("rVsPFImbProj4_8", &rVsPFImbProj4_8_, "rVsPFImbProj4_8/F");
  trackTree_p->Branch("rVsPFImbProj8_100", &rVsPFImbProj8_100_, "rVsPFImbProj8_100/F");

  //Corr. Tracks proj. onto Vs PF

  trackTree_p->Branch("rVsPFImbProjFCorr", &rVsPFImbProjFCorr_, "rVsPFImbProjFCorr/F");
  trackTree_p->Branch("rVsPFImbPerpFCorr", &rVsPFImbPerpFCorr_, "rVsPFImbPerpFCorr/F");

  trackTree_p->Branch("rVsPFImbProj0_1Corr", &rVsPFImbProj0_1Corr_, "rVsPFImbProj0_1Corr/F");
  trackTree_p->Branch("rVsPFImbProj1_2Corr", &rVsPFImbProj1_2Corr_, "rVsPFImbProj1_2Corr/F");
  trackTree_p->Branch("rVsPFImbProj2_4Corr", &rVsPFImbProj2_4Corr_, "rVsPFImbProj2_4Corr/F");
  trackTree_p->Branch("rVsPFImbProj4_8Corr", &rVsPFImbProj4_8Corr_, "rVsPFImbProj4_8Corr/F");
  trackTree_p->Branch("rVsPFImbProj8_100Corr", &rVsPFImbProj8_100Corr_, "rVsPFImbProj8_100Corr/F");

  //Tracks proj. onto Vs Calo

  trackTree_p->Branch("rVsCaloImbProjF", &rVsCaloImbProjF_, "rVsCaloImbProjF/F");
  trackTree_p->Branch("rVsCaloImbPerpF", &rVsCaloImbPerpF_, "rVsCaloImbPerpF/F");

  trackTree_p->Branch("rVsCaloImbProj0_1", &rVsCaloImbProj0_1_, "rVsCaloImbProj0_1/F");
  trackTree_p->Branch("rVsCaloImbProj1_2", &rVsCaloImbProj1_2_, "rVsCaloImbProj1_2/F");
  trackTree_p->Branch("rVsCaloImbProj2_4", &rVsCaloImbProj2_4_, "rVsCaloImbProj2_4/F");
  trackTree_p->Branch("rVsCaloImbProj4_8", &rVsCaloImbProj4_8_, "rVsCaloImbProj4_8/F");
  trackTree_p->Branch("rVsCaloImbProj8_100", &rVsCaloImbProj8_100_, "rVsCaloImbProj8_100/F");

  //Corr. Tracks proj. onto Vs Calo

  trackTree_p->Branch("rVsCaloImbProjFCorr", &rVsCaloImbProjFCorr_, "rVsCaloImbProjFCorr/F");
  trackTree_p->Branch("rVsCaloImbPerpFCorr", &rVsCaloImbPerpFCorr_, "rVsCaloImbPerpFCorr/F");

  trackTree_p->Branch("rVsCaloImbProj0_1Corr", &rVsCaloImbProj0_1Corr_, "rVsCaloImbProj0_1Corr/F");
  trackTree_p->Branch("rVsCaloImbProj1_2Corr", &rVsCaloImbProj1_2Corr_, "rVsCaloImbProj1_2Corr/F");
  trackTree_p->Branch("rVsCaloImbProj2_4Corr", &rVsCaloImbProj2_4Corr_, "rVsCaloImbProj2_4Corr/F");
  trackTree_p->Branch("rVsCaloImbProj4_8Corr", &rVsCaloImbProj4_8Corr_, "rVsCaloImbProj4_8Corr/F");
  trackTree_p->Branch("rVsCaloImbProj8_100Corr", &rVsCaloImbProj8_100Corr_, "rVsCaloImbProj8_100Corr/F");

  if(montecarlo){
    //Track tree branches iff truth avail.

    //Tracks proj. onto Truth

    trackTree_p->Branch("rTImbProjF", &rTImbProjF_, "rTImbProjF/F");
    trackTree_p->Branch("rTImbPerpF", &rTImbPerpF_, "rTImbPerpF/F");

    trackTree_p->Branch("rTImbProj0_1", &rTImbProj0_1_, "rTImbProj0_1/F");
    trackTree_p->Branch("rTImbProj1_2", &rTImbProj1_2_, "rTImbProj1_2/F");
    trackTree_p->Branch("rTImbProj2_4", &rTImbProj2_4_, "rTImbProj2_4/F");
    trackTree_p->Branch("rTImbProj4_8", &rTImbProj4_8_, "rTImbProj4_8/F");
    trackTree_p->Branch("rTImbProj8_100", &rTImbProj8_100_, "rTImbProj8_100/F");

    //Corr. Tracks proj. onto Truth

    trackTree_p->Branch("rTImbProjFCorr", &rTImbProjFCorr_, "rTImbProjFCorr/F");
    trackTree_p->Branch("rTImbPerpFCorr", &rTImbPerpFCorr_, "rTImbPerpFCorr/F");

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
  jetTree_p->Branch("recoVsPFSet", &recoVsPFSet_, "recoVsPFSet/O");
  jetTree_p->Branch("recoCaloSet", &recoCaloSet_, "recoCaloSet/O");
  jetTree_p->Branch("recoVsCaloSet", &recoVsCaloSet_, "recoVsCaloSet/O");

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

  jetTree_p->Branch("VsPFLeadJtPt", &VsPFLeadJtPt_, "VsPFLeadJtPt/F");
  jetTree_p->Branch("VsPFLeadJtPhi", &VsPFLeadJtPhi_, "VsPFLeadJtPhi/F");
  jetTree_p->Branch("VsPFLeadJtEta", &VsPFLeadJtEta_, "VsPFLeadJtEta/F");
  jetTree_p->Branch("VsPFSubLeadJtPt", &VsPFSubLeadJtPt_, "VsPFSubLeadJtPt/F");
  jetTree_p->Branch("VsPFSubLeadJtPhi", &VsPFSubLeadJtPhi_, "VsPFSubLeadJtPhi/F");
  jetTree_p->Branch("VsPFSubLeadJtEta", &VsPFSubLeadJtEta_, "VsPFSubLeadJtEta/F");
  jetTree_p->Branch("VsPFJtAsymm", &VsPFJtAsymm_, "VsPFJtAsymm/F");
    
  jetTree_p->Branch("VsCaloLeadJtPt", &VsCaloLeadJtPt_, "VsCaloLeadJtPt/F");
  jetTree_p->Branch("VsCaloLeadJtPhi", &VsCaloLeadJtPhi_, "VsCaloLeadJtPhi/F");
  jetTree_p->Branch("VsCaloLeadJtEta", &VsCaloLeadJtEta_, "VsCaloLeadJtEta/F");
  jetTree_p->Branch("VsCaloSubLeadJtPt", &VsCaloSubLeadJtPt_, "VsCaloSubLeadJtPt/F");
  jetTree_p->Branch("VsCaloSubLeadJtPhi", &VsCaloSubLeadJtPhi_, "VsCaloSubLeadJtPhi/F");
  jetTree_p->Branch("VsCaloSubLeadJtEta", &VsCaloSubLeadJtEta_, "VsCaloSubLeadJtEta/F");
  jetTree_p->Branch("VsCaloJtAsymm", &VsCaloJtAsymm_, "VsCaloJtAsymm/F");

  if(montecarlo){
    //Gen Tree Branches

    genTree_p->Branch("nGen", &nGen_, "nGen/I");
    genTree_p->Branch("genPt", &genPt_, "genPt[nGen]/F");
    genTree_p->Branch("genPtPF", &genPtPF_, "genPtPF[nGen]/F");
    genTree_p->Branch("genPtCalo", &genPtCalo_, "genPtCalo[nGen]/F");
    genTree_p->Branch("genPtT", &genPtT_, "genPtT[nGen]/F");
    genTree_p->Branch("genPtVsPF", &genPtVsPF_, "genPtVsPF[nGen]/F");
    genTree_p->Branch("genPtVsCalo", &genPtVsCalo_, "genPtVsCalo[nGen]/F");
    genTree_p->Branch("genPhi", &genPhi_, "genPhi[nGen]/F");
    genTree_p->Branch("genEta", &genEta_, "genEta[nGen]/F");
    genTree_p->Branch("genLeadDelPhi", &genLeadDelPhi_, "genLeadDelPhi[nGen]/F");

    //Gen. proj. onto Truth

    genTree_p->Branch("gTImbProjF", &gTImbProjF_, "gTImbProjF/F");
    genTree_p->Branch("gTImbPerpF", &gTImbPerpF_, "gTImbPerpF/F");

    genTree_p->Branch("gTImbProj0_1", &gTImbProj0_1_, "gTImbProj0_1/F");
    genTree_p->Branch("gTImbProj1_2", &gTImbProj1_2_, "gTImbProj1_2/F");
    genTree_p->Branch("gTImbProj2_4", &gTImbProj2_4_, "gTImbProj2_4/F");
    genTree_p->Branch("gTImbProj4_8", &gTImbProj4_8_, "gTImbProj4_8/F");
    genTree_p->Branch("gTImbProj8_100", &gTImbProj8_100_, "gTImbProj8_100/F");

    //Gen. proj. onto PF

    genTree_p->Branch("gPFImbProjF", &gPFImbProjF_, "gPFImbProjF/F");
    genTree_p->Branch("gPFImbPerpF", &gPFImbPerpF_, "gPFImbPerpF/F");

    genTree_p->Branch("gPFImbProj0_1", &gPFImbProj0_1_, "gPFImbProj0_1/F");
    genTree_p->Branch("gPFImbProj1_2", &gPFImbProj1_2_, "gPFImbProj1_2/F");
    genTree_p->Branch("gPFImbProj2_4", &gPFImbProj2_4_, "gPFImbProj2_4/F");
    genTree_p->Branch("gPFImbProj4_8", &gPFImbProj4_8_, "gPFImbProj4_8/F");
    genTree_p->Branch("gPFImbProj8_100", &gPFImbProj8_100_, "gPFImbProj8_100/F");

    //Gen. proj. onto Calo

    genTree_p->Branch("gCaloImbProjF", &gCaloImbProjF_, "gCaloImbProjF/F");
    genTree_p->Branch("gCaloImbPerpF", &gCaloImbPerpF_, "gCaloImbPerpF/F");

    genTree_p->Branch("gCaloImbProj0_1", &gCaloImbProj0_1_, "gCaloImbProj0_1/F");
    genTree_p->Branch("gCaloImbProj1_2", &gCaloImbProj1_2_, "gCaloImbProj1_2/F");
    genTree_p->Branch("gCaloImbProj2_4", &gCaloImbProj2_4_, "gCaloImbProj2_4/F");
    genTree_p->Branch("gCaloImbProj4_8", &gCaloImbProj4_8_, "gCaloImbProj4_8/F");
    genTree_p->Branch("gCaloImbProj8_100", &gCaloImbProj8_100_, "gCaloImbProj8_100/F");

    //Gen. proj. onto Vs PF

    genTree_p->Branch("gVsPFImbProjF", &gVsPFImbProjF_, "gVsPFImbProjF/F");
    genTree_p->Branch("gVsPFImbPerpF", &gVsPFImbPerpF_, "gVsPFImbPerpF/F");

    genTree_p->Branch("gVsPFImbProj0_1", &gVsPFImbProj0_1_, "gVsPFImbProj0_1/F");
    genTree_p->Branch("gVsPFImbProj1_2", &gVsPFImbProj1_2_, "gVsPFImbProj1_2/F");
    genTree_p->Branch("gVsPFImbProj2_4", &gVsPFImbProj2_4_, "gVsPFImbProj2_4/F");
    genTree_p->Branch("gVsPFImbProj4_8", &gVsPFImbProj4_8_, "gVsPFImbProj4_8/F");
    genTree_p->Branch("gVsPFImbProj8_100", &gVsPFImbProj8_100_, "gVsPFImbProj8_100/F");

    //Gen. proj. onto Vs Calo

    genTree_p->Branch("gVsCaloImbProjF", &gVsCaloImbProjF_, "gVsCaloImbProjF/F");
    genTree_p->Branch("gVsCaloImbPerpF", &gVsCaloImbPerpF_, "gVsCaloImbPerpF/F");

    genTree_p->Branch("gVsCaloImbProj0_1", &gVsCaloImbProj0_1_, "gVsCaloImbProj0_1/F");
    genTree_p->Branch("gVsCaloImbProj1_2", &gVsCaloImbProj1_2_, "gVsCaloImbProj1_2/F");
    genTree_p->Branch("gVsCaloImbProj2_4", &gVsCaloImbProj2_4_, "gVsCaloImbProj2_4/F");
    genTree_p->Branch("gVsCaloImbProj4_8", &gVsCaloImbProj4_8_, "gVsCaloImbProj4_8/F");
    genTree_p->Branch("gVsCaloImbProj8_100", &gVsCaloImbProj8_100_, "gVsCaloImbProj8_100/F");
  }
}


void GetBranches(bool montecarlo)
{
  //Track Tree Branches

  trackTree_p->SetBranchAddress("nTrk", &nTrk_);
  trackTree_p->SetBranchAddress("trkPt", &trkPt_);
  trackTree_p->SetBranchAddress("trkPtPF", &trkPtPF_);
  trackTree_p->SetBranchAddress("trkPtCalo", &trkPtCalo_);
  trackTree_p->SetBranchAddress("trkPtVsPF", &trkPtVsPF_);
  trackTree_p->SetBranchAddress("trkPtVsCalo", &trkPtVsCalo_);

  if(montecarlo)
    trackTree_p->SetBranchAddress("trkPtT", &trkPtT_);

  trackTree_p->SetBranchAddress("trkPhi", &trkPhi_);
  trackTree_p->SetBranchAddress("trkEta", &trkEta_);
  trackTree_p->SetBranchAddress("trkPFLeadDelPhi", &trkPFLeadDelPhi_);
  trackTree_p->SetBranchAddress("trkCaloLeadDelPhi", &trkCaloLeadDelPhi_);

  trackTree_p->SetBranchAddress("trkRLeadPF", &trkRLeadPF_);
  trackTree_p->SetBranchAddress("trkRMinPF", &trkRMinPF_);
  trackTree_p->SetBranchAddress("trkPtCorrPF", &trkPtCorrPF_);
  trackTree_p->SetBranchAddress("trkPtFactPF", &trkPtFactPF_);
  trackTree_p->SetBranchAddress("trkRMinCalo", &trkRMinCalo_);
  trackTree_p->SetBranchAddress("trkPtCorrCalo", &trkPtCorrCalo_);
  trackTree_p->SetBranchAddress("trkPtFactCalo", &trkPtFactCalo_);

  if(montecarlo){
  trackTree_p->SetBranchAddress("trkRMinT", &trkRMinT_);
  trackTree_p->SetBranchAddress("trkPtCorrT", &trkPtCorrT_);
  trackTree_p->SetBranchAddress("trkPtFactT", &trkPtFactT_);
  }

  trackTree_p->SetBranchAddress("trkRMinVsPF", &trkRMinVsPF_);
  trackTree_p->SetBranchAddress("trkPtCorrVsPF", &trkPtCorrVsPF_);
  trackTree_p->SetBranchAddress("trkPtFactVsPF", &trkPtFactVsPF_);
  trackTree_p->SetBranchAddress("trkRMinVsCalo", &trkRMinVsCalo_);
  trackTree_p->SetBranchAddress("trkPtCorrVsCalo", &trkPtCorrVsCalo_);
  trackTree_p->SetBranchAddress("trkPtFactVsCalo", &trkPtFactVsCalo_);

  //Tracks proj. onto PF, All, Cone, and NotCone

  trackTree_p->SetBranchAddress("rPFImbProjF", &rPFImbProjF_);
  trackTree_p->SetBranchAddress("rPFImbPerpF", &rPFImbPerpF_);
  trackTree_p->SetBranchAddress("rPFImbProj0_1", &rPFImbProj0_1_);
  trackTree_p->SetBranchAddress("rPFImbProj1_2", &rPFImbProj1_2_);
  trackTree_p->SetBranchAddress("rPFImbProj2_4", &rPFImbProj2_4_);
  trackTree_p->SetBranchAddress("rPFImbProj4_8", &rPFImbProj4_8_);
  trackTree_p->SetBranchAddress("rPFImbProj8_100", &rPFImbProj8_100_);
  trackTree_p->SetBranchAddress("rPFImbProjCF", &rPFImbProjCF_);
  trackTree_p->SetBranchAddress("rPFImbPerpCF", &rPFImbPerpCF_);
  trackTree_p->SetBranchAddress("rPFImbProjC0_1", &rPFImbProjC0_1_);
  trackTree_p->SetBranchAddress("rPFImbProjC1_2", &rPFImbProjC1_2_);
  trackTree_p->SetBranchAddress("rPFImbProjC2_4", &rPFImbProjC2_4_);
  trackTree_p->SetBranchAddress("rPFImbProjC4_8", &rPFImbProjC4_8_);
  trackTree_p->SetBranchAddress("rPFImbProjC8_100", &rPFImbProjC8_100_);
  trackTree_p->SetBranchAddress("rPFImbProjNCF", &rPFImbProjNCF_);
  trackTree_p->SetBranchAddress("rPFImbPerpNCF", &rPFImbPerpNCF_);
  trackTree_p->SetBranchAddress("rPFImbProjNC0_1", &rPFImbProjNC0_1_);
  trackTree_p->SetBranchAddress("rPFImbProjNC1_2", &rPFImbProjNC1_2_);
  trackTree_p->SetBranchAddress("rPFImbProjNC2_4", &rPFImbProjNC2_4_);
  trackTree_p->SetBranchAddress("rPFImbProjNC4_8", &rPFImbProjNC4_8_);
  trackTree_p->SetBranchAddress("rPFImbProjNC8_100", &rPFImbProjNC8_100_);

  //Corr. Tracks proj. onto PF, All, Cone, and NotCone

  trackTree_p->SetBranchAddress("rPFImbProjFCorr", &rPFImbProjFCorr_);
  trackTree_p->SetBranchAddress("rPFImbPerpFCorr", &rPFImbPerpFCorr_);
  trackTree_p->SetBranchAddress("rPFImbProj0_1Corr", &rPFImbProj0_1Corr_);
  trackTree_p->SetBranchAddress("rPFImbProj1_2Corr", &rPFImbProj1_2Corr_);
  trackTree_p->SetBranchAddress("rPFImbProj2_4Corr", &rPFImbProj2_4Corr_);
  trackTree_p->SetBranchAddress("rPFImbProj4_8Corr", &rPFImbProj4_8Corr_);
  trackTree_p->SetBranchAddress("rPFImbProj8_100Corr", &rPFImbProj8_100Corr_);

  trackTree_p->SetBranchAddress("rPFImbProjCFCorr", &rPFImbProjCFCorr_);
  trackTree_p->SetBranchAddress("rPFImbPerpCFCorr", &rPFImbPerpCFCorr_);
  trackTree_p->SetBranchAddress("rPFImbProjC0_1Corr", &rPFImbProjC0_1Corr_);
  trackTree_p->SetBranchAddress("rPFImbProjC1_2Corr", &rPFImbProjC1_2Corr_);
  trackTree_p->SetBranchAddress("rPFImbProjC2_4Corr", &rPFImbProjC2_4Corr_);
  trackTree_p->SetBranchAddress("rPFImbProjC4_8Corr", &rPFImbProjC4_8Corr_);
  trackTree_p->SetBranchAddress("rPFImbProjC8_100Corr", &rPFImbProjC8_100Corr_);
  trackTree_p->SetBranchAddress("rPFImbProjNCFCorr", &rPFImbProjNCFCorr_);
  trackTree_p->SetBranchAddress("rPFImbPerpNCFCorr", &rPFImbPerpNCFCorr_);
  trackTree_p->SetBranchAddress("rPFImbProjNC0_1Corr", &rPFImbProjNC0_1Corr_);
  trackTree_p->SetBranchAddress("rPFImbProjNC1_2Corr", &rPFImbProjNC1_2Corr_);
  trackTree_p->SetBranchAddress("rPFImbProjNC2_4Corr", &rPFImbProjNC2_4Corr_);
  trackTree_p->SetBranchAddress("rPFImbProjNC4_8Corr", &rPFImbProjNC4_8Corr_);
  trackTree_p->SetBranchAddress("rPFImbProjNC8_100Corr", &rPFImbProjNC8_100Corr_);

  //Tracks proj. onto Calo

  trackTree_p->SetBranchAddress("rCaloImbProjF", &rCaloImbProjF_);
  trackTree_p->SetBranchAddress("rCaloImbPerpF", &rCaloImbPerpF_);

  trackTree_p->SetBranchAddress("rCaloImbProj0_1", &rCaloImbProj0_1_);
  trackTree_p->SetBranchAddress("rCaloImbProj1_2", &rCaloImbProj1_2_);
  trackTree_p->SetBranchAddress("rCaloImbProj2_4", &rCaloImbProj2_4_);
  trackTree_p->SetBranchAddress("rCaloImbProj4_8", &rCaloImbProj4_8_);
  trackTree_p->SetBranchAddress("rCaloImbProj8_100", &rCaloImbProj8_100_);

  //Corr. Tracks proj. onto Calo

  trackTree_p->SetBranchAddress("rCaloImbProjFCorr", &rCaloImbProjFCorr_);
  trackTree_p->SetBranchAddress("rCaloImbPerpFCorr", &rCaloImbPerpFCorr_);

  trackTree_p->SetBranchAddress("rCaloImbProj0_1Corr", &rCaloImbProj0_1Corr_);
  trackTree_p->SetBranchAddress("rCaloImbProj1_2Corr", &rCaloImbProj1_2Corr_);
  trackTree_p->SetBranchAddress("rCaloImbProj2_4Corr", &rCaloImbProj2_4Corr_);
  trackTree_p->SetBranchAddress("rCaloImbProj4_8Corr", &rCaloImbProj4_8Corr_);
  trackTree_p->SetBranchAddress("rCaloImbProj8_100Corr", &rCaloImbProj8_100Corr_);

  //Tracks proj. onto Vs PF

  trackTree_p->SetBranchAddress("rVsPFImbProjF", &rVsPFImbProjF_);
  trackTree_p->SetBranchAddress("rVsPFImbPerpF", &rVsPFImbPerpF_);

  trackTree_p->SetBranchAddress("rVsPFImbProj0_1", &rVsPFImbProj0_1_);
  trackTree_p->SetBranchAddress("rVsPFImbProj1_2", &rVsPFImbProj1_2_);
  trackTree_p->SetBranchAddress("rVsPFImbProj2_4", &rVsPFImbProj2_4_);
  trackTree_p->SetBranchAddress("rVsPFImbProj4_8", &rVsPFImbProj4_8_);
  trackTree_p->SetBranchAddress("rVsPFImbProj8_100", &rVsPFImbProj8_100_);

  //Corr. Tracks proj. onto Vs PF

  trackTree_p->SetBranchAddress("rVsPFImbProjFCorr", &rVsPFImbProjFCorr_);
  trackTree_p->SetBranchAddress("rVsPFImbPerpFCorr", &rVsPFImbPerpFCorr_);

  trackTree_p->SetBranchAddress("rVsPFImbProj0_1Corr", &rVsPFImbProj0_1Corr_);
  trackTree_p->SetBranchAddress("rVsPFImbProj1_2Corr", &rVsPFImbProj1_2Corr_);
  trackTree_p->SetBranchAddress("rVsPFImbProj2_4Corr", &rVsPFImbProj2_4Corr_);
  trackTree_p->SetBranchAddress("rVsPFImbProj4_8Corr", &rVsPFImbProj4_8Corr_);
  trackTree_p->SetBranchAddress("rVsPFImbProj8_100Corr", &rVsPFImbProj8_100Corr_);

  //Tracks proj. onto Vs Calo

  trackTree_p->SetBranchAddress("rVsCaloImbProjF", &rVsCaloImbProjF_);
  trackTree_p->SetBranchAddress("rVsCaloImbPerpF", &rVsCaloImbPerpF_);

  trackTree_p->SetBranchAddress("rVsCaloImbProj0_1", &rVsCaloImbProj0_1_);
  trackTree_p->SetBranchAddress("rVsCaloImbProj1_2", &rVsCaloImbProj1_2_);
  trackTree_p->SetBranchAddress("rVsCaloImbProj2_4", &rVsCaloImbProj2_4_);
  trackTree_p->SetBranchAddress("rVsCaloImbProj4_8", &rVsCaloImbProj4_8_);
  trackTree_p->SetBranchAddress("rVsCaloImbProj8_100", &rVsCaloImbProj8_100_);

  //Corr. Tracks proj. onto Vs Calo

  trackTree_p->SetBranchAddress("rVsCaloImbProjFCorr", &rVsCaloImbProjFCorr_);
  trackTree_p->SetBranchAddress("rVsCaloImbPerpFCorr", &rVsCaloImbPerpFCorr_);

  trackTree_p->SetBranchAddress("rVsCaloImbProj0_1Corr", &rVsCaloImbProj0_1Corr_);
  trackTree_p->SetBranchAddress("rVsCaloImbProj1_2Corr", &rVsCaloImbProj1_2Corr_);
  trackTree_p->SetBranchAddress("rVsCaloImbProj2_4Corr", &rVsCaloImbProj2_4Corr_);
  trackTree_p->SetBranchAddress("rVsCaloImbProj4_8Corr", &rVsCaloImbProj4_8Corr_);
  trackTree_p->SetBranchAddress("rVsCaloImbProj8_100Corr", &rVsCaloImbProj8_100Corr_);

  if(montecarlo){
    //Track Tree Branches iff. Truth avail.

    //Tracks proj. onto Truth
    trackTree_p->SetBranchAddress("rTImbProjF", &rTImbProjF_);
    trackTree_p->SetBranchAddress("rTImbPerpF", &rTImbPerpF_);

    trackTree_p->SetBranchAddress("rTImbProj0_1", &rTImbProj0_1_);
    trackTree_p->SetBranchAddress("rTImbProj1_2", &rTImbProj1_2_);
    trackTree_p->SetBranchAddress("rTImbProj2_4", &rTImbProj2_4_);
    trackTree_p->SetBranchAddress("rTImbProj4_8", &rTImbProj4_8_);
    trackTree_p->SetBranchAddress("rTImbProj8_100", &rTImbProj8_100_);

    //Corr. Tracks proj. onto Truth

    trackTree_p->SetBranchAddress("rTImbProjFCorr", &rTImbProjFCorr_);
    trackTree_p->SetBranchAddress("rTImbPerpFCorr", &rTImbPerpFCorr_);

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
  jetTree_p->SetBranchAddress("recoVsPFSet", &recoVsPFSet_);
  jetTree_p->SetBranchAddress("recoCaloSet", &recoCaloSet_);
  jetTree_p->SetBranchAddress("recoVsCaloSet", &recoVsCaloSet_);

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

  jetTree_p->SetBranchAddress("VsPFLeadJtPt", &VsPFLeadJtPt_);
  jetTree_p->SetBranchAddress("VsPFLeadJtPhi", &VsPFLeadJtPhi_);
  jetTree_p->SetBranchAddress("VsPFLeadJtEta", &VsPFLeadJtEta_);
  jetTree_p->SetBranchAddress("VsPFSubLeadJtPt", &VsPFSubLeadJtPt_);
  jetTree_p->SetBranchAddress("VsPFSubLeadJtPhi", &VsPFSubLeadJtPhi_);
  jetTree_p->SetBranchAddress("VsPFSubLeadJtEta", &VsPFSubLeadJtEta_);
  jetTree_p->SetBranchAddress("VsPFJtAsymm", &VsPFJtAsymm_);

  jetTree_p->SetBranchAddress("VsCaloLeadJtPt", &VsCaloLeadJtPt_);
  jetTree_p->SetBranchAddress("VsCaloLeadJtPhi", &VsCaloLeadJtPhi_);
  jetTree_p->SetBranchAddress("VsCaloLeadJtEta", &VsCaloLeadJtEta_);
  jetTree_p->SetBranchAddress("VsCaloSubLeadJtPt", &VsCaloSubLeadJtPt_);
  jetTree_p->SetBranchAddress("VsCaloSubLeadJtPhi", &VsCaloSubLeadJtPhi_);
  jetTree_p->SetBranchAddress("VsCaloSubLeadJtEta", &VsCaloSubLeadJtEta_);
  jetTree_p->SetBranchAddress("VsCaloJtAsymm", &VsCaloJtAsymm_);

  if(montecarlo){
    //Gen Tree Branches

    genTree_p->SetBranchAddress("nGen", &nGen_);
    genTree_p->SetBranchAddress("genPt", &genPt_);
    genTree_p->SetBranchAddress("genPtPF", &genPtPF_);
    genTree_p->SetBranchAddress("genPtCalo", &genPtCalo_);
    genTree_p->SetBranchAddress("genPtT", &genPtT_);
    genTree_p->SetBranchAddress("genPtVsPF", &genPtVsPF_);
    genTree_p->SetBranchAddress("genPtVsCalo", &genPtVsCalo_);
    genTree_p->SetBranchAddress("genPhi", &genPhi_);
    genTree_p->SetBranchAddress("genEta", &genEta_);
    genTree_p->SetBranchAddress("genLeadDelPhi", &genLeadDelPhi_);

    //Gen. proj. onto Truth

    genTree_p->SetBranchAddress("gTImbProjF", &gTImbProjF_);
    genTree_p->SetBranchAddress("gTImbPerpF", &gTImbPerpF_);

    genTree_p->SetBranchAddress("gTImbProj0_1", &gTImbProj0_1_);
    genTree_p->SetBranchAddress("gTImbProj1_2", &gTImbProj1_2_);
    genTree_p->SetBranchAddress("gTImbProj2_4", &gTImbProj2_4_);
    genTree_p->SetBranchAddress("gTImbProj4_8", &gTImbProj4_8_);
    genTree_p->SetBranchAddress("gTImbProj8_100", &gTImbProj8_100_);

    //Gen. proj. onto PF

    genTree_p->SetBranchAddress("gPFImbProjF", &gPFImbProjF_);
    genTree_p->SetBranchAddress("gPFImbPerpF", &gPFImbPerpF_);

    genTree_p->SetBranchAddress("gPFImbProj0_1", &gPFImbProj0_1_);
    genTree_p->SetBranchAddress("gPFImbProj1_2", &gPFImbProj1_2_);
    genTree_p->SetBranchAddress("gPFImbProj2_4", &gPFImbProj2_4_);
    genTree_p->SetBranchAddress("gPFImbProj4_8", &gPFImbProj4_8_);
    genTree_p->SetBranchAddress("gPFImbProj8_100", &gPFImbProj8_100_);

    //Gen. proj. onto Calo

    genTree_p->SetBranchAddress("gCaloImbProjF", &gCaloImbProjF_);
    genTree_p->SetBranchAddress("gCaloImbPerpF", &gCaloImbPerpF_);

    genTree_p->SetBranchAddress("gCaloImbProj0_1", &gCaloImbProj0_1_);
    genTree_p->SetBranchAddress("gCaloImbProj1_2", &gCaloImbProj1_2_);
    genTree_p->SetBranchAddress("gCaloImbProj2_4", &gCaloImbProj2_4_);
    genTree_p->SetBranchAddress("gCaloImbProj4_8", &gCaloImbProj4_8_);
    genTree_p->SetBranchAddress("gCaloImbProj8_100", &gCaloImbProj8_100_);

    //Gen. proj. onto Vs PF

    genTree_p->SetBranchAddress("gVsPFImbProjF", &gVsPFImbProjF_);
    genTree_p->SetBranchAddress("gVsPFImbPerpF", &gVsPFImbPerpF_);

    genTree_p->SetBranchAddress("gVsPFImbProj0_1", &gVsPFImbProj0_1_);
    genTree_p->SetBranchAddress("gVsPFImbProj1_2", &gVsPFImbProj1_2_);
    genTree_p->SetBranchAddress("gVsPFImbProj2_4", &gVsPFImbProj2_4_);
    genTree_p->SetBranchAddress("gVsPFImbProj4_8", &gVsPFImbProj4_8_);
    genTree_p->SetBranchAddress("gVsPFImbProj8_100", &gVsPFImbProj8_100_);

    //Gen. proj. onto Vs Calo

    genTree_p->SetBranchAddress("gVsCaloImbProjF", &gVsCaloImbProjF_);
    genTree_p->SetBranchAddress("gVsCaloImbPerpF", &gVsCaloImbPerpF_);

    genTree_p->SetBranchAddress("gVsCaloImbProj0_1", &gVsCaloImbProj0_1_);
    genTree_p->SetBranchAddress("gVsCaloImbProj1_2", &gVsCaloImbProj1_2_);
    genTree_p->SetBranchAddress("gVsCaloImbProj2_4", &gVsCaloImbProj2_4_);
    genTree_p->SetBranchAddress("gVsCaloImbProj4_8", &gVsCaloImbProj4_8_);
    genTree_p->SetBranchAddress("gVsCaloImbProj8_100", &gVsCaloImbProj8_100_);
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
  std::cout << "DiJet Skim Init" << std::endl;

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
  recoVsPFSet_ = false;
  recoCaloSet_ = false;
  recoVsCaloSet_ = false;

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

  VsPFLeadJtPt_ = -10;
  VsPFSubLeadJtPt_ = -10;
  VsPFLeadJtPhi_ = -10;
  VsPFSubLeadJtPhi_ = -10;
  VsPFLeadJtEta_ = -10;
  VsPFSubLeadJtEta_ = -10;
  VsPFJtAsymm_ = -10;

  VsCaloLeadJtPt_ = -10;
  VsCaloSubLeadJtPt_ = -10;
  VsCaloLeadJtPhi_ = -10;
  VsCaloSubLeadJtPhi_ = -10;
  VsCaloLeadJtEta_ = -10;
  VsCaloSubLeadJtEta_ = -10;
  VsCaloJtAsymm_ = -10;
}

void InitProjPerp(bool montecarlo = false)
{
  //Tracks proj. onto PF, All, Cone, and NotCone

  rPFImbProjF_ = 0;
  rPFImbPerpF_ = 0;
  rPFImbProj0_1_ = 0;
  rPFImbProj1_2_ = 0;
  rPFImbProj2_4_ = 0;
  rPFImbProj4_8_ = 0;
  rPFImbProj8_100_ = 0;
  rPFImbProjCF_ = 0;
  rPFImbPerpCF_ = 0;
  rPFImbProjC0_1_ = 0;
  rPFImbProjC1_2_ = 0;
  rPFImbProjC2_4_ = 0;
  rPFImbProjC4_8_ = 0;
  rPFImbProjC8_100_ = 0;
  rPFImbProjNCF_ = 0;
  rPFImbPerpNCF_ = 0;
  rPFImbProjNC0_1_ = 0;
  rPFImbProjNC1_2_ = 0;
  rPFImbProjNC2_4_ = 0;
  rPFImbProjNC4_8_ = 0;
  rPFImbProjNC8_100_ = 0;

  //Corr. Tracks proj. onto PF, All, Cone, and NotCone

  rPFImbProjFCorr_ = 0;
  rPFImbPerpFCorr_ = 0;
  rPFImbProj0_1Corr_ = 0;
  rPFImbProj1_2Corr_ = 0;
  rPFImbProj2_4Corr_ = 0;
  rPFImbProj4_8Corr_ = 0;
  rPFImbProj8_100Corr_ = 0;
  rPFImbProjCFCorr_ = 0;
  rPFImbPerpCFCorr_ = 0;
  rPFImbProjC0_1Corr_ = 0;
  rPFImbProjC1_2Corr_ = 0;
  rPFImbProjC2_4Corr_ = 0;
  rPFImbProjC4_8Corr_ = 0;
  rPFImbProjC8_100Corr_ = 0;
  rPFImbProjNCFCorr_ = 0;
  rPFImbPerpNCFCorr_ = 0;
  rPFImbProjNC0_1Corr_ = 0;
  rPFImbProjNC1_2Corr_ = 0;
  rPFImbProjNC2_4Corr_ = 0;
  rPFImbProjNC4_8Corr_ = 0;
  rPFImbProjNC8_100Corr_ = 0;

  //Tracks proj. onto Calo

  rCaloImbProjF_ = 0;
  rCaloImbPerpF_ = 0;

  rCaloImbProj0_1_ = 0;
  rCaloImbProj1_2_ = 0;
  rCaloImbProj2_4_ = 0;
  rCaloImbProj4_8_ = 0;
  rCaloImbProj8_100_ = 0;

  //Corr. Tracks proj. onto Calo

  rCaloImbProjFCorr_ = 0;
  rCaloImbPerpFCorr_ = 0;

  rCaloImbProj0_1Corr_ = 0;
  rCaloImbProj1_2Corr_ = 0;
  rCaloImbProj2_4Corr_ = 0;
  rCaloImbProj4_8Corr_ = 0;
  rCaloImbProj8_100Corr_ = 0;

  //Tracks proj. onto Vs PF

  rVsPFImbProjF_ = 0;
  rVsPFImbPerpF_ = 0;

  rVsPFImbProj0_1_ = 0;
  rVsPFImbProj1_2_ = 0;
  rVsPFImbProj2_4_ = 0;
  rVsPFImbProj4_8_ = 0;
  rVsPFImbProj8_100_ = 0;

  //Corr. Tracks proj. onto Vs PF 

  rVsPFImbProjFCorr_ = 0;
  rVsPFImbPerpFCorr_ = 0;

  rVsPFImbProj0_1Corr_ = 0;
  rVsPFImbProj1_2Corr_ = 0;
  rVsPFImbProj2_4Corr_ = 0;
  rVsPFImbProj4_8Corr_ = 0;
  rVsPFImbProj8_100Corr_ = 0;

  //Tracks proj. onto Vs Calo

  rVsCaloImbProjF_ = 0;
  rVsCaloImbPerpF_ = 0;

  rVsCaloImbProj0_1_ = 0;
  rVsCaloImbProj1_2_ = 0;
  rVsCaloImbProj2_4_ = 0;
  rVsCaloImbProj4_8_ = 0;
  rVsCaloImbProj8_100_ = 0;

  //Corr. Tracks proj. onto Vs Calo

  rVsCaloImbProjFCorr_ = 0;
  rVsCaloImbPerpFCorr_ = 0;

  rVsCaloImbProj0_1Corr_ = 0;
  rVsCaloImbProj1_2Corr_ = 0;
  rVsCaloImbProj2_4Corr_ = 0;
  rVsCaloImbProj4_8Corr_ = 0;
  rVsCaloImbProj8_100Corr_ = 0;

  if(montecarlo){
    //Gen. proj. onto Truth

    gTImbProjF_ = 0;
    gTImbPerpF_ = 0;

    gTImbProj0_1_ = 0;
    gTImbProj1_2_ = 0;
    gTImbProj2_4_ = 0;
    gTImbProj4_8_ = 0;
    gTImbProj8_100_ = 0;

    //Gen. proj. onto PF

    gPFImbProjF_ = 0;
    gPFImbPerpF_ = 0;

    gPFImbProj0_1_ = 0;
    gPFImbProj1_2_ = 0;
    gPFImbProj2_4_ = 0;
    gPFImbProj4_8_ = 0;
    gPFImbProj8_100_ = 0;

    //Gen. proj. onto Calo

    gCaloImbProjF_ = 0;
    gCaloImbPerpF_ = 0;

    gCaloImbProj0_1_ = 0;
    gCaloImbProj1_2_ = 0;
    gCaloImbProj2_4_ = 0;
    gCaloImbProj4_8_ = 0;
    gCaloImbProj8_100_ = 0;

    //Gen. proj. onto Vs PF

    gVsPFImbProjF_ = 0;
    gVsPFImbPerpF_ = 0;

    gVsPFImbProj0_1_ = 0;
    gVsPFImbProj1_2_ = 0;
    gVsPFImbProj2_4_ = 0;
    gVsPFImbProj4_8_ = 0;
    gVsPFImbProj8_100_ = 0;

    //Gen. proj. onto Vs Calo

    gVsCaloImbProjF_ = 0;
    gVsCaloImbPerpF_ = 0;

    gVsCaloImbProj0_1_ = 0;
    gVsCaloImbProj1_2_ = 0;
    gVsCaloImbProj2_4_ = 0;
    gVsCaloImbProj4_8_ = 0;
    gVsCaloImbProj8_100_ = 0;

    //Tracks proj. onto Truth

    rTImbProjF_ = 0;
    rTImbPerpF_ = 0;

    rTImbProj0_1_ = 0;
    rTImbProj1_2_ = 0;
    rTImbProj2_4_ = 0;
    rTImbProj4_8_ = 0;
    rTImbProj8_100_ = 0;

    //Corr. Tracks proj. onto Truth

    rTImbProjFCorr_ = 0;
    rTImbPerpFCorr_ = 0;

    rTImbProj0_1Corr_ = 0;
    rTImbProj1_2Corr_ = 0;
    rTImbProj2_4Corr_ = 0;
    rTImbProj4_8Corr_ = 0;
    rTImbProj8_100Corr_ = 0;
  }
}

#endif
