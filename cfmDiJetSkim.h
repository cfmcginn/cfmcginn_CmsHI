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


enum AlgoType_PbPb {
  PuPF,
  PuCalo,
  VsPF,
  VsCalo,
  T,
  PuPFCorr,
  PuCaloCorr,
  VsPFCorr,
  VsCaloCorr,
  TCorr
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

const int MAXTRKS = 20000; //From SetupTrackTree.h

Int_t nTrk_;
Float_t trkPt_[MAXTRKS];
Float_t trkPtPuPF_[MAXTRKS];
Float_t trkPtPuCalo_[MAXTRKS];
Float_t trkPtT_[MAXTRKS];
Float_t trkPtVsPF_[MAXTRKS];
Float_t trkPtVsCalo_[MAXTRKS];
Float_t trkPhi_[MAXTRKS];
Float_t trkEta_[MAXTRKS];
Float_t trkPuPFLeadDelPhi_[MAXTRKS];
Float_t trkPuCaloLeadDelPhi_[MAXTRKS];
Float_t trkTLeadDelPhi_[MAXTRKS];
Float_t trkVsPFLeadDelPhi_[MAXTRKS];
Float_t trkVsCaloLeadDelPhi_[MAXTRKS];

Float_t trkRLeadPuPF_[MAXTRKS];
Float_t trkRSubLeadPuPF_[MAXTRKS];
Float_t trkRMinPuPF_[MAXTRKS];
Float_t trkPtCorrPuPF_[MAXTRKS];
Float_t trkPtFactPuPF_[MAXTRKS];

Float_t trkRLeadPuCalo_[MAXTRKS];
Float_t trkRSubLeadPuCalo_[MAXTRKS];
Float_t trkRMinPuCalo_[MAXTRKS];
Float_t trkPtCorrPuCalo_[MAXTRKS];
Float_t trkPtFactPuCalo_[MAXTRKS];

Float_t trkRLeadT_[MAXTRKS];
Float_t trkRSubLeadT_[MAXTRKS];
Float_t trkRMinT_[MAXTRKS];
Float_t trkPtCorrT_[MAXTRKS];
Float_t trkPtFactT_[MAXTRKS];

Float_t trkRLeadVsPF_[MAXTRKS];
Float_t trkRSubLeadVsPF_[MAXTRKS];
Float_t trkRMinVsPF_[MAXTRKS];
Float_t trkPtCorrVsPF_[MAXTRKS];
Float_t trkPtFactVsPF_[MAXTRKS];

Float_t trkRLeadVsCalo_[MAXTRKS];
Float_t trkRSubLeadVsCalo_[MAXTRKS];
Float_t trkRMinVsCalo_[MAXTRKS];
Float_t trkPtCorrVsCalo_[MAXTRKS];
Float_t trkPtFactVsCalo_[MAXTRKS];


//Tracks proj. onto Alg (enum ordered above, w/ corrected in back 5), All, Cone, and NotCone

Float_t rAlgImbProjAF_[10];
Float_t rAlgImbProjA0_1_[10];
Float_t rAlgImbProjA1_2_[10];
Float_t rAlgImbProjA2_4_[10];
Float_t rAlgImbProjA4_8_[10];
Float_t rAlgImbProjA8_100_[10];
Float_t rAlgImbProjACF_[10];
Float_t rAlgImbProjAC0_1_[10];
Float_t rAlgImbProjAC1_2_[10];
Float_t rAlgImbProjAC2_4_[10];
Float_t rAlgImbProjAC4_8_[10];
Float_t rAlgImbProjAC8_100_[10];
Float_t rAlgImbProjANCF_[10];
Float_t rAlgImbProjANC0_1_[10];
Float_t rAlgImbProjANC1_2_[10];
Float_t rAlgImbProjANC2_4_[10];
Float_t rAlgImbProjANC4_8_[10];
Float_t rAlgImbProjANC8_100_[10];
Float_t rAlgImbProjANCCutF_[10];
Float_t rAlgImbProjANCCut0_1_[10];
Float_t rAlgImbProjANCCut1_2_[10];
Float_t rAlgImbProjANCCut2_4_[10];
Float_t rAlgImbProjANCCut4_8_[10];
Float_t rAlgImbProjANCCut8_100_[10];

//DelR ProjA

Float_t rAlgImbProjA1CF_[10];
Float_t rAlgImbProjA1C0_1_[10];
Float_t rAlgImbProjA1C1_2_[10];
Float_t rAlgImbProjA1C2_4_[10];
Float_t rAlgImbProjA1C4_8_[10];
Float_t rAlgImbProjA1C8_100_[10];

Float_t rAlgImbProjA2CF_[10];
Float_t rAlgImbProjA2C0_1_[10];
Float_t rAlgImbProjA2C1_2_[10];
Float_t rAlgImbProjA2C2_4_[10];
Float_t rAlgImbProjA2C4_8_[10];
Float_t rAlgImbProjA2C8_100_[10];

Float_t rAlgImbProjA3CF_[10];
Float_t rAlgImbProjA3C0_1_[10];
Float_t rAlgImbProjA3C1_2_[10];
Float_t rAlgImbProjA3C2_4_[10];
Float_t rAlgImbProjA3C4_8_[10];
Float_t rAlgImbProjA3C8_100_[10];

Float_t rAlgImbProjA4CF_[10];
Float_t rAlgImbProjA4C0_1_[10];
Float_t rAlgImbProjA4C1_2_[10];
Float_t rAlgImbProjA4C2_4_[10];
Float_t rAlgImbProjA4C4_8_[10];
Float_t rAlgImbProjA4C8_100_[10];

Float_t rAlgImbProjA5CF_[10];
Float_t rAlgImbProjA5C0_1_[10];
Float_t rAlgImbProjA5C1_2_[10];
Float_t rAlgImbProjA5C2_4_[10];
Float_t rAlgImbProjA5C4_8_[10];
Float_t rAlgImbProjA5C8_100_[10];

Float_t rAlgImbProjA6CF_[10];
Float_t rAlgImbProjA6C0_1_[10];
Float_t rAlgImbProjA6C1_2_[10];
Float_t rAlgImbProjA6C2_4_[10];
Float_t rAlgImbProjA6C4_8_[10];
Float_t rAlgImbProjA6C8_100_[10];

Float_t rAlgImbProjA7CF_[10];
Float_t rAlgImbProjA7C0_1_[10];
Float_t rAlgImbProjA7C1_2_[10];
Float_t rAlgImbProjA7C2_4_[10];
Float_t rAlgImbProjA7C4_8_[10];
Float_t rAlgImbProjA7C8_100_[10];

Float_t rAlgImbProjA8CF_[10];
Float_t rAlgImbProjA8C0_1_[10];
Float_t rAlgImbProjA8C1_2_[10];
Float_t rAlgImbProjA8C2_4_[10];
Float_t rAlgImbProjA8C4_8_[10];
Float_t rAlgImbProjA8C8_100_[10];

Float_t rAlgImbProjA9CF_[10];
Float_t rAlgImbProjA9C0_1_[10];
Float_t rAlgImbProjA9C1_2_[10];
Float_t rAlgImbProjA9C2_4_[10];
Float_t rAlgImbProjA9C4_8_[10];
Float_t rAlgImbProjA9C8_100_[10];

Float_t rAlgImbProjA10CF_[10];
Float_t rAlgImbProjA10C0_1_[10];
Float_t rAlgImbProjA10C1_2_[10];
Float_t rAlgImbProjA10C2_4_[10];
Float_t rAlgImbProjA10C4_8_[10];
Float_t rAlgImbProjA10C8_100_[10];

//DelR tighter

Float_t rAlgImbProjA21CF_[10];
Float_t rAlgImbProjA21C0_1_[10];
Float_t rAlgImbProjA21C1_2_[10];
Float_t rAlgImbProjA21C2_4_[10];
Float_t rAlgImbProjA21C4_8_[10];
Float_t rAlgImbProjA21C8_100_[10];

Float_t rAlgImbProjA22CF_[10];
Float_t rAlgImbProjA22C0_1_[10];
Float_t rAlgImbProjA22C1_2_[10];
Float_t rAlgImbProjA22C2_4_[10];
Float_t rAlgImbProjA22C4_8_[10];
Float_t rAlgImbProjA22C8_100_[10];

Float_t rAlgImbProjA23CF_[10];
Float_t rAlgImbProjA23C0_1_[10];
Float_t rAlgImbProjA23C1_2_[10];
Float_t rAlgImbProjA23C2_4_[10];
Float_t rAlgImbProjA23C4_8_[10];
Float_t rAlgImbProjA23C8_100_[10];

Float_t rAlgImbProjA24CF_[10];
Float_t rAlgImbProjA24C0_1_[10];
Float_t rAlgImbProjA24C1_2_[10];
Float_t rAlgImbProjA24C2_4_[10];
Float_t rAlgImbProjA24C4_8_[10];
Float_t rAlgImbProjA24C8_100_[10];

//Jet Tree Variables

const int MAXJETS = 504; //From SetupJetTree.h

Int_t run_;
Int_t evt_;
Int_t lumi_;
Int_t hiBin_;
Int_t nref_;

Float_t hiEvtPlane_;
Float_t psin_;

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

//Event Set Bool array, [0] == PuPF, [1] == PuCalo, .etc according to enum

Bool_t eventSet_[5];
Float_t centWeight_[5];
Float_t centWeight_2pi3_[5];
Float_t centWeight_120_5pi6_[5];

Bool_t isQuarkJet_[5];
Bool_t isGluonJet_[5];

Float_t pthat_;
Float_t pthatWeight_;

//Jet Set, Array by algorithm, according to enum above

Float_t AlgLeadJtPt_[5];
Float_t AlgLeadJtPhi_[5];
Float_t AlgLeadJtEta_[5];
Float_t AlgSubLeadJtPt_[5];
Float_t AlgSubLeadJtPhi_[5];
Float_t AlgSubLeadJtEta_[5];
Float_t AlgThirdJtPt_[5];
Float_t AlgThirdJtPhi_[5];
Float_t AlgThirdJtEta_[5];
Float_t AlgFourthJtPt_[5];
Float_t AlgFourthJtPhi_[5];
Float_t AlgFourthJtEta_[5];

Float_t AlgJtAvePhi_[5];
Float_t AlgJtDelPhi_[5];
Float_t AlgJtAsymm_[5];

Float_t AlgLeadRefPt_[5];
Float_t AlgLeadRefPhi_[5];
Float_t AlgLeadRefEta_[5];
Float_t AlgSubLeadRefPt_[5];
Float_t AlgSubLeadRefPhi_[5];
Float_t AlgSubLeadRefEta_[5];
Float_t AlgThirdRefPt_[5];
Float_t AlgThirdRefPhi_[5];
Float_t AlgThirdRefEta_[5];
Float_t AlgFourthRefPt_[5];
Float_t AlgFourthRefPhi_[5];
Float_t AlgFourthRefEta_[5];

Float_t AlgRefDelPhi_[5];
Float_t AlgRefAsymm_[5];

//Gen Tree Variables

const int MAXGEN = 50000; //From SetupGenParticleTree.h

Int_t nGen_;
Float_t genPt_[MAXGEN];
Float_t genPtPuPF_[MAXGEN];
Float_t genPtPuCalo_[MAXGEN];
Float_t genPtT_[MAXGEN];
Float_t genPtVsPF_[MAXGEN];
Float_t genPtVsCalo_[MAXGEN];
Float_t genPhi_[MAXGEN];
Float_t genEta_[MAXGEN];
Float_t genLeadDelPhi_[MAXGEN];

//Gen. proj. onto Jets, ordered by algorithm according to enum, PuPF == [0], PuCalo == [1], etc.

Float_t gAlgImbProjAF_[5];
Float_t gAlgImbProjA0_1_[5];
Float_t gAlgImbProjA1_2_[5];
Float_t gAlgImbProjA2_4_[5];
Float_t gAlgImbProjA4_8_[5];
Float_t gAlgImbProjA8_100_[5];
Float_t gAlgImbProjACF_[5];
Float_t gAlgImbProjAC0_1_[5];
Float_t gAlgImbProjAC1_2_[5];
Float_t gAlgImbProjAC2_4_[5];
Float_t gAlgImbProjAC4_8_[5];
Float_t gAlgImbProjAC8_100_[5];
Float_t gAlgImbProjANCF_[5];
Float_t gAlgImbProjANC0_1_[5];
Float_t gAlgImbProjANC1_2_[5];
Float_t gAlgImbProjANC2_4_[5];
Float_t gAlgImbProjANC4_8_[5];
Float_t gAlgImbProjANC8_100_[5];
Float_t gAlgImbProjANCCutF_[5];
Float_t gAlgImbProjANCCut0_1_[5];
Float_t gAlgImbProjANCCut1_2_[5];
Float_t gAlgImbProjANCCut2_4_[5];
Float_t gAlgImbProjANCCut4_8_[5];
Float_t gAlgImbProjANCCut8_100_[5];

Float_t gAlgImbProjSigACF_[5];
Float_t gAlgImbProjSigAC0_1_[5];
Float_t gAlgImbProjSigAC1_2_[5];
Float_t gAlgImbProjSigAC2_4_[5];
Float_t gAlgImbProjSigAC4_8_[5];
Float_t gAlgImbProjSigAC8_100_[5];
Float_t gAlgImbProjSigANCF_[5];
Float_t gAlgImbProjSigANC0_1_[5];
Float_t gAlgImbProjSigANC1_2_[5];
Float_t gAlgImbProjSigANC2_4_[5];
Float_t gAlgImbProjSigANC4_8_[5];
Float_t gAlgImbProjSigANC8_100_[5];

Float_t gAlgImbProjUEACF_[5];
Float_t gAlgImbProjUEAC0_1_[5];
Float_t gAlgImbProjUEAC1_2_[5];
Float_t gAlgImbProjUEAC2_4_[5];
Float_t gAlgImbProjUEAC4_8_[5];
Float_t gAlgImbProjUEAC8_100_[5];
Float_t gAlgImbProjUEANCF_[5];
Float_t gAlgImbProjUEANC0_1_[5];
Float_t gAlgImbProjUEANC1_2_[5];
Float_t gAlgImbProjUEANC2_4_[5];
Float_t gAlgImbProjUEANC4_8_[5];
Float_t gAlgImbProjUEANC8_100_[5];

//truth delRs
//ProjA delR

Float_t gAlgImbProjA1CF_[5];
Float_t gAlgImbProjA1C0_1_[5];
Float_t gAlgImbProjA1C1_2_[5];
Float_t gAlgImbProjA1C2_4_[5];
Float_t gAlgImbProjA1C4_8_[5];
Float_t gAlgImbProjA1C8_100_[5];

Float_t gAlgImbProjA2CF_[5];
Float_t gAlgImbProjA2C0_1_[5];
Float_t gAlgImbProjA2C1_2_[5];
Float_t gAlgImbProjA2C2_4_[5];
Float_t gAlgImbProjA2C4_8_[5];
Float_t gAlgImbProjA2C8_100_[5];

Float_t gAlgImbProjA3CF_[5];
Float_t gAlgImbProjA3C0_1_[5];
Float_t gAlgImbProjA3C1_2_[5];
Float_t gAlgImbProjA3C2_4_[5];
Float_t gAlgImbProjA3C4_8_[5];
Float_t gAlgImbProjA3C8_100_[5];

Float_t gAlgImbProjA4CF_[5];
Float_t gAlgImbProjA4C0_1_[5];
Float_t gAlgImbProjA4C1_2_[5];
Float_t gAlgImbProjA4C2_4_[5];
Float_t gAlgImbProjA4C4_8_[5];
Float_t gAlgImbProjA4C8_100_[5];

Float_t gAlgImbProjA5CF_[5];
Float_t gAlgImbProjA5C0_1_[5];
Float_t gAlgImbProjA5C1_2_[5];
Float_t gAlgImbProjA5C2_4_[5];
Float_t gAlgImbProjA5C4_8_[5];
Float_t gAlgImbProjA5C8_100_[5];

Float_t gAlgImbProjA6CF_[5];
Float_t gAlgImbProjA6C0_1_[5];
Float_t gAlgImbProjA6C1_2_[5];
Float_t gAlgImbProjA6C2_4_[5];
Float_t gAlgImbProjA6C4_8_[5];
Float_t gAlgImbProjA6C8_100_[5];

Float_t gAlgImbProjA7CF_[5];
Float_t gAlgImbProjA7C0_1_[5];
Float_t gAlgImbProjA7C1_2_[5];
Float_t gAlgImbProjA7C2_4_[5];
Float_t gAlgImbProjA7C4_8_[5];
Float_t gAlgImbProjA7C8_100_[5];

Float_t gAlgImbProjA8CF_[5];
Float_t gAlgImbProjA8C0_1_[5];
Float_t gAlgImbProjA8C1_2_[5];
Float_t gAlgImbProjA8C2_4_[5];
Float_t gAlgImbProjA8C4_8_[5];
Float_t gAlgImbProjA8C8_100_[5];

Float_t gAlgImbProjA9CF_[5];
Float_t gAlgImbProjA9C0_1_[5];
Float_t gAlgImbProjA9C1_2_[5];
Float_t gAlgImbProjA9C2_4_[5];
Float_t gAlgImbProjA9C4_8_[5];
Float_t gAlgImbProjA9C8_100_[5];

Float_t gAlgImbProjA10CF_[5];
Float_t gAlgImbProjA10C0_1_[5];
Float_t gAlgImbProjA10C1_2_[5];
Float_t gAlgImbProjA10C2_4_[5];
Float_t gAlgImbProjA10C4_8_[5];
Float_t gAlgImbProjA10C8_100_[5];



//New DelRs

Float_t gAlgImbProjA21CF_[5];
Float_t gAlgImbProjA21C0_1_[5];
Float_t gAlgImbProjA21C1_2_[5];
Float_t gAlgImbProjA21C2_4_[5];
Float_t gAlgImbProjA21C4_8_[5];
Float_t gAlgImbProjA21C8_100_[5];

Float_t gAlgImbProjA22CF_[5];
Float_t gAlgImbProjA22C0_1_[5];
Float_t gAlgImbProjA22C1_2_[5];
Float_t gAlgImbProjA22C2_4_[5];
Float_t gAlgImbProjA22C4_8_[5];
Float_t gAlgImbProjA22C8_100_[5];

Float_t gAlgImbProjA23CF_[5];
Float_t gAlgImbProjA23C0_1_[5];
Float_t gAlgImbProjA23C1_2_[5];
Float_t gAlgImbProjA23C2_4_[5];
Float_t gAlgImbProjA23C4_8_[5];
Float_t gAlgImbProjA23C8_100_[5];

Float_t gAlgImbProjA24CF_[5];
Float_t gAlgImbProjA24C0_1_[5];
Float_t gAlgImbProjA24C1_2_[5];
Float_t gAlgImbProjA24C2_4_[5];
Float_t gAlgImbProjA24C4_8_[5];
Float_t gAlgImbProjA24C8_100_[5];

void SetBranches(bool montecarlo)
{
  //Track Tree Branches

  std::cout << "Branches Set" << std::endl;

  /*  
  trackTree_p->Branch("nTrk", &nTrk_, "nTrk/I");
  trackTree_p->Branch("trkPt", &trkPt_, "trkPt[nTrk]/F");
  trackTree_p->Branch("trkPtPuPF", &trkPtPuPF_, "trkPtPuPF[nTrk]/F");
  trackTree_p->Branch("trkPtPuCalo", &trkPtPuCalo_, "trkPtPuCalo[nTrk]/F");
  trackTree_p->Branch("trkPtVsPF", &trkPtVsPF_, "trkPtVsPF[nTrk]/F");
  trackTree_p->Branch("trkPtVsCalo", &trkPtVsCalo_, "trkPtVsCalo[nTrk]/F");

  if(montecarlo)
    trackTree_p->Branch("trkPtT", &trkPtT_, "trkPtT[nTrk]/F");

  trackTree_p->Branch("trkPhi", &trkPhi_, "trkPhi[nTrk]/F");
  trackTree_p->Branch("trkEta", &trkEta_, "trkEta[nTrk]/F");
  trackTree_p->Branch("trkPuPFLeadDelPhi", &trkPuPFLeadDelPhi_, "trkPuPFLeadDelPhi[nTrk]/F");
  trackTree_p->Branch("trkPuCaloLeadDelPhi", &trkPuCaloLeadDelPhi_, "trkPuCaloLeadDelPhi[nTrk]/F");
  trackTree_p->Branch("trkTLeadDelPhi", &trkTLeadDelPhi_, "trkTLeadDelPhi[nTrk]/F");
  trackTree_p->Branch("trkVsPFLeadDelPhi", &trkVsPFLeadDelPhi_, "trkVsPFLeadDelPhi[nTrk]/F");
  trackTree_p->Branch("trkVsCaloLeadDelPhi", &trkVsCaloLeadDelPhi_, "trkVsCaloLeadDelPhi[nTrk]/F");

  trackTree_p->Branch("trkRLeadPuPF", &trkRLeadPuPF_, "trkRLeadPuPF[nTrk]/F");
  trackTree_p->Branch("trkRSubLeadPuPF", &trkRSubLeadPuPF_, "trkRSubLeadPuPF[nTrk]/F");
  trackTree_p->Branch("trkRMinPuPF", &trkRMinPuPF_, "trkRMinPuPF[nTrk]/F");
  trackTree_p->Branch("trkPtCorrPuPF", &trkPtCorrPuPF_, "trkPtCorrPuPF[nTrk]/F");
  trackTree_p->Branch("trkPtFactPuPF", &trkPtFactPuPF_, "trkPtFactPuPF[nTrk]/F");

  trackTree_p->Branch("trkRLeadPuCalo", &trkRLeadPuCalo_, "trkRLeadPuCalo[nTrk]/F");
  trackTree_p->Branch("trkRSubLeadPuCalo", &trkRSubLeadPuCalo_, "trkRSubLeadPuCalo[nTrk]/F");
  trackTree_p->Branch("trkRMinPuCalo", &trkRMinPuCalo_, "trkRMinPuCalo[nTrk]/F");
  trackTree_p->Branch("trkPtCorrPuCalo", &trkPtCorrPuCalo_, "trkPtCorrPuCalo[nTrk]/F");
  trackTree_p->Branch("trkPtFactPuCalo", &trkPtFactPuCalo_, "trkPtFactPuCalo[nTrk]/F");

  if(montecarlo){
    trackTree_p->Branch("trkRLeadT", &trkRLeadT_, "trkRLeadT[nTrk]/F");
    trackTree_p->Branch("trkRSubLeadT", &trkRSubLeadT_, "trkRSubLeadT[nTrk]/F");
    trackTree_p->Branch("trkRMinT", &trkRMinT_, "trkRMinT[nTrk]/F");
    trackTree_p->Branch("trkPtCorrT", &trkPtCorrT_, "trkPtCorrT[nTrk]/F");
    trackTree_p->Branch("trkPtFactT", &trkPtFactT_, "trkPtFactT[nTrk]/F");
  }

  trackTree_p->Branch("trkRLeadVsPF", &trkRLeadVsPF_, "trkRLeadVsPF[nTrk]/F");
  trackTree_p->Branch("trkRSubLeadVsPF", &trkRSubLeadVsPF_, "trkRSubLeadVsPF[nTrk]/F");
  trackTree_p->Branch("trkRMinVsPF", &trkRMinVsPF_, "trkRMinVsPF[nTrk]/F");
  trackTree_p->Branch("trkPtCorrVsPF", &trkPtCorrVsPF_, "trkPtCorrVsPF[nTrk]/F");
  trackTree_p->Branch("trkPtFactVsPF", &trkPtFactVsPF_, "trkPtFactVsPF[nTrk]/F");

  trackTree_p->Branch("trkRLeadVsCalo", &trkRLeadVsCalo_, "trkRLeadVsCalo[nTrk]/F");
  trackTree_p->Branch("trkRSubLeadVsCalo", &trkRSubLeadVsCalo_, "trkRSubLeadVsCalo[nTrk]/F");
  trackTree_p->Branch("trkRMinVsCalo", &trkRMinVsCalo_, "trkRMinVsCalo[nTrk]/F");
  trackTree_p->Branch("trkPtCorrVsCalo", &trkPtCorrVsCalo_, "trkPtCorrVsCalo[nTrk]/F");
  trackTree_p->Branch("trkPtFactVsCalo", &trkPtFactVsCalo_, "trkPtFactVsCalo[nTrk]/F");  
  */  
  //Tracks proj. onto Alg, ordered according to enum above, All, Cone, and NotCone

  trackTree_p->Branch("rAlgImbProjAF", &rAlgImbProjAF_, "rAlgImbProjAF[10]/F");
  trackTree_p->Branch("rAlgImbProjA0_1", &rAlgImbProjA0_1_, "rAlgImbProjA0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA1_2", &rAlgImbProjA1_2_, "rAlgImbProjA1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA2_4", &rAlgImbProjA2_4_, "rAlgImbProjA2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA4_8", &rAlgImbProjA4_8_, "rAlgImbProjA4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA8_100", &rAlgImbProjA8_100_, "rAlgImbProjA8_100[10]/F");
  trackTree_p->Branch("rAlgImbProjACF", &rAlgImbProjACF_, "rAlgImbProjACF[10]/F");
  trackTree_p->Branch("rAlgImbProjAC0_1", &rAlgImbProjAC0_1_, "rAlgImbProjAC0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjAC1_2", &rAlgImbProjAC1_2_, "rAlgImbProjAC1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjAC2_4", &rAlgImbProjAC2_4_, "rAlgImbProjAC2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjAC4_8", &rAlgImbProjAC4_8_, "rAlgImbProjAC4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjAC8_100", &rAlgImbProjAC8_100_, "rAlgImbProjAC8_100[10]/F");
  trackTree_p->Branch("rAlgImbProjANCF", &rAlgImbProjANCF_, "rAlgImbProjANCF[10]/F");
  trackTree_p->Branch("rAlgImbProjANC0_1", &rAlgImbProjANC0_1_, "rAlgImbProjANC0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjANC1_2", &rAlgImbProjANC1_2_, "rAlgImbProjANC1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjANC2_4", &rAlgImbProjANC2_4_, "rAlgImbProjANC2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjANC4_8", &rAlgImbProjANC4_8_, "rAlgImbProjANC4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjANC8_100", &rAlgImbProjANC8_100_, "rAlgImbProjANC8_100[10]/F");
  trackTree_p->Branch("rAlgImbProjANCCutF", &rAlgImbProjANCCutF_, "rAlgImbProjANCCutF[10]/F");
  trackTree_p->Branch("rAlgImbProjANCCut0_1", &rAlgImbProjANCCut0_1_, "rAlgImbProjANCCut0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjANCCut1_2", &rAlgImbProjANCCut1_2_, "rAlgImbProjANCCut1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjANCCut2_4", &rAlgImbProjANCCut2_4_, "rAlgImbProjANCCut2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjANCCut4_8", &rAlgImbProjANCCut4_8_, "rAlgImbProjANCCut4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjANCCut8_100", &rAlgImbProjANCCut8_100_, "rAlgImbProjANCCut8_100[10]/F");

  //ProjA DelRs

  trackTree_p->Branch("rAlgImbProjA1CF", &rAlgImbProjA1CF_, "rAlgImbProjA1CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA1C0_1", &rAlgImbProjA1C0_1_, "rAlgImbProjA1C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA1C1_2", &rAlgImbProjA1C1_2_, "rAlgImbProjA1C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA1C2_4", &rAlgImbProjA1C2_4_, "rAlgImbProjA1C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA1C4_8", &rAlgImbProjA1C4_8_, "rAlgImbProjA1C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA1C8_100", &rAlgImbProjA1C8_100_, "rAlgImbProjA1C8_100[10]/F");

  trackTree_p->Branch("rAlgImbProjA2CF", &rAlgImbProjA2CF_, "rAlgImbProjA2CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA2C0_1", &rAlgImbProjA2C0_1_, "rAlgImbProjA2C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA2C1_2", &rAlgImbProjA2C1_2_, "rAlgImbProjA2C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA2C2_4", &rAlgImbProjA2C2_4_, "rAlgImbProjA2C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA2C4_8", &rAlgImbProjA2C4_8_, "rAlgImbProjA2C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA2C8_100", &rAlgImbProjA2C8_100_, "rAlgImbProjA2C8_100[10]/F");

  trackTree_p->Branch("rAlgImbProjA3CF", &rAlgImbProjA3CF_, "rAlgImbProjA3CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA3C0_1", &rAlgImbProjA3C0_1_, "rAlgImbProjA3C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA3C1_2", &rAlgImbProjA3C1_2_, "rAlgImbProjA3C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA3C2_4", &rAlgImbProjA3C2_4_, "rAlgImbProjA3C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA3C4_8", &rAlgImbProjA3C4_8_, "rAlgImbProjA3C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA3C8_100", &rAlgImbProjA3C8_100_, "rAlgImbProjA3C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProjA4CF", &rAlgImbProjA4CF_, "rAlgImbProjA4CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA4C0_1", &rAlgImbProjA4C0_1_, "rAlgImbProjA4C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA4C1_2", &rAlgImbProjA4C1_2_, "rAlgImbProjA4C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA4C2_4", &rAlgImbProjA4C2_4_, "rAlgImbProjA4C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA4C4_8", &rAlgImbProjA4C4_8_, "rAlgImbProjA4C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA4C8_100", &rAlgImbProjA4C8_100_, "rAlgImbProjA4C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProjA5CF", &rAlgImbProjA5CF_, "rAlgImbProjA5CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA5C0_1", &rAlgImbProjA5C0_1_, "rAlgImbProjA5C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA5C1_2", &rAlgImbProjA5C1_2_, "rAlgImbProjA5C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA5C2_4", &rAlgImbProjA5C2_4_, "rAlgImbProjA5C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA5C4_8", &rAlgImbProjA5C4_8_, "rAlgImbProjA5C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA5C8_100", &rAlgImbProjA5C8_100_, "rAlgImbProjA5C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProjA6CF", &rAlgImbProjA6CF_, "rAlgImbProjA6CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA6C0_1", &rAlgImbProjA6C0_1_, "rAlgImbProjA6C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA6C1_2", &rAlgImbProjA6C1_2_, "rAlgImbProjA6C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA6C2_4", &rAlgImbProjA6C2_4_, "rAlgImbProjA6C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA6C4_8", &rAlgImbProjA6C4_8_, "rAlgImbProjA6C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA6C8_100", &rAlgImbProjA6C8_100_, "rAlgImbProjA6C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProjA7CF", &rAlgImbProjA7CF_, "rAlgImbProjA7CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA7C0_1", &rAlgImbProjA7C0_1_, "rAlgImbProjA7C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA7C1_2", &rAlgImbProjA7C1_2_, "rAlgImbProjA7C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA7C2_4", &rAlgImbProjA7C2_4_, "rAlgImbProjA7C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA7C4_8", &rAlgImbProjA7C4_8_, "rAlgImbProjA7C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA7C8_100", &rAlgImbProjA7C8_100_, "rAlgImbProjA7C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProjA8CF", &rAlgImbProjA8CF_, "rAlgImbProjA8CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA8C0_1", &rAlgImbProjA8C0_1_, "rAlgImbProjA8C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA8C1_2", &rAlgImbProjA8C1_2_, "rAlgImbProjA8C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA8C2_4", &rAlgImbProjA8C2_4_, "rAlgImbProjA8C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA8C4_8", &rAlgImbProjA8C4_8_, "rAlgImbProjA8C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA8C8_100", &rAlgImbProjA8C8_100_, "rAlgImbProjA8C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProjA9CF", &rAlgImbProjA9CF_, "rAlgImbProjA9CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA9C0_1", &rAlgImbProjA9C0_1_, "rAlgImbProjA9C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA9C1_2", &rAlgImbProjA9C1_2_, "rAlgImbProjA9C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA9C2_4", &rAlgImbProjA9C2_4_, "rAlgImbProjA9C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA9C4_8", &rAlgImbProjA9C4_8_, "rAlgImbProjA9C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA9C8_100", &rAlgImbProjA9C8_100_, "rAlgImbProjA9C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProjA10CF", &rAlgImbProjA10CF_, "rAlgImbProjA10CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA10C0_1", &rAlgImbProjA10C0_1_, "rAlgImbProjA10C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA10C1_2", &rAlgImbProjA10C1_2_, "rAlgImbProjA10C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA10C2_4", &rAlgImbProjA10C2_4_, "rAlgImbProjA10C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA10C4_8", &rAlgImbProjA10C4_8_, "rAlgImbProjA10C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA10C8_100", &rAlgImbProjA10C8_100_, "rAlgImbProjA10C8_100[10]/F");  


  //Delrs

  trackTree_p->Branch("rAlgImbProjA21CF", &rAlgImbProjA21CF_, "rAlgImbProjA21CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA21C0_1", &rAlgImbProjA21C0_1_, "rAlgImbProjA21C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA21C1_2", &rAlgImbProjA21C1_2_, "rAlgImbProjA21C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA21C2_4", &rAlgImbProjA21C2_4_, "rAlgImbProjA21C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA21C4_8", &rAlgImbProjA21C4_8_, "rAlgImbProjA21C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA21C8_100", &rAlgImbProjA21C8_100_, "rAlgImbProjA21C8_100[10]/F");

  trackTree_p->Branch("rAlgImbProjA22CF", &rAlgImbProjA22CF_, "rAlgImbProjA22CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA22C0_1", &rAlgImbProjA22C0_1_, "rAlgImbProjA22C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA22C1_2", &rAlgImbProjA22C1_2_, "rAlgImbProjA22C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA22C2_4", &rAlgImbProjA22C2_4_, "rAlgImbProjA22C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA22C4_8", &rAlgImbProjA22C4_8_, "rAlgImbProjA22C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA22C8_100", &rAlgImbProjA22C8_100_, "rAlgImbProjA22C8_100[10]/F");

  trackTree_p->Branch("rAlgImbProjA23CF", &rAlgImbProjA23CF_, "rAlgImbProjA23CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA23C0_1", &rAlgImbProjA23C0_1_, "rAlgImbProjA23C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA23C1_2", &rAlgImbProjA23C1_2_, "rAlgImbProjA23C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA23C2_4", &rAlgImbProjA23C2_4_, "rAlgImbProjA23C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA23C4_8", &rAlgImbProjA23C4_8_, "rAlgImbProjA23C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA23C8_100", &rAlgImbProjA23C8_100_, "rAlgImbProjA23C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProjA24CF", &rAlgImbProjA24CF_, "rAlgImbProjA24CF[10]/F");
  trackTree_p->Branch("rAlgImbProjA24C0_1", &rAlgImbProjA24C0_1_, "rAlgImbProjA24C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjA24C1_2", &rAlgImbProjA24C1_2_, "rAlgImbProjA24C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjA24C2_4", &rAlgImbProjA24C2_4_, "rAlgImbProjA24C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjA24C4_8", &rAlgImbProjA24C4_8_, "rAlgImbProjA24C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjA24C8_100", &rAlgImbProjA24C8_100_, "rAlgImbProjA24C8_100[10]/F");  


  //Jet Tree Branches

  jetTree_p->Branch("run", &run_, "run/I");
  jetTree_p->Branch("evt", &evt_, "evt/I");
  jetTree_p->Branch("lumi", &lumi_, "lumi/I");
  jetTree_p->Branch("hiBin", &hiBin_, "hiBin/I");

  jetTree_p->Branch("hiEvtPlane", &hiEvtPlane_, "hiEvtPlane/F");
  jetTree_p->Branch("psin", &psin_, "psin/F");

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

  jetTree_p->Branch("eventSet", &eventSet_, "eventSet[5]/O");
  jetTree_p->Branch("centWeight", &centWeight_, "centWeight[5]/F");
  jetTree_p->Branch("centWeight_2pi3", &centWeight_2pi3_, "centWeight_2pi3[5]/F");
  jetTree_p->Branch("centWeight_120_5pi6", &centWeight_120_5pi6_, "centWeight_120_5pi6[5]/F");

  jetTree_p->Branch("AlgLeadJtPt", &AlgLeadJtPt_, "AlgLeadJtPt[5]/F");
  jetTree_p->Branch("AlgLeadJtPhi", &AlgLeadJtPhi_, "AlgLeadJtPhi[5]/F");
  jetTree_p->Branch("AlgLeadJtEta", &AlgLeadJtEta_, "AlgLeadJtEta[5]/F");
  jetTree_p->Branch("AlgSubLeadJtPt", &AlgSubLeadJtPt_, "AlgSubLeadJtPt[5]/F");
  jetTree_p->Branch("AlgSubLeadJtPhi", &AlgSubLeadJtPhi_, "AlgSubLeadJtPhi[5]/F");
  jetTree_p->Branch("AlgSubLeadJtEta", &AlgSubLeadJtEta_, "AlgSubLeadJtEta[5]/F");
  jetTree_p->Branch("AlgThirdJtPt", &AlgThirdJtPt_, "AlgThirdJtPt[5]/F");
  jetTree_p->Branch("AlgThirdJtPhi", &AlgThirdJtPhi_, "AlgThirdJtPhi[5]/F");
  jetTree_p->Branch("AlgThirdJtEta", &AlgThirdJtEta_, "AlgThirdJtEta[5]/F");
  jetTree_p->Branch("AlgFourthJtPt", &AlgFourthJtPt_, "AlgFourthJtPt[5]/F");
  jetTree_p->Branch("AlgFourthJtPhi", &AlgFourthJtPhi_, "AlgFourthJtPhi[5]/F");
  jetTree_p->Branch("AlgFourthJtEta", &AlgFourthJtEta_, "AlgFourthJtEta[5]/F");
  jetTree_p->Branch("AlgJtAvePhi", &AlgJtAvePhi_, "AlgJtAvePhi[5]/F");
  jetTree_p->Branch("AlgJtDelPhi", &AlgJtDelPhi_, "AlgJtDelPhi[5]/F");
  jetTree_p->Branch("AlgJtAsymm", &AlgJtAsymm_, "AlgJtAsymm[5]/F");

  if(montecarlo){
    jetTree_p->Branch("isQuarkJet", &isQuarkJet_, "isQuarkJet[5]/O");
    jetTree_p->Branch("isGluonJet", &isGluonJet_, "isGluonJet[5]/O");

    jetTree_p->Branch("pthatWeight", &pthatWeight_, "pthatWeight/F");
    jetTree_p->Branch("pthat", &pthat_, "pthat/F");

    //refpt for jets immediately above
    jetTree_p->Branch("AlgLeadRefPt", &AlgLeadRefPt_, "AlgLeadRefPt[5]/F");
    jetTree_p->Branch("AlgLeadRefPhi", &AlgLeadRefPhi_, "AlgLeadRefPhi[5]/F");
    jetTree_p->Branch("AlgLeadRefEta", &AlgLeadRefEta_, "AlgLeadRefEta[5]/F");
    jetTree_p->Branch("AlgSubLeadRefPt", &AlgSubLeadRefPt_, "AlgSubLeadRefPt[5]/F");
    jetTree_p->Branch("AlgSubLeadRefPhi", &AlgSubLeadRefPhi_, "AlgSubLeadRefPhi[5]/F");
    jetTree_p->Branch("AlgSubLeadRefEta", &AlgSubLeadRefEta_, "AlgSubLeadRefEta[5]/F");
    jetTree_p->Branch("AlgThirdRefPt", &AlgThirdRefPt_, "AlgThirdRefPt[5]/F");
    jetTree_p->Branch("AlgThirdRefPhi", &AlgThirdRefPhi_, "AlgThirdRefPhi[5]/F");
    jetTree_p->Branch("AlgThirdRefEta", &AlgThirdRefEta_, "AlgThirdRefEta[5]/F");
    jetTree_p->Branch("AlgFourthRefPt", &AlgFourthRefPt_, "AlgFourthRefPt[5]/F");
    jetTree_p->Branch("AlgFourthRefPhi", &AlgFourthRefPhi_, "AlgFourthRefPhi[5]/F");
    jetTree_p->Branch("AlgFourthRefEta", &AlgFourthRefEta_, "AlgFourthRefEta[5]/F");

    jetTree_p->Branch("AlgRefDelPhi", &AlgRefDelPhi_, "AlgRefDelPhi[5]/F");
    jetTree_p->Branch("AlgRefAsymm", &AlgRefAsymm_, "AlgRefAsymm[5]/F");

    //Gen Tree Branches

    /*
    genTree_p->Branch("nGen", &nGen_, "nGen/I");
    genTree_p->Branch("genPt", &genPt_, "genPt[nGen]/F");
    genTree_p->Branch("genPtPuPF", &genPtPuPF_, "genPtPuPF[nGen]/F");
    genTree_p->Branch("genPtPuCalo", &genPtPuCalo_, "genPtPuCalo[nGen]/F");
    genTree_p->Branch("genPtT", &genPtT_, "genPtT[nGen]/F");
    genTree_p->Branch("genPtVsPF", &genPtVsPF_, "genPtVsPF[nGen]/F");
    genTree_p->Branch("genPtVsCalo", &genPtVsCalo_, "genPtVsCalo[nGen]/F");
    genTree_p->Branch("genPhi", &genPhi_, "genPhi[nGen]/F");
    genTree_p->Branch("genEta", &genEta_, "genEta[nGen]/F");
    genTree_p->Branch("genLeadDelPhi", &genLeadDelPhi_, "genLeadDelPhi[nGen]/F");
    */    

    //Gen. proj. onto jetAlg, array ordered according to enum

    genTree_p->Branch("gAlgImbProjAF", &gAlgImbProjAF_, "gAlgImbProjAF[5]/F");
    genTree_p->Branch("gAlgImbProjA0_1", &gAlgImbProjA0_1_, "gAlgImbProjA0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA1_2", &gAlgImbProjA1_2_, "gAlgImbProjA1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA2_4", &gAlgImbProjA2_4_, "gAlgImbProjA2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA4_8", &gAlgImbProjA4_8_, "gAlgImbProjA4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA8_100", &gAlgImbProjA8_100_, "gAlgImbProjA8_100[5]/F");
    genTree_p->Branch("gAlgImbProjACF", &gAlgImbProjACF_, "gAlgImbProjACF[5]/F");
    genTree_p->Branch("gAlgImbProjAC0_1", &gAlgImbProjAC0_1_, "gAlgImbProjAC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjAC1_2", &gAlgImbProjAC1_2_, "gAlgImbProjAC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjAC2_4", &gAlgImbProjAC2_4_, "gAlgImbProjAC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjAC4_8", &gAlgImbProjAC4_8_, "gAlgImbProjAC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjAC8_100", &gAlgImbProjAC8_100_, "gAlgImbProjAC8_100[5]/F");
    genTree_p->Branch("gAlgImbProjANCF", &gAlgImbProjANCF_, "gAlgImbProjANCF[5]/F");
    genTree_p->Branch("gAlgImbProjANC0_1", &gAlgImbProjANC0_1_, "gAlgImbProjANC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjANC1_2", &gAlgImbProjANC1_2_, "gAlgImbProjANC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjANC2_4", &gAlgImbProjANC2_4_, "gAlgImbProjANC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjANC4_8", &gAlgImbProjANC4_8_, "gAlgImbProjANC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjANC8_100", &gAlgImbProjANC8_100_, "gAlgImbProjANC8_100[5]/F");
    genTree_p->Branch("gAlgImbProjANCCutF", &gAlgImbProjANCCutF_, "gAlgImbProjANCCutF[5]/F");
    genTree_p->Branch("gAlgImbProjANCCut0_1", &gAlgImbProjANCCut0_1_, "gAlgImbProjANCCut0_1[5]/F");
    genTree_p->Branch("gAlgImbProjANCCut1_2", &gAlgImbProjANCCut1_2_, "gAlgImbProjANCCut1_2[5]/F");
    genTree_p->Branch("gAlgImbProjANCCut2_4", &gAlgImbProjANCCut2_4_, "gAlgImbProjANCCut2_4[5]/F");
    genTree_p->Branch("gAlgImbProjANCCut4_8", &gAlgImbProjANCCut4_8_, "gAlgImbProjANCCut4_8[5]/F");
    genTree_p->Branch("gAlgImbProjANCCut8_100", &gAlgImbProjANCCut8_100_, "gAlgImbProjANCCut8_100[5]/F");

    genTree_p->Branch("gAlgImbProjSigACF", &gAlgImbProjSigACF_, "gAlgImbProjSigACF[5]/F");
    genTree_p->Branch("gAlgImbProjSigAC0_1", &gAlgImbProjSigAC0_1_, "gAlgImbProjSigAC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjSigAC1_2", &gAlgImbProjSigAC1_2_, "gAlgImbProjSigAC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjSigAC2_4", &gAlgImbProjSigAC2_4_, "gAlgImbProjSigAC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjSigAC4_8", &gAlgImbProjSigAC4_8_, "gAlgImbProjSigAC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjSigAC8_100", &gAlgImbProjSigAC8_100_, "gAlgImbProjSigAC8_100[5]/F");
    genTree_p->Branch("gAlgImbProjSigANCF", &gAlgImbProjSigANCF_, "gAlgImbProjSigANCF[5]/F");
    genTree_p->Branch("gAlgImbProjSigANC0_1", &gAlgImbProjSigANC0_1_, "gAlgImbProjSigANC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjSigANC1_2", &gAlgImbProjSigANC1_2_, "gAlgImbProjSigANC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjSigANC2_4", &gAlgImbProjSigANC2_4_, "gAlgImbProjSigANC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjSigANC4_8", &gAlgImbProjSigANC4_8_, "gAlgImbProjSigANC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjSigANC8_100", &gAlgImbProjSigANC8_100_, "gAlgImbProjSigANC8_100[5]/F");

    genTree_p->Branch("gAlgImbProjUEACF", &gAlgImbProjUEACF_, "gAlgImbProjUEACF[5]/F");
    genTree_p->Branch("gAlgImbProjUEAC0_1", &gAlgImbProjUEAC0_1_, "gAlgImbProjUEAC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjUEAC1_2", &gAlgImbProjUEAC1_2_, "gAlgImbProjUEAC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjUEAC2_4", &gAlgImbProjUEAC2_4_, "gAlgImbProjUEAC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjUEAC4_8", &gAlgImbProjUEAC4_8_, "gAlgImbProjUEAC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjUEAC8_100", &gAlgImbProjUEAC8_100_, "gAlgImbProjUEAC8_100[5]/F");
    genTree_p->Branch("gAlgImbProjUEANCF", &gAlgImbProjUEANCF_, "gAlgImbProjUEANCF[5]/F");
    genTree_p->Branch("gAlgImbProjUEANC0_1", &gAlgImbProjUEANC0_1_, "gAlgImbProjUEANC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjUEANC1_2", &gAlgImbProjUEANC1_2_, "gAlgImbProjUEANC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjUEANC2_4", &gAlgImbProjUEANC2_4_, "gAlgImbProjUEANC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjUEANC4_8", &gAlgImbProjUEANC4_8_, "gAlgImbProjUEANC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjUEANC8_100", &gAlgImbProjUEANC8_100_, "gAlgImbProjUEANC8_100[5]/F");

    //Truth Del Rs, Proj
    //Proj A's

    genTree_p->Branch("gAlgImbProjA1CF", &gAlgImbProjA1CF_, "gAlgImbProjA1CF[5]/F");
    genTree_p->Branch("gAlgImbProjA1C0_1", &gAlgImbProjA1C0_1_, "gAlgImbProjA1C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA1C1_2", &gAlgImbProjA1C1_2_, "gAlgImbProjA1C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA1C2_4", &gAlgImbProjA1C2_4_, "gAlgImbProjA1C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA1C4_8", &gAlgImbProjA1C4_8_, "gAlgImbProjA1C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA1C8_100", &gAlgImbProjA1C8_100_, "gAlgImbProjA1C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA2CF", &gAlgImbProjA2CF_, "gAlgImbProjA2CF[5]/F");
    genTree_p->Branch("gAlgImbProjA2C0_1", &gAlgImbProjA2C0_1_, "gAlgImbProjA2C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA2C1_2", &gAlgImbProjA2C1_2_, "gAlgImbProjA2C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA2C2_4", &gAlgImbProjA2C2_4_, "gAlgImbProjA2C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA2C4_8", &gAlgImbProjA2C4_8_, "gAlgImbProjA2C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA2C8_100", &gAlgImbProjA2C8_100_, "gAlgImbProjA2C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA3CF", &gAlgImbProjA3CF_, "gAlgImbProjA3CF[5]/F");
    genTree_p->Branch("gAlgImbProjA3C0_1", &gAlgImbProjA3C0_1_, "gAlgImbProjA3C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA3C1_2", &gAlgImbProjA3C1_2_, "gAlgImbProjA3C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA3C2_4", &gAlgImbProjA3C2_4_, "gAlgImbProjA3C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA3C4_8", &gAlgImbProjA3C4_8_, "gAlgImbProjA3C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA3C8_100", &gAlgImbProjA3C8_100_, "gAlgImbProjA3C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA4CF", &gAlgImbProjA4CF_, "gAlgImbProjA4CF[5]/F");
    genTree_p->Branch("gAlgImbProjA4C0_1", &gAlgImbProjA4C0_1_, "gAlgImbProjA4C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA4C1_2", &gAlgImbProjA4C1_2_, "gAlgImbProjA4C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA4C2_4", &gAlgImbProjA4C2_4_, "gAlgImbProjA4C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA4C4_8", &gAlgImbProjA4C4_8_, "gAlgImbProjA4C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA4C8_100", &gAlgImbProjA4C8_100_, "gAlgImbProjA4C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA5CF", &gAlgImbProjA5CF_, "gAlgImbProjA5CF[5]/F");
    genTree_p->Branch("gAlgImbProjA5C0_1", &gAlgImbProjA5C0_1_, "gAlgImbProjA5C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA5C1_2", &gAlgImbProjA5C1_2_, "gAlgImbProjA5C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA5C2_4", &gAlgImbProjA5C2_4_, "gAlgImbProjA5C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA5C4_8", &gAlgImbProjA5C4_8_, "gAlgImbProjA5C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA5C8_100", &gAlgImbProjA5C8_100_, "gAlgImbProjA5C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA6CF", &gAlgImbProjA6CF_, "gAlgImbProjA6CF[5]/F");
    genTree_p->Branch("gAlgImbProjA6C0_1", &gAlgImbProjA6C0_1_, "gAlgImbProjA6C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA6C1_2", &gAlgImbProjA6C1_2_, "gAlgImbProjA6C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA6C2_4", &gAlgImbProjA6C2_4_, "gAlgImbProjA6C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA6C4_8", &gAlgImbProjA6C4_8_, "gAlgImbProjA6C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA6C8_100", &gAlgImbProjA6C8_100_, "gAlgImbProjA6C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA7CF", &gAlgImbProjA7CF_, "gAlgImbProjA7CF[5]/F");
    genTree_p->Branch("gAlgImbProjA7C0_1", &gAlgImbProjA7C0_1_, "gAlgImbProjA7C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA7C1_2", &gAlgImbProjA7C1_2_, "gAlgImbProjA7C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA7C2_4", &gAlgImbProjA7C2_4_, "gAlgImbProjA7C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA7C4_8", &gAlgImbProjA7C4_8_, "gAlgImbProjA7C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA7C8_100", &gAlgImbProjA7C8_100_, "gAlgImbProjA7C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA8CF", &gAlgImbProjA8CF_, "gAlgImbProjA8CF[5]/F");
    genTree_p->Branch("gAlgImbProjA8C0_1", &gAlgImbProjA8C0_1_, "gAlgImbProjA8C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA8C1_2", &gAlgImbProjA8C1_2_, "gAlgImbProjA8C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA8C2_4", &gAlgImbProjA8C2_4_, "gAlgImbProjA8C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA8C4_8", &gAlgImbProjA8C4_8_, "gAlgImbProjA8C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA8C8_100", &gAlgImbProjA8C8_100_, "gAlgImbProjA8C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA9CF", &gAlgImbProjA9CF_, "gAlgImbProjA9CF[5]/F");
    genTree_p->Branch("gAlgImbProjA9C0_1", &gAlgImbProjA9C0_1_, "gAlgImbProjA9C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA9C1_2", &gAlgImbProjA9C1_2_, "gAlgImbProjA9C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA9C2_4", &gAlgImbProjA9C2_4_, "gAlgImbProjA9C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA9C4_8", &gAlgImbProjA9C4_8_, "gAlgImbProjA9C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA9C8_100", &gAlgImbProjA9C8_100_, "gAlgImbProjA9C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA10CF", &gAlgImbProjA10CF_, "gAlgImbProjA10CF[5]/F");
    genTree_p->Branch("gAlgImbProjA10C0_1", &gAlgImbProjA10C0_1_, "gAlgImbProjA10C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA10C1_2", &gAlgImbProjA10C1_2_, "gAlgImbProjA10C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA10C2_4", &gAlgImbProjA10C2_4_, "gAlgImbProjA10C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA10C4_8", &gAlgImbProjA10C4_8_, "gAlgImbProjA10C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA10C8_100", &gAlgImbProjA10C8_100_, "gAlgImbProjA10C8_100[5]/F");


    genTree_p->Branch("gAlgImbProjA21CF", &gAlgImbProjA21CF_, "gAlgImbProjA21CF[5]/F");
    genTree_p->Branch("gAlgImbProjA21C0_1", &gAlgImbProjA21C0_1_, "gAlgImbProjA21C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA21C1_2", &gAlgImbProjA21C1_2_, "gAlgImbProjA21C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA21C2_4", &gAlgImbProjA21C2_4_, "gAlgImbProjA21C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA21C4_8", &gAlgImbProjA21C4_8_, "gAlgImbProjA21C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA21C8_100", &gAlgImbProjA21C8_100_, "gAlgImbProjA21C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA22CF", &gAlgImbProjA22CF_, "gAlgImbProjA22CF[5]/F");
    genTree_p->Branch("gAlgImbProjA22C0_1", &gAlgImbProjA22C0_1_, "gAlgImbProjA22C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA22C1_2", &gAlgImbProjA22C1_2_, "gAlgImbProjA22C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA22C2_4", &gAlgImbProjA22C2_4_, "gAlgImbProjA22C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA22C4_8", &gAlgImbProjA22C4_8_, "gAlgImbProjA22C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA22C8_100", &gAlgImbProjA22C8_100_, "gAlgImbProjA22C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA23CF", &gAlgImbProjA23CF_, "gAlgImbProjA23CF[5]/F");
    genTree_p->Branch("gAlgImbProjA23C0_1", &gAlgImbProjA23C0_1_, "gAlgImbProjA23C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA23C1_2", &gAlgImbProjA23C1_2_, "gAlgImbProjA23C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA23C2_4", &gAlgImbProjA23C2_4_, "gAlgImbProjA23C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA23C4_8", &gAlgImbProjA23C4_8_, "gAlgImbProjA23C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA23C8_100", &gAlgImbProjA23C8_100_, "gAlgImbProjA23C8_100[5]/F");

    genTree_p->Branch("gAlgImbProjA24CF", &gAlgImbProjA24CF_, "gAlgImbProjA24CF[5]/F");
    genTree_p->Branch("gAlgImbProjA24C0_1", &gAlgImbProjA24C0_1_, "gAlgImbProjA24C0_1[5]/F");
    genTree_p->Branch("gAlgImbProjA24C1_2", &gAlgImbProjA24C1_2_, "gAlgImbProjA24C1_2[5]/F");
    genTree_p->Branch("gAlgImbProjA24C2_4", &gAlgImbProjA24C2_4_, "gAlgImbProjA24C2_4[5]/F");
    genTree_p->Branch("gAlgImbProjA24C4_8", &gAlgImbProjA24C4_8_, "gAlgImbProjA24C4_8[5]/F");
    genTree_p->Branch("gAlgImbProjA24C8_100", &gAlgImbProjA24C8_100_, "gAlgImbProjA24C8_100[5]/F");

  }
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
  for(Int_t initIter = 0; initIter < 5; initIter++){
    eventSet_[initIter] = false;
    centWeight_[initIter] = -1;
    centWeight_2pi3_[initIter] = -1;
    centWeight_120_5pi6_[initIter] = -1;

    AlgLeadJtPt_[initIter] = -10;
    AlgSubLeadJtPt_[initIter] = -10;
    AlgThirdJtPt_[initIter] = -10;
    AlgFourthJtPt_[initIter] = -10;
    AlgLeadJtPhi_[initIter] = -10;
    AlgSubLeadJtPhi_[initIter] = -10;
    AlgThirdJtPhi_[initIter] = -10;
    AlgFourthJtPhi_[initIter] = -10;
    AlgLeadJtEta_[initIter] = -10;
    AlgSubLeadJtEta_[initIter] = -10;
    AlgThirdJtEta_[initIter] = -10;
    AlgFourthJtEta_[initIter] = -10;

    AlgJtAvePhi_[initIter] = -10;
    AlgJtDelPhi_[initIter] = -10;
    AlgJtAsymm_[initIter] = -10;

    if(montecarlo){
      isQuarkJet_[initIter] = false;
      isGluonJet_[initIter] = false;

      pthatWeight_ = -10;
      pthat_ = -10;

      AlgLeadRefPt_[initIter] = -10;
      AlgSubLeadRefPt_[initIter] = -10;
      AlgThirdRefPt_[initIter] = -10;
      AlgFourthRefPt_[initIter] = -10;

      AlgLeadRefPhi_[initIter] = -10;
      AlgSubLeadRefPhi_[initIter] = -10;
      AlgThirdRefPhi_[initIter] = -10;
      AlgFourthRefPhi_[initIter] = -10;

      AlgLeadRefEta_[initIter] = -10;
      AlgSubLeadRefEta_[initIter] = -10;
      AlgThirdRefEta_[initIter] = -10;
      AlgFourthRefEta_[initIter] = -10;

      AlgRefDelPhi_[initIter] = -10;
      AlgRefAsymm_[initIter] = -10;
    }

  }
}

void InitProjPerp(bool montecarlo = false)
{
  //Tracks proj. onto Alg, ordered according to enum above, corr in the back 5, All, Cone, and NotCone

  for(Int_t initIter = 0; initIter < 10; initIter++){
    rAlgImbProjAF_[initIter] = 0;
    rAlgImbProjA0_1_[initIter] = 0;
    rAlgImbProjA1_2_[initIter] = 0;
    rAlgImbProjA2_4_[initIter] = 0;
    rAlgImbProjA4_8_[initIter] = 0;
    rAlgImbProjA8_100_[initIter] = 0;
    rAlgImbProjACF_[initIter] = 0;
    rAlgImbProjAC0_1_[initIter] = 0;
    rAlgImbProjAC1_2_[initIter] = 0;
    rAlgImbProjAC2_4_[initIter] = 0;
    rAlgImbProjAC4_8_[initIter] = 0;
    rAlgImbProjAC8_100_[initIter] = 0;
    rAlgImbProjANCF_[initIter] = 0;
    rAlgImbProjANC0_1_[initIter] = 0;
    rAlgImbProjANC1_2_[initIter] = 0;
    rAlgImbProjANC2_4_[initIter] = 0;
    rAlgImbProjANC4_8_[initIter] = 0;
    rAlgImbProjANC8_100_[initIter] = 0;
    rAlgImbProjANCCutF_[initIter] = 0;
    rAlgImbProjANCCut0_1_[initIter] = 0;
    rAlgImbProjANCCut1_2_[initIter] = 0;
    rAlgImbProjANCCut2_4_[initIter] = 0;
    rAlgImbProjANCCut4_8_[initIter] = 0;
    rAlgImbProjANCCut8_100_[initIter] = 0;

    //DelRs Proj

    //ProjA

    rAlgImbProjA1CF_[initIter] = 0;
    rAlgImbProjA1C0_1_[initIter] = 0;
    rAlgImbProjA1C1_2_[initIter] = 0;
    rAlgImbProjA1C2_4_[initIter] = 0;
    rAlgImbProjA1C4_8_[initIter] = 0;
    rAlgImbProjA1C8_100_[initIter] = 0;

    rAlgImbProjA2CF_[initIter] = 0;
    rAlgImbProjA2C0_1_[initIter] = 0;
    rAlgImbProjA2C1_2_[initIter] = 0;
    rAlgImbProjA2C2_4_[initIter] = 0;
    rAlgImbProjA2C4_8_[initIter] = 0;
    rAlgImbProjA2C8_100_[initIter] = 0;

    rAlgImbProjA3CF_[initIter] = 0;
    rAlgImbProjA3C0_1_[initIter] = 0;
    rAlgImbProjA3C1_2_[initIter] = 0;
    rAlgImbProjA3C2_4_[initIter] = 0;
    rAlgImbProjA3C4_8_[initIter] = 0;
    rAlgImbProjA3C8_100_[initIter] = 0;

    rAlgImbProjA4CF_[initIter] = 0;
    rAlgImbProjA4C0_1_[initIter] = 0;
    rAlgImbProjA4C1_2_[initIter] = 0;
    rAlgImbProjA4C2_4_[initIter] = 0;
    rAlgImbProjA4C4_8_[initIter] = 0;
    rAlgImbProjA4C8_100_[initIter] = 0;

    rAlgImbProjA5CF_[initIter] = 0;
    rAlgImbProjA5C0_1_[initIter] = 0;
    rAlgImbProjA5C1_2_[initIter] = 0;
    rAlgImbProjA5C2_4_[initIter] = 0;
    rAlgImbProjA5C4_8_[initIter] = 0;
    rAlgImbProjA5C8_100_[initIter] = 0;

    rAlgImbProjA6CF_[initIter] = 0;
    rAlgImbProjA6C0_1_[initIter] = 0;
    rAlgImbProjA6C1_2_[initIter] = 0;
    rAlgImbProjA6C2_4_[initIter] = 0;
    rAlgImbProjA6C4_8_[initIter] = 0;
    rAlgImbProjA6C8_100_[initIter] = 0;

    rAlgImbProjA7CF_[initIter] = 0;
    rAlgImbProjA7C0_1_[initIter] = 0;
    rAlgImbProjA7C1_2_[initIter] = 0;
    rAlgImbProjA7C2_4_[initIter] = 0;
    rAlgImbProjA7C4_8_[initIter] = 0;
    rAlgImbProjA7C8_100_[initIter] = 0;

    rAlgImbProjA8CF_[initIter] = 0;
    rAlgImbProjA8C0_1_[initIter] = 0;
    rAlgImbProjA8C1_2_[initIter] = 0;
    rAlgImbProjA8C2_4_[initIter] = 0;
    rAlgImbProjA8C4_8_[initIter] = 0;
    rAlgImbProjA8C8_100_[initIter] = 0;

    rAlgImbProjA9CF_[initIter] = 0;
    rAlgImbProjA9C0_1_[initIter] = 0;
    rAlgImbProjA9C1_2_[initIter] = 0;
    rAlgImbProjA9C2_4_[initIter] = 0;
    rAlgImbProjA9C4_8_[initIter] = 0;
    rAlgImbProjA9C8_100_[initIter] = 0;

    rAlgImbProjA10CF_[initIter] = 0;
    rAlgImbProjA10C0_1_[initIter] = 0;
    rAlgImbProjA10C1_2_[initIter] = 0;
    rAlgImbProjA10C2_4_[initIter] = 0;
    rAlgImbProjA10C4_8_[initIter] = 0;
    rAlgImbProjA10C8_100_[initIter] = 0;

    //delR plots 2

    rAlgImbProjA21CF_[initIter] = 0;
    rAlgImbProjA21C0_1_[initIter] = 0;
    rAlgImbProjA21C1_2_[initIter] = 0;
    rAlgImbProjA21C2_4_[initIter] = 0;
    rAlgImbProjA21C4_8_[initIter] = 0;
    rAlgImbProjA21C8_100_[initIter] = 0;

    rAlgImbProjA22CF_[initIter] = 0;
    rAlgImbProjA22C0_1_[initIter] = 0;
    rAlgImbProjA22C1_2_[initIter] = 0;
    rAlgImbProjA22C2_4_[initIter] = 0;
    rAlgImbProjA22C4_8_[initIter] = 0;
    rAlgImbProjA22C8_100_[initIter] = 0;

    rAlgImbProjA23CF_[initIter] = 0;
    rAlgImbProjA23C0_1_[initIter] = 0;
    rAlgImbProjA23C1_2_[initIter] = 0;
    rAlgImbProjA23C2_4_[initIter] = 0;
    rAlgImbProjA23C4_8_[initIter] = 0;
    rAlgImbProjA23C8_100_[initIter] = 0;

    rAlgImbProjA24CF_[initIter] = 0;
    rAlgImbProjA24C0_1_[initIter] = 0;
    rAlgImbProjA24C1_2_[initIter] = 0;
    rAlgImbProjA24C2_4_[initIter] = 0;
    rAlgImbProjA24C4_8_[initIter] = 0;
    rAlgImbProjA24C8_100_[initIter] = 0;
  }

  if(montecarlo){
    //Gen. proj. onto Truth
    for(Int_t initIter = 0; initIter < 5; initIter++){

      gAlgImbProjAF_[initIter] = 0;
      gAlgImbProjA0_1_[initIter] = 0;
      gAlgImbProjA1_2_[initIter] = 0;
      gAlgImbProjA2_4_[initIter] = 0;
      gAlgImbProjA4_8_[initIter] = 0;
      gAlgImbProjA8_100_[initIter] = 0;
      gAlgImbProjACF_[initIter] = 0;
      gAlgImbProjAC0_1_[initIter] = 0;
      gAlgImbProjAC1_2_[initIter] = 0;
      gAlgImbProjAC2_4_[initIter] = 0;
      gAlgImbProjAC4_8_[initIter] = 0;
      gAlgImbProjAC8_100_[initIter] = 0;
      gAlgImbProjANCF_[initIter] = 0;
      gAlgImbProjANC0_1_[initIter] = 0;
      gAlgImbProjANC1_2_[initIter] = 0;
      gAlgImbProjANC2_4_[initIter] = 0;
      gAlgImbProjANC4_8_[initIter] = 0;
      gAlgImbProjANC8_100_[initIter] = 0;
      gAlgImbProjANCCutF_[initIter] = 0;
      gAlgImbProjANCCut0_1_[initIter] = 0;
      gAlgImbProjANCCut1_2_[initIter] = 0;
      gAlgImbProjANCCut2_4_[initIter] = 0;
      gAlgImbProjANCCut4_8_[initIter] = 0;
      gAlgImbProjANCCut8_100_[initIter] = 0;

      gAlgImbProjSigACF_[initIter] = 0;
      gAlgImbProjSigAC0_1_[initIter] = 0;
      gAlgImbProjSigAC1_2_[initIter] = 0;
      gAlgImbProjSigAC2_4_[initIter] = 0;
      gAlgImbProjSigAC4_8_[initIter] = 0;
      gAlgImbProjSigAC8_100_[initIter] = 0;
      gAlgImbProjSigANCF_[initIter] = 0;
      gAlgImbProjSigANC0_1_[initIter] = 0;
      gAlgImbProjSigANC1_2_[initIter] = 0;
      gAlgImbProjSigANC2_4_[initIter] = 0;
      gAlgImbProjSigANC4_8_[initIter] = 0;
      gAlgImbProjSigANC8_100_[initIter] = 0;

      gAlgImbProjUEACF_[initIter] = 0;
      gAlgImbProjUEAC0_1_[initIter] = 0;
      gAlgImbProjUEAC1_2_[initIter] = 0;
      gAlgImbProjUEAC2_4_[initIter] = 0;
      gAlgImbProjUEAC4_8_[initIter] = 0;
      gAlgImbProjUEAC8_100_[initIter] = 0;
      gAlgImbProjUEANCF_[initIter] = 0;
      gAlgImbProjUEANC0_1_[initIter] = 0;
      gAlgImbProjUEANC1_2_[initIter] = 0;
      gAlgImbProjUEANC2_4_[initIter] = 0;
      gAlgImbProjUEANC4_8_[initIter] = 0;
      gAlgImbProjUEANC8_100_[initIter] = 0;

      //DelRs
      //proj As

      gAlgImbProjA1CF_[initIter] = 0;
      gAlgImbProjA1C0_1_[initIter] = 0;
      gAlgImbProjA1C1_2_[initIter] = 0;
      gAlgImbProjA1C2_4_[initIter] = 0;
      gAlgImbProjA1C4_8_[initIter] = 0;
      gAlgImbProjA1C8_100_[initIter] = 0;

      gAlgImbProjA2CF_[initIter] = 0;
      gAlgImbProjA2C0_1_[initIter] = 0;
      gAlgImbProjA2C1_2_[initIter] = 0;
      gAlgImbProjA2C2_4_[initIter] = 0;
      gAlgImbProjA2C4_8_[initIter] = 0;
      gAlgImbProjA2C8_100_[initIter] = 0;

      gAlgImbProjA3CF_[initIter] = 0;
      gAlgImbProjA3C0_1_[initIter] = 0;
      gAlgImbProjA3C1_2_[initIter] = 0;
      gAlgImbProjA3C2_4_[initIter] = 0;
      gAlgImbProjA3C4_8_[initIter] = 0;
      gAlgImbProjA3C8_100_[initIter] = 0;

      gAlgImbProjA4CF_[initIter] = 0;
      gAlgImbProjA4C0_1_[initIter] = 0;
      gAlgImbProjA4C1_2_[initIter] = 0;
      gAlgImbProjA4C2_4_[initIter] = 0;
      gAlgImbProjA4C4_8_[initIter] = 0;
      gAlgImbProjA4C8_100_[initIter] = 0;

      gAlgImbProjA5CF_[initIter] = 0;
      gAlgImbProjA5C0_1_[initIter] = 0;
      gAlgImbProjA5C1_2_[initIter] = 0;
      gAlgImbProjA5C2_4_[initIter] = 0;
      gAlgImbProjA5C4_8_[initIter] = 0;
      gAlgImbProjA5C8_100_[initIter] = 0;

      gAlgImbProjA6CF_[initIter] = 0;
      gAlgImbProjA6C0_1_[initIter] = 0;
      gAlgImbProjA6C1_2_[initIter] = 0;
      gAlgImbProjA6C2_4_[initIter] = 0;
      gAlgImbProjA6C4_8_[initIter] = 0;
      gAlgImbProjA6C8_100_[initIter] = 0;

      gAlgImbProjA7CF_[initIter] = 0;
      gAlgImbProjA7C0_1_[initIter] = 0;
      gAlgImbProjA7C1_2_[initIter] = 0;
      gAlgImbProjA7C2_4_[initIter] = 0;
      gAlgImbProjA7C4_8_[initIter] = 0;
      gAlgImbProjA7C8_100_[initIter] = 0;

      gAlgImbProjA8CF_[initIter] = 0;
      gAlgImbProjA8C0_1_[initIter] = 0;
      gAlgImbProjA8C1_2_[initIter] = 0;
      gAlgImbProjA8C2_4_[initIter] = 0;
      gAlgImbProjA8C4_8_[initIter] = 0;
      gAlgImbProjA8C8_100_[initIter] = 0;

      gAlgImbProjA9CF_[initIter] = 0;
      gAlgImbProjA9C0_1_[initIter] = 0;
      gAlgImbProjA9C1_2_[initIter] = 0;
      gAlgImbProjA9C2_4_[initIter] = 0;
      gAlgImbProjA9C4_8_[initIter] = 0;
      gAlgImbProjA9C8_100_[initIter] = 0;

      gAlgImbProjA10CF_[initIter] = 0;
      gAlgImbProjA10C0_1_[initIter] = 0;
      gAlgImbProjA10C1_2_[initIter] = 0;
      gAlgImbProjA10C2_4_[initIter] = 0;
      gAlgImbProjA10C4_8_[initIter] = 0;
      gAlgImbProjA10C8_100_[initIter] = 0;


      gAlgImbProjA21CF_[initIter] = 0;
      gAlgImbProjA21C0_1_[initIter] = 0;
      gAlgImbProjA21C1_2_[initIter] = 0;
      gAlgImbProjA21C2_4_[initIter] = 0;
      gAlgImbProjA21C4_8_[initIter] = 0;
      gAlgImbProjA21C8_100_[initIter] = 0;

      gAlgImbProjA22CF_[initIter] = 0;
      gAlgImbProjA22C0_1_[initIter] = 0;
      gAlgImbProjA22C1_2_[initIter] = 0;
      gAlgImbProjA22C2_4_[initIter] = 0;
      gAlgImbProjA22C4_8_[initIter] = 0;
      gAlgImbProjA22C8_100_[initIter] = 0;

      gAlgImbProjA23CF_[initIter] = 0;
      gAlgImbProjA23C0_1_[initIter] = 0;
      gAlgImbProjA23C1_2_[initIter] = 0;
      gAlgImbProjA23C2_4_[initIter] = 0;
      gAlgImbProjA23C4_8_[initIter] = 0;
      gAlgImbProjA23C8_100_[initIter] = 0;

      gAlgImbProjA24CF_[initIter] = 0;
      gAlgImbProjA24C0_1_[initIter] = 0;
      gAlgImbProjA24C1_2_[initIter] = 0;
      gAlgImbProjA24C2_4_[initIter] = 0;
      gAlgImbProjA24C4_8_[initIter] = 0;
      gAlgImbProjA24C8_100_[initIter] = 0;
    }    
  }

}

#endif
