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


enum AlgoType {
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

Float_t rAlgImbProjF_[10];
Float_t rAlgImbPerpF_[10];
Float_t rAlgImbProj0_1_[10];
Float_t rAlgImbProj1_2_[10];
Float_t rAlgImbProj2_4_[10];
Float_t rAlgImbProj4_8_[10];
Float_t rAlgImbProj8_100_[10];
Float_t rAlgImbProjCF_[10];
Float_t rAlgImbPerpCF_[10];
Float_t rAlgImbProjC0_1_[10];
Float_t rAlgImbProjC1_2_[10];
Float_t rAlgImbProjC2_4_[10];
Float_t rAlgImbProjC4_8_[10];
Float_t rAlgImbProjC8_100_[10];
Float_t rAlgImbProjNCF_[10];
Float_t rAlgImbPerpNCF_[10];
Float_t rAlgImbProjNC0_1_[10];
Float_t rAlgImbProjNC1_2_[10];
Float_t rAlgImbProjNC2_4_[10];
Float_t rAlgImbProjNC4_8_[10];
Float_t rAlgImbProjNC8_100_[10];

Float_t rAlgImbHemF_[10];
Float_t rAlgImbHem0_1_[10];
Float_t rAlgImbHem1_2_[10];
Float_t rAlgImbHem2_4_[10];
Float_t rAlgImbHem4_8_[10];
Float_t rAlgImbHem8_100_[10];
Float_t rAlgImbHemCF_[10];
Float_t rAlgImbHemC0_1_[10];
Float_t rAlgImbHemC1_2_[10];
Float_t rAlgImbHemC2_4_[10];
Float_t rAlgImbHemC4_8_[10];
Float_t rAlgImbHemC8_100_[10];
Float_t rAlgImbHemNCF_[10];
Float_t rAlgImbHemNC0_1_[10];
Float_t rAlgImbHemNC1_2_[10];
Float_t rAlgImbHemNC2_4_[10];
Float_t rAlgImbHemNC4_8_[10];
Float_t rAlgImbHemNC8_100_[10];


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

//DelR's

Float_t rAlgImbProj1CF_[10];
Float_t rAlgImbProj1C0_1_[10];
Float_t rAlgImbProj1C1_2_[10];
Float_t rAlgImbProj1C2_4_[10];
Float_t rAlgImbProj1C4_8_[10];
Float_t rAlgImbProj1C8_100_[10];

Float_t rAlgImbProj2CF_[10];
Float_t rAlgImbProj2C0_1_[10];
Float_t rAlgImbProj2C1_2_[10];
Float_t rAlgImbProj2C2_4_[10];
Float_t rAlgImbProj2C4_8_[10];
Float_t rAlgImbProj2C8_100_[10];

Float_t rAlgImbProj3CF_[10];
Float_t rAlgImbProj3C0_1_[10];
Float_t rAlgImbProj3C1_2_[10];
Float_t rAlgImbProj3C2_4_[10];
Float_t rAlgImbProj3C4_8_[10];
Float_t rAlgImbProj3C8_100_[10];

Float_t rAlgImbProj4CF_[10];
Float_t rAlgImbProj4C0_1_[10];
Float_t rAlgImbProj4C1_2_[10];
Float_t rAlgImbProj4C2_4_[10];
Float_t rAlgImbProj4C4_8_[10];
Float_t rAlgImbProj4C8_100_[10];

Float_t rAlgImbProj5CF_[10];
Float_t rAlgImbProj5C0_1_[10];
Float_t rAlgImbProj5C1_2_[10];
Float_t rAlgImbProj5C2_4_[10];
Float_t rAlgImbProj5C4_8_[10];
Float_t rAlgImbProj5C8_100_[10];

Float_t rAlgImbProj6CF_[10];
Float_t rAlgImbProj6C0_1_[10];
Float_t rAlgImbProj6C1_2_[10];
Float_t rAlgImbProj6C2_4_[10];
Float_t rAlgImbProj6C4_8_[10];
Float_t rAlgImbProj6C8_100_[10];

Float_t rAlgImbProj7CF_[10];
Float_t rAlgImbProj7C0_1_[10];
Float_t rAlgImbProj7C1_2_[10];
Float_t rAlgImbProj7C2_4_[10];
Float_t rAlgImbProj7C4_8_[10];
Float_t rAlgImbProj7C8_100_[10];

Float_t rAlgImbProj8CF_[10];
Float_t rAlgImbProj8C0_1_[10];
Float_t rAlgImbProj8C1_2_[10];
Float_t rAlgImbProj8C2_4_[10];
Float_t rAlgImbProj8C4_8_[10];
Float_t rAlgImbProj8C8_100_[10];

Float_t rAlgImbProj9CF_[10];
Float_t rAlgImbProj9C0_1_[10];
Float_t rAlgImbProj9C1_2_[10];
Float_t rAlgImbProj9C2_4_[10];
Float_t rAlgImbProj9C4_8_[10];
Float_t rAlgImbProj9C8_100_[10];

Float_t rAlgImbProj10CF_[10];
Float_t rAlgImbProj10C0_1_[10];
Float_t rAlgImbProj10C1_2_[10];
Float_t rAlgImbProj10C2_4_[10];
Float_t rAlgImbProj10C4_8_[10];
Float_t rAlgImbProj10C8_100_[10];

//Del R hem

Float_t rAlgImbHem1CF_[10];
Float_t rAlgImbHem1C0_1_[10];
Float_t rAlgImbHem1C1_2_[10];
Float_t rAlgImbHem1C2_4_[10];
Float_t rAlgImbHem1C4_8_[10];
Float_t rAlgImbHem1C8_100_[10];

Float_t rAlgImbHem2CF_[10];
Float_t rAlgImbHem2C0_1_[10];
Float_t rAlgImbHem2C1_2_[10];
Float_t rAlgImbHem2C2_4_[10];
Float_t rAlgImbHem2C4_8_[10];
Float_t rAlgImbHem2C8_100_[10];

Float_t rAlgImbHem3CF_[10];
Float_t rAlgImbHem3C0_1_[10];
Float_t rAlgImbHem3C1_2_[10];
Float_t rAlgImbHem3C2_4_[10];
Float_t rAlgImbHem3C4_8_[10];
Float_t rAlgImbHem3C8_100_[10];

Float_t rAlgImbHem4CF_[10];
Float_t rAlgImbHem4C0_1_[10];
Float_t rAlgImbHem4C1_2_[10];
Float_t rAlgImbHem4C2_4_[10];
Float_t rAlgImbHem4C4_8_[10];
Float_t rAlgImbHem4C8_100_[10];

Float_t rAlgImbHem5CF_[10];
Float_t rAlgImbHem5C0_1_[10];
Float_t rAlgImbHem5C1_2_[10];
Float_t rAlgImbHem5C2_4_[10];
Float_t rAlgImbHem5C4_8_[10];
Float_t rAlgImbHem5C8_100_[10];

Float_t rAlgImbHem6CF_[10];
Float_t rAlgImbHem6C0_1_[10];
Float_t rAlgImbHem6C1_2_[10];
Float_t rAlgImbHem6C2_4_[10];
Float_t rAlgImbHem6C4_8_[10];
Float_t rAlgImbHem6C8_100_[10];

Float_t rAlgImbHem7CF_[10];
Float_t rAlgImbHem7C0_1_[10];
Float_t rAlgImbHem7C1_2_[10];
Float_t rAlgImbHem7C2_4_[10];
Float_t rAlgImbHem7C4_8_[10];
Float_t rAlgImbHem7C8_100_[10];

Float_t rAlgImbHem8CF_[10];
Float_t rAlgImbHem8C0_1_[10];
Float_t rAlgImbHem8C1_2_[10];
Float_t rAlgImbHem8C2_4_[10];
Float_t rAlgImbHem8C4_8_[10];
Float_t rAlgImbHem8C8_100_[10];

Float_t rAlgImbHem9CF_[10];
Float_t rAlgImbHem9C0_1_[10];
Float_t rAlgImbHem9C1_2_[10];
Float_t rAlgImbHem9C2_4_[10];
Float_t rAlgImbHem9C4_8_[10];
Float_t rAlgImbHem9C8_100_[10];

Float_t rAlgImbHem10CF_[10];
Float_t rAlgImbHem10C0_1_[10];
Float_t rAlgImbHem10C1_2_[10];
Float_t rAlgImbHem10C2_4_[10];
Float_t rAlgImbHem10C4_8_[10];
Float_t rAlgImbHem10C8_100_[10];

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

//Event Set Bool array, [0] == PuPF, [1] == PuCalo, .etc according to enum

Bool_t eventSet_[5];
Float_t setWeight_[5];
Float_t setWeight_2pi3_[5];

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

Float_t gAlgImbProjF_[5];
Float_t gAlgImbPerpF_[5];
Float_t gAlgImbProj0_1_[5];
Float_t gAlgImbProj1_2_[5];
Float_t gAlgImbProj2_4_[5];
Float_t gAlgImbProj4_8_[5];
Float_t gAlgImbProj8_100_[5];
Float_t gAlgImbProjCF_[5];
Float_t gAlgImbPerpCF_[5];
Float_t gAlgImbProjC0_1_[5];
Float_t gAlgImbProjC1_2_[5];
Float_t gAlgImbProjC2_4_[5];
Float_t gAlgImbProjC4_8_[5];
Float_t gAlgImbProjC8_100_[5];
Float_t gAlgImbProjNCF_[5];
Float_t gAlgImbPerpNCF_[5];
Float_t gAlgImbProjNC0_1_[5];
Float_t gAlgImbProjNC1_2_[5];
Float_t gAlgImbProjNC2_4_[5];
Float_t gAlgImbProjNC4_8_[5];
Float_t gAlgImbProjNC8_100_[5];


Float_t gAlgImbProjSigCF_[5];
Float_t gAlgImbPerpSigCF_[5];
Float_t gAlgImbProjSigC0_1_[5];
Float_t gAlgImbProjSigC1_2_[5];
Float_t gAlgImbProjSigC2_4_[5];
Float_t gAlgImbProjSigC4_8_[5];
Float_t gAlgImbProjSigC8_100_[5];
Float_t gAlgImbProjSigNCF_[5];
Float_t gAlgImbPerpSigNCF_[5];
Float_t gAlgImbProjSigNC0_1_[5];
Float_t gAlgImbProjSigNC1_2_[5];
Float_t gAlgImbProjSigNC2_4_[5];
Float_t gAlgImbProjSigNC4_8_[5];
Float_t gAlgImbProjSigNC8_100_[5];

Float_t gAlgImbProjUECF_[5];
Float_t gAlgImbPerpUECF_[5];
Float_t gAlgImbProjUEC0_1_[5];
Float_t gAlgImbProjUEC1_2_[5];
Float_t gAlgImbProjUEC2_4_[5];
Float_t gAlgImbProjUEC4_8_[5];
Float_t gAlgImbProjUEC8_100_[5];
Float_t gAlgImbProjUENCF_[5];
Float_t gAlgImbPerpUENCF_[5];
Float_t gAlgImbProjUENC0_1_[5];
Float_t gAlgImbProjUENC1_2_[5];
Float_t gAlgImbProjUENC2_4_[5];
Float_t gAlgImbProjUENC4_8_[5];
Float_t gAlgImbProjUENC8_100_[5];


Float_t gAlgImbHemF_[5];
Float_t gAlgImbHem0_1_[5];
Float_t gAlgImbHem1_2_[5];
Float_t gAlgImbHem2_4_[5];
Float_t gAlgImbHem4_8_[5];
Float_t gAlgImbHem8_100_[5];
Float_t gAlgImbHemCF_[5];
Float_t gAlgImbHemC0_1_[5];
Float_t gAlgImbHemC1_2_[5];
Float_t gAlgImbHemC2_4_[5];
Float_t gAlgImbHemC4_8_[5];
Float_t gAlgImbHemC8_100_[5];
Float_t gAlgImbHemNCF_[5];
Float_t gAlgImbHemNC0_1_[5];
Float_t gAlgImbHemNC1_2_[5];
Float_t gAlgImbHemNC2_4_[5];
Float_t gAlgImbHemNC4_8_[5];
Float_t gAlgImbHemNC8_100_[5];

Float_t gAlgImbHemSigCF_[5];
Float_t gAlgImbHemSigC0_1_[5];
Float_t gAlgImbHemSigC1_2_[5];
Float_t gAlgImbHemSigC2_4_[5];
Float_t gAlgImbHemSigC4_8_[5];
Float_t gAlgImbHemSigC8_100_[5];
Float_t gAlgImbHemSigNCF_[5];
Float_t gAlgImbHemSigNC0_1_[5];
Float_t gAlgImbHemSigNC1_2_[5];
Float_t gAlgImbHemSigNC2_4_[5];
Float_t gAlgImbHemSigNC4_8_[5];
Float_t gAlgImbHemSigNC8_100_[5];

Float_t gAlgImbHemUECF_[5];
Float_t gAlgImbHemUEC0_1_[5];
Float_t gAlgImbHemUEC1_2_[5];
Float_t gAlgImbHemUEC2_4_[5];
Float_t gAlgImbHemUEC4_8_[5];
Float_t gAlgImbHemUEC8_100_[5];
Float_t gAlgImbHemUENCF_[5];
Float_t gAlgImbHemUENC0_1_[5];
Float_t gAlgImbHemUENC1_2_[5];
Float_t gAlgImbHemUENC2_4_[5];
Float_t gAlgImbHemUENC4_8_[5];
Float_t gAlgImbHemUENC8_100_[5];

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

Float_t gAlgImbProj1CF_[5];
Float_t gAlgImbProj1C0_1_[5];
Float_t gAlgImbProj1C1_2_[5];
Float_t gAlgImbProj1C2_4_[5];
Float_t gAlgImbProj1C4_8_[5];
Float_t gAlgImbProj1C8_100_[5];

Float_t gAlgImbProj2CF_[5];
Float_t gAlgImbProj2C0_1_[5];
Float_t gAlgImbProj2C1_2_[5];
Float_t gAlgImbProj2C2_4_[5];
Float_t gAlgImbProj2C4_8_[5];
Float_t gAlgImbProj2C8_100_[5];

Float_t gAlgImbProj3CF_[5];
Float_t gAlgImbProj3C0_1_[5];
Float_t gAlgImbProj3C1_2_[5];
Float_t gAlgImbProj3C2_4_[5];
Float_t gAlgImbProj3C4_8_[5];
Float_t gAlgImbProj3C8_100_[5];

Float_t gAlgImbProj4CF_[5];
Float_t gAlgImbProj4C0_1_[5];
Float_t gAlgImbProj4C1_2_[5];
Float_t gAlgImbProj4C2_4_[5];
Float_t gAlgImbProj4C4_8_[5];
Float_t gAlgImbProj4C8_100_[5];

Float_t gAlgImbProj5CF_[5];
Float_t gAlgImbProj5C0_1_[5];
Float_t gAlgImbProj5C1_2_[5];
Float_t gAlgImbProj5C2_4_[5];
Float_t gAlgImbProj5C4_8_[5];
Float_t gAlgImbProj5C8_100_[5];

Float_t gAlgImbProj6CF_[5];
Float_t gAlgImbProj6C0_1_[5];
Float_t gAlgImbProj6C1_2_[5];
Float_t gAlgImbProj6C2_4_[5];
Float_t gAlgImbProj6C4_8_[5];
Float_t gAlgImbProj6C8_100_[5];

Float_t gAlgImbProj7CF_[5];
Float_t gAlgImbProj7C0_1_[5];
Float_t gAlgImbProj7C1_2_[5];
Float_t gAlgImbProj7C2_4_[5];
Float_t gAlgImbProj7C4_8_[5];
Float_t gAlgImbProj7C8_100_[5];

Float_t gAlgImbProj8CF_[5];
Float_t gAlgImbProj8C0_1_[5];
Float_t gAlgImbProj8C1_2_[5];
Float_t gAlgImbProj8C2_4_[5];
Float_t gAlgImbProj8C4_8_[5];
Float_t gAlgImbProj8C8_100_[5];

Float_t gAlgImbProj9CF_[5];
Float_t gAlgImbProj9C0_1_[5];
Float_t gAlgImbProj9C1_2_[5];
Float_t gAlgImbProj9C2_4_[5];
Float_t gAlgImbProj9C4_8_[5];
Float_t gAlgImbProj9C8_100_[5];

Float_t gAlgImbProj10CF_[5];
Float_t gAlgImbProj10C0_1_[5];
Float_t gAlgImbProj10C1_2_[5];
Float_t gAlgImbProj10C2_4_[5];
Float_t gAlgImbProj10C4_8_[5];
Float_t gAlgImbProj10C8_100_[5];

//Hem delR

Float_t gAlgImbHem1CF_[5];
Float_t gAlgImbHem1C0_1_[5];
Float_t gAlgImbHem1C1_2_[5];
Float_t gAlgImbHem1C2_4_[5];
Float_t gAlgImbHem1C4_8_[5];
Float_t gAlgImbHem1C8_100_[5];

Float_t gAlgImbHem2CF_[5];
Float_t gAlgImbHem2C0_1_[5];
Float_t gAlgImbHem2C1_2_[5];
Float_t gAlgImbHem2C2_4_[5];
Float_t gAlgImbHem2C4_8_[5];
Float_t gAlgImbHem2C8_100_[5];

Float_t gAlgImbHem3CF_[5];
Float_t gAlgImbHem3C0_1_[5];
Float_t gAlgImbHem3C1_2_[5];
Float_t gAlgImbHem3C2_4_[5];
Float_t gAlgImbHem3C4_8_[5];
Float_t gAlgImbHem3C8_100_[5];

Float_t gAlgImbHem4CF_[5];
Float_t gAlgImbHem4C0_1_[5];
Float_t gAlgImbHem4C1_2_[5];
Float_t gAlgImbHem4C2_4_[5];
Float_t gAlgImbHem4C4_8_[5];
Float_t gAlgImbHem4C8_100_[5];

Float_t gAlgImbHem5CF_[5];
Float_t gAlgImbHem5C0_1_[5];
Float_t gAlgImbHem5C1_2_[5];
Float_t gAlgImbHem5C2_4_[5];
Float_t gAlgImbHem5C4_8_[5];
Float_t gAlgImbHem5C8_100_[5];

Float_t gAlgImbHem6CF_[5];
Float_t gAlgImbHem6C0_1_[5];
Float_t gAlgImbHem6C1_2_[5];
Float_t gAlgImbHem6C2_4_[5];
Float_t gAlgImbHem6C4_8_[5];
Float_t gAlgImbHem6C8_100_[5];

Float_t gAlgImbHem7CF_[5];
Float_t gAlgImbHem7C0_1_[5];
Float_t gAlgImbHem7C1_2_[5];
Float_t gAlgImbHem7C2_4_[5];
Float_t gAlgImbHem7C4_8_[5];
Float_t gAlgImbHem7C8_100_[5];

Float_t gAlgImbHem8CF_[5];
Float_t gAlgImbHem8C0_1_[5];
Float_t gAlgImbHem8C1_2_[5];
Float_t gAlgImbHem8C2_4_[5];
Float_t gAlgImbHem8C4_8_[5];
Float_t gAlgImbHem8C8_100_[5];

Float_t gAlgImbHem9CF_[5];
Float_t gAlgImbHem9C0_1_[5];
Float_t gAlgImbHem9C1_2_[5];
Float_t gAlgImbHem9C2_4_[5];
Float_t gAlgImbHem9C4_8_[5];
Float_t gAlgImbHem9C8_100_[5];

Float_t gAlgImbHem10CF_[5];
Float_t gAlgImbHem10C0_1_[5];
Float_t gAlgImbHem10C1_2_[5];
Float_t gAlgImbHem10C2_4_[5];
Float_t gAlgImbHem10C4_8_[5];
Float_t gAlgImbHem10C8_100_[5];


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


//Truth delR plots


void SetBranches(bool montecarlo)
{
  //Track Tree Branches

  std::cout << "Branches Set" << std::endl;
  
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

  //Tracks proj. onto Alg, ordered according to enum above, All, Cone, and NotCone

  trackTree_p->Branch("rAlgImbProjF", &rAlgImbProjF_, "rAlgImbProjF[10]/F");
  trackTree_p->Branch("rAlgImbPerpF", &rAlgImbPerpF_, "rAlgImbPerpF[10]/F");
  trackTree_p->Branch("rAlgImbProj0_1", &rAlgImbProj0_1_, "rAlgImbProj0_1[10]/F");
  trackTree_p->Branch("rAlgImbProj1_2", &rAlgImbProj1_2_, "rAlgImbProj1_2[10]/F");
  trackTree_p->Branch("rAlgImbProj2_4", &rAlgImbProj2_4_, "rAlgImbProj2_4[10]/F");
  trackTree_p->Branch("rAlgImbProj4_8", &rAlgImbProj4_8_, "rAlgImbProj4_8[10]/F");
  trackTree_p->Branch("rAlgImbProj8_100", &rAlgImbProj8_100_, "rAlgImbProj8_100[10]/F");
  trackTree_p->Branch("rAlgImbProjCF", &rAlgImbProjCF_, "rAlgImbProjCF[10]/F");
  trackTree_p->Branch("rAlgImbPerpCF", &rAlgImbPerpCF_, "rAlgImbPerpCF[10]/F");
  trackTree_p->Branch("rAlgImbProjC0_1", &rAlgImbProjC0_1_, "rAlgImbProjC0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjC1_2", &rAlgImbProjC1_2_, "rAlgImbProjC1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjC2_4", &rAlgImbProjC2_4_, "rAlgImbProjC2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjC4_8", &rAlgImbProjC4_8_, "rAlgImbProjC4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjC8_100", &rAlgImbProjC8_100_, "rAlgImbProjC8_100[10]/F");
  trackTree_p->Branch("rAlgImbProjNCF", &rAlgImbProjNCF_, "rAlgImbProjNCF[10]/F");
  trackTree_p->Branch("rAlgImbPerpNCF", &rAlgImbPerpNCF_, "rAlgImbPerpNCF[10]/F");
  trackTree_p->Branch("rAlgImbProjNC0_1", &rAlgImbProjNC0_1_, "rAlgImbProjNC0_1[10]/F");
  trackTree_p->Branch("rAlgImbProjNC1_2", &rAlgImbProjNC1_2_, "rAlgImbProjNC1_2[10]/F");
  trackTree_p->Branch("rAlgImbProjNC2_4", &rAlgImbProjNC2_4_, "rAlgImbProjNC2_4[10]/F");
  trackTree_p->Branch("rAlgImbProjNC4_8", &rAlgImbProjNC4_8_, "rAlgImbProjNC4_8[10]/F");
  trackTree_p->Branch("rAlgImbProjNC8_100", &rAlgImbProjNC8_100_, "rAlgImbProjNC8_100[10]/F");

  trackTree_p->Branch("rAlgImbHemF", &rAlgImbHemF_, "rAlgImbHemF[10]/F");
  trackTree_p->Branch("rAlgImbHem0_1", &rAlgImbHem0_1_, "rAlgImbHem0_1[10]/F");
  trackTree_p->Branch("rAlgImbHem1_2", &rAlgImbHem1_2_, "rAlgImbHem1_2[10]/F");
  trackTree_p->Branch("rAlgImbHem2_4", &rAlgImbHem2_4_, "rAlgImbHem2_4[10]/F");
  trackTree_p->Branch("rAlgImbHem4_8", &rAlgImbHem4_8_, "rAlgImbHem4_8[10]/F");
  trackTree_p->Branch("rAlgImbHem8_100", &rAlgImbHem8_100_, "rAlgImbHem8_100[10]/F");
  trackTree_p->Branch("rAlgImbHemCF", &rAlgImbHemCF_, "rAlgImbHemCF[10]/F");
  trackTree_p->Branch("rAlgImbHemC0_1", &rAlgImbHemC0_1_, "rAlgImbHemC0_1[10]/F");
  trackTree_p->Branch("rAlgImbHemC1_2", &rAlgImbHemC1_2_, "rAlgImbHemC1_2[10]/F");
  trackTree_p->Branch("rAlgImbHemC2_4", &rAlgImbHemC2_4_, "rAlgImbHemC2_4[10]/F");
  trackTree_p->Branch("rAlgImbHemC4_8", &rAlgImbHemC4_8_, "rAlgImbHemC4_8[10]/F");
  trackTree_p->Branch("rAlgImbHemC8_100", &rAlgImbHemC8_100_, "rAlgImbHemC8_100[10]/F");
  trackTree_p->Branch("rAlgImbHemNCF", &rAlgImbHemNCF_, "rAlgImbHemNCF[10]/F");
  trackTree_p->Branch("rAlgImbHemNC0_1", &rAlgImbHemNC0_1_, "rAlgImbHemNC0_1[10]/F");
  trackTree_p->Branch("rAlgImbHemNC1_2", &rAlgImbHemNC1_2_, "rAlgImbHemNC1_2[10]/F");
  trackTree_p->Branch("rAlgImbHemNC2_4", &rAlgImbHemNC2_4_, "rAlgImbHemNC2_4[10]/F");
  trackTree_p->Branch("rAlgImbHemNC4_8", &rAlgImbHemNC4_8_, "rAlgImbHemNC4_8[10]/F");
  trackTree_p->Branch("rAlgImbHemNC8_100", &rAlgImbHemNC8_100_, "rAlgImbHemNC8_100[10]/F");


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

  //DelRs

  trackTree_p->Branch("rAlgImbProj1CF", &rAlgImbProj1CF_, "rAlgImbProj1CF[10]/F");
  trackTree_p->Branch("rAlgImbProj1C0_1", &rAlgImbProj1C0_1_, "rAlgImbProj1C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProj1C1_2", &rAlgImbProj1C1_2_, "rAlgImbProj1C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProj1C2_4", &rAlgImbProj1C2_4_, "rAlgImbProj1C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProj1C4_8", &rAlgImbProj1C4_8_, "rAlgImbProj1C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProj1C8_100", &rAlgImbProj1C8_100_, "rAlgImbProj1C8_100[10]/F");

  trackTree_p->Branch("rAlgImbProj2CF", &rAlgImbProj2CF_, "rAlgImbProj2CF[10]/F");
  trackTree_p->Branch("rAlgImbProj2C0_1", &rAlgImbProj2C0_1_, "rAlgImbProj2C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProj2C1_2", &rAlgImbProj2C1_2_, "rAlgImbProj2C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProj2C2_4", &rAlgImbProj2C2_4_, "rAlgImbProj2C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProj2C4_8", &rAlgImbProj2C4_8_, "rAlgImbProj2C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProj2C8_100", &rAlgImbProj2C8_100_, "rAlgImbProj2C8_100[10]/F");

  trackTree_p->Branch("rAlgImbProj3CF", &rAlgImbProj3CF_, "rAlgImbProj3CF[10]/F");
  trackTree_p->Branch("rAlgImbProj3C0_1", &rAlgImbProj3C0_1_, "rAlgImbProj3C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProj3C1_2", &rAlgImbProj3C1_2_, "rAlgImbProj3C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProj3C2_4", &rAlgImbProj3C2_4_, "rAlgImbProj3C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProj3C4_8", &rAlgImbProj3C4_8_, "rAlgImbProj3C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProj3C8_100", &rAlgImbProj3C8_100_, "rAlgImbProj3C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProj4CF", &rAlgImbProj4CF_, "rAlgImbProj4CF[10]/F");
  trackTree_p->Branch("rAlgImbProj4C0_1", &rAlgImbProj4C0_1_, "rAlgImbProj4C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProj4C1_2", &rAlgImbProj4C1_2_, "rAlgImbProj4C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProj4C2_4", &rAlgImbProj4C2_4_, "rAlgImbProj4C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProj4C4_8", &rAlgImbProj4C4_8_, "rAlgImbProj4C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProj4C8_100", &rAlgImbProj4C8_100_, "rAlgImbProj4C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProj5CF", &rAlgImbProj5CF_, "rAlgImbProj5CF[10]/F");
  trackTree_p->Branch("rAlgImbProj5C0_1", &rAlgImbProj5C0_1_, "rAlgImbProj5C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProj5C1_2", &rAlgImbProj5C1_2_, "rAlgImbProj5C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProj5C2_4", &rAlgImbProj5C2_4_, "rAlgImbProj5C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProj5C4_8", &rAlgImbProj5C4_8_, "rAlgImbProj5C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProj5C8_100", &rAlgImbProj5C8_100_, "rAlgImbProj5C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProj6CF", &rAlgImbProj6CF_, "rAlgImbProj6CF[10]/F");
  trackTree_p->Branch("rAlgImbProj6C0_1", &rAlgImbProj6C0_1_, "rAlgImbProj6C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProj6C1_2", &rAlgImbProj6C1_2_, "rAlgImbProj6C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProj6C2_4", &rAlgImbProj6C2_4_, "rAlgImbProj6C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProj6C4_8", &rAlgImbProj6C4_8_, "rAlgImbProj6C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProj6C8_100", &rAlgImbProj6C8_100_, "rAlgImbProj6C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProj7CF", &rAlgImbProj7CF_, "rAlgImbProj7CF[10]/F");
  trackTree_p->Branch("rAlgImbProj7C0_1", &rAlgImbProj7C0_1_, "rAlgImbProj7C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProj7C1_2", &rAlgImbProj7C1_2_, "rAlgImbProj7C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProj7C2_4", &rAlgImbProj7C2_4_, "rAlgImbProj7C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProj7C4_8", &rAlgImbProj7C4_8_, "rAlgImbProj7C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProj7C8_100", &rAlgImbProj7C8_100_, "rAlgImbProj7C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProj8CF", &rAlgImbProj8CF_, "rAlgImbProj8CF[10]/F");
  trackTree_p->Branch("rAlgImbProj8C0_1", &rAlgImbProj8C0_1_, "rAlgImbProj8C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProj8C1_2", &rAlgImbProj8C1_2_, "rAlgImbProj8C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProj8C2_4", &rAlgImbProj8C2_4_, "rAlgImbProj8C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProj8C4_8", &rAlgImbProj8C4_8_, "rAlgImbProj8C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProj8C8_100", &rAlgImbProj8C8_100_, "rAlgImbProj8C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProj9CF", &rAlgImbProj9CF_, "rAlgImbProj9CF[10]/F");
  trackTree_p->Branch("rAlgImbProj9C0_1", &rAlgImbProj9C0_1_, "rAlgImbProj9C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProj9C1_2", &rAlgImbProj9C1_2_, "rAlgImbProj9C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProj9C2_4", &rAlgImbProj9C2_4_, "rAlgImbProj9C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProj9C4_8", &rAlgImbProj9C4_8_, "rAlgImbProj9C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProj9C8_100", &rAlgImbProj9C8_100_, "rAlgImbProj9C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbProj10CF", &rAlgImbProj10CF_, "rAlgImbProj10CF[10]/F");
  trackTree_p->Branch("rAlgImbProj10C0_1", &rAlgImbProj10C0_1_, "rAlgImbProj10C0_1[10]/F");
  trackTree_p->Branch("rAlgImbProj10C1_2", &rAlgImbProj10C1_2_, "rAlgImbProj10C1_2[10]/F");
  trackTree_p->Branch("rAlgImbProj10C2_4", &rAlgImbProj10C2_4_, "rAlgImbProj10C2_4[10]/F");
  trackTree_p->Branch("rAlgImbProj10C4_8", &rAlgImbProj10C4_8_, "rAlgImbProj10C4_8[10]/F");
  trackTree_p->Branch("rAlgImbProj10C8_100", &rAlgImbProj10C8_100_, "rAlgImbProj10C8_100[10]/F");  

  //Hem DelRs

  trackTree_p->Branch("rAlgImbHem1CF", &rAlgImbHem1CF_, "rAlgImbHem1CF[10]/F");
  trackTree_p->Branch("rAlgImbHem1C0_1", &rAlgImbHem1C0_1_, "rAlgImbHem1C0_1[10]/F");
  trackTree_p->Branch("rAlgImbHem1C1_2", &rAlgImbHem1C1_2_, "rAlgImbHem1C1_2[10]/F");
  trackTree_p->Branch("rAlgImbHem1C2_4", &rAlgImbHem1C2_4_, "rAlgImbHem1C2_4[10]/F");
  trackTree_p->Branch("rAlgImbHem1C4_8", &rAlgImbHem1C4_8_, "rAlgImbHem1C4_8[10]/F");
  trackTree_p->Branch("rAlgImbHem1C8_100", &rAlgImbHem1C8_100_, "rAlgImbHem1C8_100[10]/F");

  trackTree_p->Branch("rAlgImbHem2CF", &rAlgImbHem2CF_, "rAlgImbHem2CF[10]/F");
  trackTree_p->Branch("rAlgImbHem2C0_1", &rAlgImbHem2C0_1_, "rAlgImbHem2C0_1[10]/F");
  trackTree_p->Branch("rAlgImbHem2C1_2", &rAlgImbHem2C1_2_, "rAlgImbHem2C1_2[10]/F");
  trackTree_p->Branch("rAlgImbHem2C2_4", &rAlgImbHem2C2_4_, "rAlgImbHem2C2_4[10]/F");
  trackTree_p->Branch("rAlgImbHem2C4_8", &rAlgImbHem2C4_8_, "rAlgImbHem2C4_8[10]/F");
  trackTree_p->Branch("rAlgImbHem2C8_100", &rAlgImbHem2C8_100_, "rAlgImbHem2C8_100[10]/F");

  trackTree_p->Branch("rAlgImbHem3CF", &rAlgImbHem3CF_, "rAlgImbHem3CF[10]/F");
  trackTree_p->Branch("rAlgImbHem3C0_1", &rAlgImbHem3C0_1_, "rAlgImbHem3C0_1[10]/F");
  trackTree_p->Branch("rAlgImbHem3C1_2", &rAlgImbHem3C1_2_, "rAlgImbHem3C1_2[10]/F");
  trackTree_p->Branch("rAlgImbHem3C2_4", &rAlgImbHem3C2_4_, "rAlgImbHem3C2_4[10]/F");
  trackTree_p->Branch("rAlgImbHem3C4_8", &rAlgImbHem3C4_8_, "rAlgImbHem3C4_8[10]/F");
  trackTree_p->Branch("rAlgImbHem3C8_100", &rAlgImbHem3C8_100_, "rAlgImbHem3C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbHem4CF", &rAlgImbHem4CF_, "rAlgImbHem4CF[10]/F");
  trackTree_p->Branch("rAlgImbHem4C0_1", &rAlgImbHem4C0_1_, "rAlgImbHem4C0_1[10]/F");
  trackTree_p->Branch("rAlgImbHem4C1_2", &rAlgImbHem4C1_2_, "rAlgImbHem4C1_2[10]/F");
  trackTree_p->Branch("rAlgImbHem4C2_4", &rAlgImbHem4C2_4_, "rAlgImbHem4C2_4[10]/F");
  trackTree_p->Branch("rAlgImbHem4C4_8", &rAlgImbHem4C4_8_, "rAlgImbHem4C4_8[10]/F");
  trackTree_p->Branch("rAlgImbHem4C8_100", &rAlgImbHem4C8_100_, "rAlgImbHem4C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbHem5CF", &rAlgImbHem5CF_, "rAlgImbHem5CF[10]/F");
  trackTree_p->Branch("rAlgImbHem5C0_1", &rAlgImbHem5C0_1_, "rAlgImbHem5C0_1[10]/F");
  trackTree_p->Branch("rAlgImbHem5C1_2", &rAlgImbHem5C1_2_, "rAlgImbHem5C1_2[10]/F");
  trackTree_p->Branch("rAlgImbHem5C2_4", &rAlgImbHem5C2_4_, "rAlgImbHem5C2_4[10]/F");
  trackTree_p->Branch("rAlgImbHem5C4_8", &rAlgImbHem5C4_8_, "rAlgImbHem5C4_8[10]/F");
  trackTree_p->Branch("rAlgImbHem5C8_100", &rAlgImbHem5C8_100_, "rAlgImbHem5C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbHem6CF", &rAlgImbHem6CF_, "rAlgImbHem6CF[10]/F");
  trackTree_p->Branch("rAlgImbHem6C0_1", &rAlgImbHem6C0_1_, "rAlgImbHem6C0_1[10]/F");
  trackTree_p->Branch("rAlgImbHem6C1_2", &rAlgImbHem6C1_2_, "rAlgImbHem6C1_2[10]/F");
  trackTree_p->Branch("rAlgImbHem6C2_4", &rAlgImbHem6C2_4_, "rAlgImbHem6C2_4[10]/F");
  trackTree_p->Branch("rAlgImbHem6C4_8", &rAlgImbHem6C4_8_, "rAlgImbHem6C4_8[10]/F");
  trackTree_p->Branch("rAlgImbHem6C8_100", &rAlgImbHem6C8_100_, "rAlgImbHem6C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbHem7CF", &rAlgImbHem7CF_, "rAlgImbHem7CF[10]/F");
  trackTree_p->Branch("rAlgImbHem7C0_1", &rAlgImbHem7C0_1_, "rAlgImbHem7C0_1[10]/F");
  trackTree_p->Branch("rAlgImbHem7C1_2", &rAlgImbHem7C1_2_, "rAlgImbHem7C1_2[10]/F");
  trackTree_p->Branch("rAlgImbHem7C2_4", &rAlgImbHem7C2_4_, "rAlgImbHem7C2_4[10]/F");
  trackTree_p->Branch("rAlgImbHem7C4_8", &rAlgImbHem7C4_8_, "rAlgImbHem7C4_8[10]/F");
  trackTree_p->Branch("rAlgImbHem7C8_100", &rAlgImbHem7C8_100_, "rAlgImbHem7C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbHem8CF", &rAlgImbHem8CF_, "rAlgImbHem8CF[10]/F");
  trackTree_p->Branch("rAlgImbHem8C0_1", &rAlgImbHem8C0_1_, "rAlgImbHem8C0_1[10]/F");
  trackTree_p->Branch("rAlgImbHem8C1_2", &rAlgImbHem8C1_2_, "rAlgImbHem8C1_2[10]/F");
  trackTree_p->Branch("rAlgImbHem8C2_4", &rAlgImbHem8C2_4_, "rAlgImbHem8C2_4[10]/F");
  trackTree_p->Branch("rAlgImbHem8C4_8", &rAlgImbHem8C4_8_, "rAlgImbHem8C4_8[10]/F");
  trackTree_p->Branch("rAlgImbHem8C8_100", &rAlgImbHem8C8_100_, "rAlgImbHem8C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbHem9CF", &rAlgImbHem9CF_, "rAlgImbHem9CF[10]/F");
  trackTree_p->Branch("rAlgImbHem9C0_1", &rAlgImbHem9C0_1_, "rAlgImbHem9C0_1[10]/F");
  trackTree_p->Branch("rAlgImbHem9C1_2", &rAlgImbHem9C1_2_, "rAlgImbHem9C1_2[10]/F");
  trackTree_p->Branch("rAlgImbHem9C2_4", &rAlgImbHem9C2_4_, "rAlgImbHem9C2_4[10]/F");
  trackTree_p->Branch("rAlgImbHem9C4_8", &rAlgImbHem9C4_8_, "rAlgImbHem9C4_8[10]/F");
  trackTree_p->Branch("rAlgImbHem9C8_100", &rAlgImbHem9C8_100_, "rAlgImbHem9C8_100[10]/F");  

  trackTree_p->Branch("rAlgImbHem10CF", &rAlgImbHem10CF_, "rAlgImbHem10CF[10]/F");
  trackTree_p->Branch("rAlgImbHem10C0_1", &rAlgImbHem10C0_1_, "rAlgImbHem10C0_1[10]/F");
  trackTree_p->Branch("rAlgImbHem10C1_2", &rAlgImbHem10C1_2_, "rAlgImbHem10C1_2[10]/F");
  trackTree_p->Branch("rAlgImbHem10C2_4", &rAlgImbHem10C2_4_, "rAlgImbHem10C2_4[10]/F");
  trackTree_p->Branch("rAlgImbHem10C4_8", &rAlgImbHem10C4_8_, "rAlgImbHem10C4_8[10]/F");
  trackTree_p->Branch("rAlgImbHem10C8_100", &rAlgImbHem10C8_100_, "rAlgImbHem10C8_100[10]/F");  

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

  jetTree_p->Branch("eventSet", &eventSet_, "eventSet[5]/O");
  jetTree_p->Branch("setWeight", &setWeight_, "setWeight[5]/F");
  jetTree_p->Branch("setWeight_2pi3", &setWeight_2pi3_, "setWeight_2pi3[5]/F");

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

    //Gen. proj. onto jetAlg, array ordered according to enum

    genTree_p->Branch("gAlgImbProjF", &gAlgImbProjF_, "gAlgImbProjF[5]/F");
    genTree_p->Branch("gAlgImbPerpF", &gAlgImbPerpF_, "gAlgImbPerpF[5]/F");
    genTree_p->Branch("gAlgImbProj0_1", &gAlgImbProj0_1_, "gAlgImbProj0_1[5]/F");
    genTree_p->Branch("gAlgImbProj1_2", &gAlgImbProj1_2_, "gAlgImbProj1_2[5]/F");
    genTree_p->Branch("gAlgImbProj2_4", &gAlgImbProj2_4_, "gAlgImbProj2_4[5]/F");
    genTree_p->Branch("gAlgImbProj4_8", &gAlgImbProj4_8_, "gAlgImbProj4_8[5]/F");
    genTree_p->Branch("gAlgImbProj8_100", &gAlgImbProj8_100_, "gAlgImbProj8_100[5]/F");
    genTree_p->Branch("gAlgImbProjCF", &gAlgImbProjCF_, "gAlgImbProjCF[5]/F");
    genTree_p->Branch("gAlgImbPerpCF", &gAlgImbPerpCF_, "gAlgImbPerpCF[5]/F");
    genTree_p->Branch("gAlgImbProjC0_1", &gAlgImbProjC0_1_, "gAlgImbProjC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjC1_2", &gAlgImbProjC1_2_, "gAlgImbProjC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjC2_4", &gAlgImbProjC2_4_, "gAlgImbProjC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjC4_8", &gAlgImbProjC4_8_, "gAlgImbProjC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjC8_100", &gAlgImbProjC8_100_, "gAlgImbProjC8_100[5]/F");
    genTree_p->Branch("gAlgImbProjNCF", &gAlgImbProjNCF_, "gAlgImbProjNCF[5]/F");
    genTree_p->Branch("gAlgImbPerpNCF", &gAlgImbPerpNCF_, "gAlgImbPerpNCF[5]/F");
    genTree_p->Branch("gAlgImbProjNC0_1", &gAlgImbProjNC0_1_, "gAlgImbProjNC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjNC1_2", &gAlgImbProjNC1_2_, "gAlgImbProjNC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjNC2_4", &gAlgImbProjNC2_4_, "gAlgImbProjNC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjNC4_8", &gAlgImbProjNC4_8_, "gAlgImbProjNC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjNC8_100", &gAlgImbProjNC8_100_, "gAlgImbProjNC8_100[5]/F");

    genTree_p->Branch("gAlgImbProjSigCF", &gAlgImbProjSigCF_, "gAlgImbProjSigCF[5]/F");
    genTree_p->Branch("gAlgImbPerpSigCF", &gAlgImbPerpSigCF_, "gAlgImbPerpSigCF[5]/F");
    genTree_p->Branch("gAlgImbProjSigC0_1", &gAlgImbProjSigC0_1_, "gAlgImbProjSigC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjSigC1_2", &gAlgImbProjSigC1_2_, "gAlgImbProjSigC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjSigC2_4", &gAlgImbProjSigC2_4_, "gAlgImbProjSigC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjSigC4_8", &gAlgImbProjSigC4_8_, "gAlgImbProjSigC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjSigC8_100", &gAlgImbProjSigC8_100_, "gAlgImbProjSigC8_100[5]/F");
    genTree_p->Branch("gAlgImbProjSigNCF", &gAlgImbProjSigNCF_, "gAlgImbProjSigNCF[5]/F");
    genTree_p->Branch("gAlgImbPerpSigNCF", &gAlgImbPerpSigNCF_, "gAlgImbPerpSigNCF[5]/F");
    genTree_p->Branch("gAlgImbProjSigNC0_1", &gAlgImbProjSigNC0_1_, "gAlgImbProjSigNC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjSigNC1_2", &gAlgImbProjSigNC1_2_, "gAlgImbProjSigNC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjSigNC2_4", &gAlgImbProjSigNC2_4_, "gAlgImbProjSigNC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjSigNC4_8", &gAlgImbProjSigNC4_8_, "gAlgImbProjSigNC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjSigNC8_100", &gAlgImbProjSigNC8_100_, "gAlgImbProjSigNC8_100[5]/F");

    genTree_p->Branch("gAlgImbProjUECF", &gAlgImbProjUECF_, "gAlgImbProjUECF[5]/F");
    genTree_p->Branch("gAlgImbPerpUECF", &gAlgImbPerpUECF_, "gAlgImbPerpUECF[5]/F");
    genTree_p->Branch("gAlgImbProjUEC0_1", &gAlgImbProjUEC0_1_, "gAlgImbProjUEC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjUEC1_2", &gAlgImbProjUEC1_2_, "gAlgImbProjUEC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjUEC2_4", &gAlgImbProjUEC2_4_, "gAlgImbProjUEC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjUEC4_8", &gAlgImbProjUEC4_8_, "gAlgImbProjUEC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjUEC8_100", &gAlgImbProjUEC8_100_, "gAlgImbProjUEC8_100[5]/F");
    genTree_p->Branch("gAlgImbProjUENCF", &gAlgImbProjUENCF_, "gAlgImbProjUENCF[5]/F");
    genTree_p->Branch("gAlgImbPerpUENCF", &gAlgImbPerpUENCF_, "gAlgImbPerpUENCF[5]/F");
    genTree_p->Branch("gAlgImbProjUENC0_1", &gAlgImbProjUENC0_1_, "gAlgImbProjUENC0_1[5]/F");
    genTree_p->Branch("gAlgImbProjUENC1_2", &gAlgImbProjUENC1_2_, "gAlgImbProjUENC1_2[5]/F");
    genTree_p->Branch("gAlgImbProjUENC2_4", &gAlgImbProjUENC2_4_, "gAlgImbProjUENC2_4[5]/F");
    genTree_p->Branch("gAlgImbProjUENC4_8", &gAlgImbProjUENC4_8_, "gAlgImbProjUENC4_8[5]/F");
    genTree_p->Branch("gAlgImbProjUENC8_100", &gAlgImbProjUENC8_100_, "gAlgImbProjUENC8_100[5]/F");

    genTree_p->Branch("gAlgImbHemF", &gAlgImbHemF_, "gAlgImbHemF[5]/F");
    genTree_p->Branch("gAlgImbHem0_1", &gAlgImbHem0_1_, "gAlgImbHem0_1[5]/F");
    genTree_p->Branch("gAlgImbHem1_2", &gAlgImbHem1_2_, "gAlgImbHem1_2[5]/F");
    genTree_p->Branch("gAlgImbHem2_4", &gAlgImbHem2_4_, "gAlgImbHem2_4[5]/F");
    genTree_p->Branch("gAlgImbHem4_8", &gAlgImbHem4_8_, "gAlgImbHem4_8[5]/F");
    genTree_p->Branch("gAlgImbHem8_100", &gAlgImbHem8_100_, "gAlgImbHem8_100[5]/F");
    genTree_p->Branch("gAlgImbHemCF", &gAlgImbHemCF_, "gAlgImbHemCF[5]/F");
    genTree_p->Branch("gAlgImbHemC0_1", &gAlgImbHemC0_1_, "gAlgImbHemC0_1[5]/F");
    genTree_p->Branch("gAlgImbHemC1_2", &gAlgImbHemC1_2_, "gAlgImbHemC1_2[5]/F");
    genTree_p->Branch("gAlgImbHemC2_4", &gAlgImbHemC2_4_, "gAlgImbHemC2_4[5]/F");
    genTree_p->Branch("gAlgImbHemC4_8", &gAlgImbHemC4_8_, "gAlgImbHemC4_8[5]/F");
    genTree_p->Branch("gAlgImbHemC8_100", &gAlgImbHemC8_100_, "gAlgImbHemC8_100[5]/F");
    genTree_p->Branch("gAlgImbHemNCF", &gAlgImbHemNCF_, "gAlgImbHemNCF[5]/F");
    genTree_p->Branch("gAlgImbHemNC0_1", &gAlgImbHemNC0_1_, "gAlgImbHemNC0_1[5]/F");
    genTree_p->Branch("gAlgImbHemNC1_2", &gAlgImbHemNC1_2_, "gAlgImbHemNC1_2[5]/F");
    genTree_p->Branch("gAlgImbHemNC2_4", &gAlgImbHemNC2_4_, "gAlgImbHemNC2_4[5]/F");
    genTree_p->Branch("gAlgImbHemNC4_8", &gAlgImbHemNC4_8_, "gAlgImbHemNC4_8[5]/F");
    genTree_p->Branch("gAlgImbHemNC8_100", &gAlgImbHemNC8_100_, "gAlgImbHemNC8_100[5]/F");

    genTree_p->Branch("gAlgImbHemSigCF", &gAlgImbHemSigCF_, "gAlgImbHemSigCF[5]/F");
    genTree_p->Branch("gAlgImbHemSigC0_1", &gAlgImbHemSigC0_1_, "gAlgImbHemSigC0_1[5]/F");
    genTree_p->Branch("gAlgImbHemSigC1_2", &gAlgImbHemSigC1_2_, "gAlgImbHemSigC1_2[5]/F");
    genTree_p->Branch("gAlgImbHemSigC2_4", &gAlgImbHemSigC2_4_, "gAlgImbHemSigC2_4[5]/F");
    genTree_p->Branch("gAlgImbHemSigC4_8", &gAlgImbHemSigC4_8_, "gAlgImbHemSigC4_8[5]/F");
    genTree_p->Branch("gAlgImbHemSigC8_100", &gAlgImbHemSigC8_100_, "gAlgImbHemSigC8_100[5]/F");
    genTree_p->Branch("gAlgImbHemSigNCF", &gAlgImbHemSigNCF_, "gAlgImbHemSigNCF[5]/F");
    genTree_p->Branch("gAlgImbHemSigNC0_1", &gAlgImbHemSigNC0_1_, "gAlgImbHemSigNC0_1[5]/F");
    genTree_p->Branch("gAlgImbHemSigNC1_2", &gAlgImbHemSigNC1_2_, "gAlgImbHemSigNC1_2[5]/F");
    genTree_p->Branch("gAlgImbHemSigNC2_4", &gAlgImbHemSigNC2_4_, "gAlgImbHemSigNC2_4[5]/F");
    genTree_p->Branch("gAlgImbHemSigNC4_8", &gAlgImbHemSigNC4_8_, "gAlgImbHemSigNC4_8[5]/F");
    genTree_p->Branch("gAlgImbHemSigNC8_100", &gAlgImbHemSigNC8_100_, "gAlgImbHemSigNC8_100[5]/F");

    genTree_p->Branch("gAlgImbHemUECF", &gAlgImbHemUECF_, "gAlgImbHemUECF[5]/F");
    genTree_p->Branch("gAlgImbHemUEC0_1", &gAlgImbHemUEC0_1_, "gAlgImbHemUEC0_1[5]/F");
    genTree_p->Branch("gAlgImbHemUEC1_2", &gAlgImbHemUEC1_2_, "gAlgImbHemUEC1_2[5]/F");
    genTree_p->Branch("gAlgImbHemUEC2_4", &gAlgImbHemUEC2_4_, "gAlgImbHemUEC2_4[5]/F");
    genTree_p->Branch("gAlgImbHemUEC4_8", &gAlgImbHemUEC4_8_, "gAlgImbHemUEC4_8[5]/F");
    genTree_p->Branch("gAlgImbHemUEC8_100", &gAlgImbHemUEC8_100_, "gAlgImbHemUEC8_100[5]/F");
    genTree_p->Branch("gAlgImbHemUENCF", &gAlgImbHemUENCF_, "gAlgImbHemUENCF[5]/F");
    genTree_p->Branch("gAlgImbHemUENC0_1", &gAlgImbHemUENC0_1_, "gAlgImbHemUENC0_1[5]/F");
    genTree_p->Branch("gAlgImbHemUENC1_2", &gAlgImbHemUENC1_2_, "gAlgImbHemUENC1_2[5]/F");
    genTree_p->Branch("gAlgImbHemUENC2_4", &gAlgImbHemUENC2_4_, "gAlgImbHemUENC2_4[5]/F");
    genTree_p->Branch("gAlgImbHemUENC4_8", &gAlgImbHemUENC4_8_, "gAlgImbHemUENC4_8[5]/F");
    genTree_p->Branch("gAlgImbHemUENC8_100", &gAlgImbHemUENC8_100_, "gAlgImbHemUENC8_100[5]/F");

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

    genTree_p->Branch("gAlgImbProj1CF", &gAlgImbProj1CF_, "gAlgImbProj1CF[5]/F");
    genTree_p->Branch("gAlgImbProj1C0_1", &gAlgImbProj1C0_1_, "gAlgImbProj1C0_1[5]/F");
    genTree_p->Branch("gAlgImbProj1C1_2", &gAlgImbProj1C1_2_, "gAlgImbProj1C1_2[5]/F");
    genTree_p->Branch("gAlgImbProj1C2_4", &gAlgImbProj1C2_4_, "gAlgImbProj1C2_4[5]/F");
    genTree_p->Branch("gAlgImbProj1C4_8", &gAlgImbProj1C4_8_, "gAlgImbProj1C4_8[5]/F");
    genTree_p->Branch("gAlgImbProj1C8_100", &gAlgImbProj1C8_100_, "gAlgImbProj1C8_100[5]/F");

    genTree_p->Branch("gAlgImbProj2CF", &gAlgImbProj2CF_, "gAlgImbProj2CF[5]/F");
    genTree_p->Branch("gAlgImbProj2C0_1", &gAlgImbProj2C0_1_, "gAlgImbProj2C0_1[5]/F");
    genTree_p->Branch("gAlgImbProj2C1_2", &gAlgImbProj2C1_2_, "gAlgImbProj2C1_2[5]/F");
    genTree_p->Branch("gAlgImbProj2C2_4", &gAlgImbProj2C2_4_, "gAlgImbProj2C2_4[5]/F");
    genTree_p->Branch("gAlgImbProj2C4_8", &gAlgImbProj2C4_8_, "gAlgImbProj2C4_8[5]/F");
    genTree_p->Branch("gAlgImbProj2C8_100", &gAlgImbProj2C8_100_, "gAlgImbProj2C8_100[5]/F");

    genTree_p->Branch("gAlgImbProj3CF", &gAlgImbProj3CF_, "gAlgImbProj3CF[5]/F");
    genTree_p->Branch("gAlgImbProj3C0_1", &gAlgImbProj3C0_1_, "gAlgImbProj3C0_1[5]/F");
    genTree_p->Branch("gAlgImbProj3C1_2", &gAlgImbProj3C1_2_, "gAlgImbProj3C1_2[5]/F");
    genTree_p->Branch("gAlgImbProj3C2_4", &gAlgImbProj3C2_4_, "gAlgImbProj3C2_4[5]/F");
    genTree_p->Branch("gAlgImbProj3C4_8", &gAlgImbProj3C4_8_, "gAlgImbProj3C4_8[5]/F");
    genTree_p->Branch("gAlgImbProj3C8_100", &gAlgImbProj3C8_100_, "gAlgImbProj3C8_100[5]/F");

    genTree_p->Branch("gAlgImbProj4CF", &gAlgImbProj4CF_, "gAlgImbProj4CF[5]/F");
    genTree_p->Branch("gAlgImbProj4C0_1", &gAlgImbProj4C0_1_, "gAlgImbProj4C0_1[5]/F");
    genTree_p->Branch("gAlgImbProj4C1_2", &gAlgImbProj4C1_2_, "gAlgImbProj4C1_2[5]/F");
    genTree_p->Branch("gAlgImbProj4C2_4", &gAlgImbProj4C2_4_, "gAlgImbProj4C2_4[5]/F");
    genTree_p->Branch("gAlgImbProj4C4_8", &gAlgImbProj4C4_8_, "gAlgImbProj4C4_8[5]/F");
    genTree_p->Branch("gAlgImbProj4C8_100", &gAlgImbProj4C8_100_, "gAlgImbProj4C8_100[5]/F");

    genTree_p->Branch("gAlgImbProj5CF", &gAlgImbProj5CF_, "gAlgImbProj5CF[5]/F");
    genTree_p->Branch("gAlgImbProj5C0_1", &gAlgImbProj5C0_1_, "gAlgImbProj5C0_1[5]/F");
    genTree_p->Branch("gAlgImbProj5C1_2", &gAlgImbProj5C1_2_, "gAlgImbProj5C1_2[5]/F");
    genTree_p->Branch("gAlgImbProj5C2_4", &gAlgImbProj5C2_4_, "gAlgImbProj5C2_4[5]/F");
    genTree_p->Branch("gAlgImbProj5C4_8", &gAlgImbProj5C4_8_, "gAlgImbProj5C4_8[5]/F");
    genTree_p->Branch("gAlgImbProj5C8_100", &gAlgImbProj5C8_100_, "gAlgImbProj5C8_100[5]/F");

    genTree_p->Branch("gAlgImbProj6CF", &gAlgImbProj6CF_, "gAlgImbProj6CF[5]/F");
    genTree_p->Branch("gAlgImbProj6C0_1", &gAlgImbProj6C0_1_, "gAlgImbProj6C0_1[5]/F");
    genTree_p->Branch("gAlgImbProj6C1_2", &gAlgImbProj6C1_2_, "gAlgImbProj6C1_2[5]/F");
    genTree_p->Branch("gAlgImbProj6C2_4", &gAlgImbProj6C2_4_, "gAlgImbProj6C2_4[5]/F");
    genTree_p->Branch("gAlgImbProj6C4_8", &gAlgImbProj6C4_8_, "gAlgImbProj6C4_8[5]/F");
    genTree_p->Branch("gAlgImbProj6C8_100", &gAlgImbProj6C8_100_, "gAlgImbProj6C8_100[5]/F");

    genTree_p->Branch("gAlgImbProj7CF", &gAlgImbProj7CF_, "gAlgImbProj7CF[5]/F");
    genTree_p->Branch("gAlgImbProj7C0_1", &gAlgImbProj7C0_1_, "gAlgImbProj7C0_1[5]/F");
    genTree_p->Branch("gAlgImbProj7C1_2", &gAlgImbProj7C1_2_, "gAlgImbProj7C1_2[5]/F");
    genTree_p->Branch("gAlgImbProj7C2_4", &gAlgImbProj7C2_4_, "gAlgImbProj7C2_4[5]/F");
    genTree_p->Branch("gAlgImbProj7C4_8", &gAlgImbProj7C4_8_, "gAlgImbProj7C4_8[5]/F");
    genTree_p->Branch("gAlgImbProj7C8_100", &gAlgImbProj7C8_100_, "gAlgImbProj7C8_100[5]/F");

    genTree_p->Branch("gAlgImbProj8CF", &gAlgImbProj8CF_, "gAlgImbProj8CF[5]/F");
    genTree_p->Branch("gAlgImbProj8C0_1", &gAlgImbProj8C0_1_, "gAlgImbProj8C0_1[5]/F");
    genTree_p->Branch("gAlgImbProj8C1_2", &gAlgImbProj8C1_2_, "gAlgImbProj8C1_2[5]/F");
    genTree_p->Branch("gAlgImbProj8C2_4", &gAlgImbProj8C2_4_, "gAlgImbProj8C2_4[5]/F");
    genTree_p->Branch("gAlgImbProj8C4_8", &gAlgImbProj8C4_8_, "gAlgImbProj8C4_8[5]/F");
    genTree_p->Branch("gAlgImbProj8C8_100", &gAlgImbProj8C8_100_, "gAlgImbProj8C8_100[5]/F");

    genTree_p->Branch("gAlgImbProj9CF", &gAlgImbProj9CF_, "gAlgImbProj9CF[5]/F");
    genTree_p->Branch("gAlgImbProj9C0_1", &gAlgImbProj9C0_1_, "gAlgImbProj9C0_1[5]/F");
    genTree_p->Branch("gAlgImbProj9C1_2", &gAlgImbProj9C1_2_, "gAlgImbProj9C1_2[5]/F");
    genTree_p->Branch("gAlgImbProj9C2_4", &gAlgImbProj9C2_4_, "gAlgImbProj9C2_4[5]/F");
    genTree_p->Branch("gAlgImbProj9C4_8", &gAlgImbProj9C4_8_, "gAlgImbProj9C4_8[5]/F");
    genTree_p->Branch("gAlgImbProj9C8_100", &gAlgImbProj9C8_100_, "gAlgImbProj9C8_100[5]/F");

    genTree_p->Branch("gAlgImbProj10CF", &gAlgImbProj10CF_, "gAlgImbProj10CF[5]/F");
    genTree_p->Branch("gAlgImbProj10C0_1", &gAlgImbProj10C0_1_, "gAlgImbProj10C0_1[5]/F");
    genTree_p->Branch("gAlgImbProj10C1_2", &gAlgImbProj10C1_2_, "gAlgImbProj10C1_2[5]/F");
    genTree_p->Branch("gAlgImbProj10C2_4", &gAlgImbProj10C2_4_, "gAlgImbProj10C2_4[5]/F");
    genTree_p->Branch("gAlgImbProj10C4_8", &gAlgImbProj10C4_8_, "gAlgImbProj10C4_8[5]/F");
    genTree_p->Branch("gAlgImbProj10C8_100", &gAlgImbProj10C8_100_, "gAlgImbProj10C8_100[5]/F");

    //Hems

    genTree_p->Branch("gAlgImbHem1CF", &gAlgImbHem1CF_, "gAlgImbHem1CF[5]/F");
    genTree_p->Branch("gAlgImbHem1C0_1", &gAlgImbHem1C0_1_, "gAlgImbHem1C0_1[5]/F");
    genTree_p->Branch("gAlgImbHem1C1_2", &gAlgImbHem1C1_2_, "gAlgImbHem1C1_2[5]/F");
    genTree_p->Branch("gAlgImbHem1C2_4", &gAlgImbHem1C2_4_, "gAlgImbHem1C2_4[5]/F");
    genTree_p->Branch("gAlgImbHem1C4_8", &gAlgImbHem1C4_8_, "gAlgImbHem1C4_8[5]/F");
    genTree_p->Branch("gAlgImbHem1C8_100", &gAlgImbHem1C8_100_, "gAlgImbHem1C8_100[5]/F");

    genTree_p->Branch("gAlgImbHem2CF", &gAlgImbHem2CF_, "gAlgImbHem2CF[5]/F");
    genTree_p->Branch("gAlgImbHem2C0_1", &gAlgImbHem2C0_1_, "gAlgImbHem2C0_1[5]/F");
    genTree_p->Branch("gAlgImbHem2C1_2", &gAlgImbHem2C1_2_, "gAlgImbHem2C1_2[5]/F");
    genTree_p->Branch("gAlgImbHem2C2_4", &gAlgImbHem2C2_4_, "gAlgImbHem2C2_4[5]/F");
    genTree_p->Branch("gAlgImbHem2C4_8", &gAlgImbHem2C4_8_, "gAlgImbHem2C4_8[5]/F");
    genTree_p->Branch("gAlgImbHem2C8_100", &gAlgImbHem2C8_100_, "gAlgImbHem2C8_100[5]/F");

    genTree_p->Branch("gAlgImbHem3CF", &gAlgImbHem3CF_, "gAlgImbHem3CF[5]/F");
    genTree_p->Branch("gAlgImbHem3C0_1", &gAlgImbHem3C0_1_, "gAlgImbHem3C0_1[5]/F");
    genTree_p->Branch("gAlgImbHem3C1_2", &gAlgImbHem3C1_2_, "gAlgImbHem3C1_2[5]/F");
    genTree_p->Branch("gAlgImbHem3C2_4", &gAlgImbHem3C2_4_, "gAlgImbHem3C2_4[5]/F");
    genTree_p->Branch("gAlgImbHem3C4_8", &gAlgImbHem3C4_8_, "gAlgImbHem3C4_8[5]/F");
    genTree_p->Branch("gAlgImbHem3C8_100", &gAlgImbHem3C8_100_, "gAlgImbHem3C8_100[5]/F");

    genTree_p->Branch("gAlgImbHem4CF", &gAlgImbHem4CF_, "gAlgImbHem4CF[5]/F");
    genTree_p->Branch("gAlgImbHem4C0_1", &gAlgImbHem4C0_1_, "gAlgImbHem4C0_1[5]/F");
    genTree_p->Branch("gAlgImbHem4C1_2", &gAlgImbHem4C1_2_, "gAlgImbHem4C1_2[5]/F");
    genTree_p->Branch("gAlgImbHem4C2_4", &gAlgImbHem4C2_4_, "gAlgImbHem4C2_4[5]/F");
    genTree_p->Branch("gAlgImbHem4C4_8", &gAlgImbHem4C4_8_, "gAlgImbHem4C4_8[5]/F");
    genTree_p->Branch("gAlgImbHem4C8_100", &gAlgImbHem4C8_100_, "gAlgImbHem4C8_100[5]/F");

    genTree_p->Branch("gAlgImbHem5CF", &gAlgImbHem5CF_, "gAlgImbHem5CF[5]/F");
    genTree_p->Branch("gAlgImbHem5C0_1", &gAlgImbHem5C0_1_, "gAlgImbHem5C0_1[5]/F");
    genTree_p->Branch("gAlgImbHem5C1_2", &gAlgImbHem5C1_2_, "gAlgImbHem5C1_2[5]/F");
    genTree_p->Branch("gAlgImbHem5C2_4", &gAlgImbHem5C2_4_, "gAlgImbHem5C2_4[5]/F");
    genTree_p->Branch("gAlgImbHem5C4_8", &gAlgImbHem5C4_8_, "gAlgImbHem5C4_8[5]/F");
    genTree_p->Branch("gAlgImbHem5C8_100", &gAlgImbHem5C8_100_, "gAlgImbHem5C8_100[5]/F");

    genTree_p->Branch("gAlgImbHem6CF", &gAlgImbHem6CF_, "gAlgImbHem6CF[5]/F");
    genTree_p->Branch("gAlgImbHem6C0_1", &gAlgImbHem6C0_1_, "gAlgImbHem6C0_1[5]/F");
    genTree_p->Branch("gAlgImbHem6C1_2", &gAlgImbHem6C1_2_, "gAlgImbHem6C1_2[5]/F");
    genTree_p->Branch("gAlgImbHem6C2_4", &gAlgImbHem6C2_4_, "gAlgImbHem6C2_4[5]/F");
    genTree_p->Branch("gAlgImbHem6C4_8", &gAlgImbHem6C4_8_, "gAlgImbHem6C4_8[5]/F");
    genTree_p->Branch("gAlgImbHem6C8_100", &gAlgImbHem6C8_100_, "gAlgImbHem6C8_100[5]/F");

    genTree_p->Branch("gAlgImbHem7CF", &gAlgImbHem7CF_, "gAlgImbHem7CF[5]/F");
    genTree_p->Branch("gAlgImbHem7C0_1", &gAlgImbHem7C0_1_, "gAlgImbHem7C0_1[5]/F");
    genTree_p->Branch("gAlgImbHem7C1_2", &gAlgImbHem7C1_2_, "gAlgImbHem7C1_2[5]/F");
    genTree_p->Branch("gAlgImbHem7C2_4", &gAlgImbHem7C2_4_, "gAlgImbHem7C2_4[5]/F");
    genTree_p->Branch("gAlgImbHem7C4_8", &gAlgImbHem7C4_8_, "gAlgImbHem7C4_8[5]/F");
    genTree_p->Branch("gAlgImbHem7C8_100", &gAlgImbHem7C8_100_, "gAlgImbHem7C8_100[5]/F");

    genTree_p->Branch("gAlgImbHem8CF", &gAlgImbHem8CF_, "gAlgImbHem8CF[5]/F");
    genTree_p->Branch("gAlgImbHem8C0_1", &gAlgImbHem8C0_1_, "gAlgImbHem8C0_1[5]/F");
    genTree_p->Branch("gAlgImbHem8C1_2", &gAlgImbHem8C1_2_, "gAlgImbHem8C1_2[5]/F");
    genTree_p->Branch("gAlgImbHem8C2_4", &gAlgImbHem8C2_4_, "gAlgImbHem8C2_4[5]/F");
    genTree_p->Branch("gAlgImbHem8C4_8", &gAlgImbHem8C4_8_, "gAlgImbHem8C4_8[5]/F");
    genTree_p->Branch("gAlgImbHem8C8_100", &gAlgImbHem8C8_100_, "gAlgImbHem8C8_100[5]/F");

    genTree_p->Branch("gAlgImbHem9CF", &gAlgImbHem9CF_, "gAlgImbHem9CF[5]/F");
    genTree_p->Branch("gAlgImbHem9C0_1", &gAlgImbHem9C0_1_, "gAlgImbHem9C0_1[5]/F");
    genTree_p->Branch("gAlgImbHem9C1_2", &gAlgImbHem9C1_2_, "gAlgImbHem9C1_2[5]/F");
    genTree_p->Branch("gAlgImbHem9C2_4", &gAlgImbHem9C2_4_, "gAlgImbHem9C2_4[5]/F");
    genTree_p->Branch("gAlgImbHem9C4_8", &gAlgImbHem9C4_8_, "gAlgImbHem9C4_8[5]/F");
    genTree_p->Branch("gAlgImbHem9C8_100", &gAlgImbHem9C8_100_, "gAlgImbHem9C8_100[5]/F");

    genTree_p->Branch("gAlgImbHem10CF", &gAlgImbHem10CF_, "gAlgImbHem10CF[5]/F");
    genTree_p->Branch("gAlgImbHem10C0_1", &gAlgImbHem10C0_1_, "gAlgImbHem10C0_1[5]/F");
    genTree_p->Branch("gAlgImbHem10C1_2", &gAlgImbHem10C1_2_, "gAlgImbHem10C1_2[5]/F");
    genTree_p->Branch("gAlgImbHem10C2_4", &gAlgImbHem10C2_4_, "gAlgImbHem10C2_4[5]/F");
    genTree_p->Branch("gAlgImbHem10C4_8", &gAlgImbHem10C4_8_, "gAlgImbHem10C4_8[5]/F");
    genTree_p->Branch("gAlgImbHem10C8_100", &gAlgImbHem10C8_100_, "gAlgImbHem10C8_100[5]/F");


    //Hem's

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
    setWeight_[initIter] = -1;
    setWeight_2pi3_[initIter] = -1;

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
    rAlgImbProjF_[initIter] = 0;
    rAlgImbPerpF_[initIter] = 0;
    rAlgImbProj0_1_[initIter] = 0;
    rAlgImbProj1_2_[initIter] = 0;
    rAlgImbProj2_4_[initIter] = 0;
    rAlgImbProj4_8_[initIter] = 0;
    rAlgImbProj8_100_[initIter] = 0;
    rAlgImbProjCF_[initIter] = 0;
    rAlgImbPerpCF_[initIter] = 0;
    rAlgImbProjC0_1_[initIter] = 0;
    rAlgImbProjC1_2_[initIter] = 0;
    rAlgImbProjC2_4_[initIter] = 0;
    rAlgImbProjC4_8_[initIter] = 0;
    rAlgImbProjC8_100_[initIter] = 0;
    rAlgImbProjNCF_[initIter] = 0;
    rAlgImbPerpNCF_[initIter] = 0;
    rAlgImbProjNC0_1_[initIter] = 0;
    rAlgImbProjNC1_2_[initIter] = 0;
    rAlgImbProjNC2_4_[initIter] = 0;
    rAlgImbProjNC4_8_[initIter] = 0;
    rAlgImbProjNC8_100_[initIter] = 0;

    rAlgImbHemF_[initIter] = 0;
    rAlgImbHem0_1_[initIter] = 0;
    rAlgImbHem1_2_[initIter] = 0;
    rAlgImbHem2_4_[initIter] = 0;
    rAlgImbHem4_8_[initIter] = 0;
    rAlgImbHem8_100_[initIter] = 0;
    rAlgImbHemCF_[initIter] = 0;
    rAlgImbHemC0_1_[initIter] = 0;
    rAlgImbHemC1_2_[initIter] = 0;
    rAlgImbHemC2_4_[initIter] = 0;
    rAlgImbHemC4_8_[initIter] = 0;
    rAlgImbHemC8_100_[initIter] = 0;
    rAlgImbHemNCF_[initIter] = 0;
    rAlgImbHemNC0_1_[initIter] = 0;
    rAlgImbHemNC1_2_[initIter] = 0;
    rAlgImbHemNC2_4_[initIter] = 0;
    rAlgImbHemNC4_8_[initIter] = 0;
    rAlgImbHemNC8_100_[initIter] = 0;

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

    //DelRs Proj

    rAlgImbProj1CF_[initIter] = 0;
    rAlgImbProj1C0_1_[initIter] = 0;
    rAlgImbProj1C1_2_[initIter] = 0;
    rAlgImbProj1C2_4_[initIter] = 0;
    rAlgImbProj1C4_8_[initIter] = 0;
    rAlgImbProj1C8_100_[initIter] = 0;

    rAlgImbProj2CF_[initIter] = 0;
    rAlgImbProj2C0_1_[initIter] = 0;
    rAlgImbProj2C1_2_[initIter] = 0;
    rAlgImbProj2C2_4_[initIter] = 0;
    rAlgImbProj2C4_8_[initIter] = 0;
    rAlgImbProj2C8_100_[initIter] = 0;

    rAlgImbProj3CF_[initIter] = 0;
    rAlgImbProj3C0_1_[initIter] = 0;
    rAlgImbProj3C1_2_[initIter] = 0;
    rAlgImbProj3C2_4_[initIter] = 0;
    rAlgImbProj3C4_8_[initIter] = 0;
    rAlgImbProj3C8_100_[initIter] = 0;

    rAlgImbProj4CF_[initIter] = 0;
    rAlgImbProj4C0_1_[initIter] = 0;
    rAlgImbProj4C1_2_[initIter] = 0;
    rAlgImbProj4C2_4_[initIter] = 0;
    rAlgImbProj4C4_8_[initIter] = 0;
    rAlgImbProj4C8_100_[initIter] = 0;

    rAlgImbProj5CF_[initIter] = 0;
    rAlgImbProj5C0_1_[initIter] = 0;
    rAlgImbProj5C1_2_[initIter] = 0;
    rAlgImbProj5C2_4_[initIter] = 0;
    rAlgImbProj5C4_8_[initIter] = 0;
    rAlgImbProj5C8_100_[initIter] = 0;

    rAlgImbProj6CF_[initIter] = 0;
    rAlgImbProj6C0_1_[initIter] = 0;
    rAlgImbProj6C1_2_[initIter] = 0;
    rAlgImbProj6C2_4_[initIter] = 0;
    rAlgImbProj6C4_8_[initIter] = 0;
    rAlgImbProj6C8_100_[initIter] = 0;

    rAlgImbProj7CF_[initIter] = 0;
    rAlgImbProj7C0_1_[initIter] = 0;
    rAlgImbProj7C1_2_[initIter] = 0;
    rAlgImbProj7C2_4_[initIter] = 0;
    rAlgImbProj7C4_8_[initIter] = 0;
    rAlgImbProj7C8_100_[initIter] = 0;

    rAlgImbProj8CF_[initIter] = 0;
    rAlgImbProj8C0_1_[initIter] = 0;
    rAlgImbProj8C1_2_[initIter] = 0;
    rAlgImbProj8C2_4_[initIter] = 0;
    rAlgImbProj8C4_8_[initIter] = 0;
    rAlgImbProj8C8_100_[initIter] = 0;

    rAlgImbProj9CF_[initIter] = 0;
    rAlgImbProj9C0_1_[initIter] = 0;
    rAlgImbProj9C1_2_[initIter] = 0;
    rAlgImbProj9C2_4_[initIter] = 0;
    rAlgImbProj9C4_8_[initIter] = 0;
    rAlgImbProj9C8_100_[initIter] = 0;

    rAlgImbProj10CF_[initIter] = 0;
    rAlgImbProj10C0_1_[initIter] = 0;
    rAlgImbProj10C1_2_[initIter] = 0;
    rAlgImbProj10C2_4_[initIter] = 0;
    rAlgImbProj10C4_8_[initIter] = 0;
    rAlgImbProj10C8_100_[initIter] = 0;

    // Hem

    rAlgImbHem1CF_[initIter] = 0;
    rAlgImbHem1C0_1_[initIter] = 0;
    rAlgImbHem1C1_2_[initIter] = 0;
    rAlgImbHem1C2_4_[initIter] = 0;
    rAlgImbHem1C4_8_[initIter] = 0;
    rAlgImbHem1C8_100_[initIter] = 0;

    rAlgImbHem2CF_[initIter] = 0;
    rAlgImbHem2C0_1_[initIter] = 0;
    rAlgImbHem2C1_2_[initIter] = 0;
    rAlgImbHem2C2_4_[initIter] = 0;
    rAlgImbHem2C4_8_[initIter] = 0;
    rAlgImbHem2C8_100_[initIter] = 0;

    rAlgImbHem3CF_[initIter] = 0;
    rAlgImbHem3C0_1_[initIter] = 0;
    rAlgImbHem3C1_2_[initIter] = 0;
    rAlgImbHem3C2_4_[initIter] = 0;
    rAlgImbHem3C4_8_[initIter] = 0;
    rAlgImbHem3C8_100_[initIter] = 0;

    rAlgImbHem4CF_[initIter] = 0;
    rAlgImbHem4C0_1_[initIter] = 0;
    rAlgImbHem4C1_2_[initIter] = 0;
    rAlgImbHem4C2_4_[initIter] = 0;
    rAlgImbHem4C4_8_[initIter] = 0;
    rAlgImbHem4C8_100_[initIter] = 0;

    rAlgImbHem5CF_[initIter] = 0;
    rAlgImbHem5C0_1_[initIter] = 0;
    rAlgImbHem5C1_2_[initIter] = 0;
    rAlgImbHem5C2_4_[initIter] = 0;
    rAlgImbHem5C4_8_[initIter] = 0;
    rAlgImbHem5C8_100_[initIter] = 0;

    rAlgImbHem6CF_[initIter] = 0;
    rAlgImbHem6C0_1_[initIter] = 0;
    rAlgImbHem6C1_2_[initIter] = 0;
    rAlgImbHem6C2_4_[initIter] = 0;
    rAlgImbHem6C4_8_[initIter] = 0;
    rAlgImbHem6C8_100_[initIter] = 0;

    rAlgImbHem7CF_[initIter] = 0;
    rAlgImbHem7C0_1_[initIter] = 0;
    rAlgImbHem7C1_2_[initIter] = 0;
    rAlgImbHem7C2_4_[initIter] = 0;
    rAlgImbHem7C4_8_[initIter] = 0;
    rAlgImbHem7C8_100_[initIter] = 0;

    rAlgImbHem8CF_[initIter] = 0;
    rAlgImbHem8C0_1_[initIter] = 0;
    rAlgImbHem8C1_2_[initIter] = 0;
    rAlgImbHem8C2_4_[initIter] = 0;
    rAlgImbHem8C4_8_[initIter] = 0;
    rAlgImbHem8C8_100_[initIter] = 0;

    rAlgImbHem9CF_[initIter] = 0;
    rAlgImbHem9C0_1_[initIter] = 0;
    rAlgImbHem9C1_2_[initIter] = 0;
    rAlgImbHem9C2_4_[initIter] = 0;
    rAlgImbHem9C4_8_[initIter] = 0;
    rAlgImbHem9C8_100_[initIter] = 0;

    rAlgImbHem10CF_[initIter] = 0;
    rAlgImbHem10C0_1_[initIter] = 0;
    rAlgImbHem10C1_2_[initIter] = 0;
    rAlgImbHem10C2_4_[initIter] = 0;
    rAlgImbHem10C4_8_[initIter] = 0;
    rAlgImbHem10C8_100_[initIter] = 0;

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
  }

  if(montecarlo){
    //Gen. proj. onto Truth
    for(Int_t initIter = 0; initIter < 5; initIter++){
      gAlgImbProjF_[initIter] = 0;
      gAlgImbPerpF_[initIter] = 0;
      gAlgImbProj0_1_[initIter] = 0;
      gAlgImbProj1_2_[initIter] = 0;
      gAlgImbProj2_4_[initIter] = 0;
      gAlgImbProj4_8_[initIter] = 0;
      gAlgImbProj8_100_[initIter] = 0;
      gAlgImbProjCF_[initIter] = 0;
      gAlgImbPerpCF_[initIter] = 0;
      gAlgImbProjC0_1_[initIter] = 0;
      gAlgImbProjC1_2_[initIter] = 0;
      gAlgImbProjC2_4_[initIter] = 0;
      gAlgImbProjC4_8_[initIter] = 0;
      gAlgImbProjC8_100_[initIter] = 0;
      gAlgImbProjNCF_[initIter] = 0;
      gAlgImbPerpNCF_[initIter] = 0;
      gAlgImbProjNC0_1_[initIter] = 0;
      gAlgImbProjNC1_2_[initIter] = 0;
      gAlgImbProjNC2_4_[initIter] = 0;
      gAlgImbProjNC4_8_[initIter] = 0;
      gAlgImbProjNC8_100_[initIter] = 0;

      gAlgImbProjSigCF_[initIter] = 0;
      gAlgImbPerpSigCF_[initIter] = 0;
      gAlgImbProjSigC0_1_[initIter] = 0;
      gAlgImbProjSigC1_2_[initIter] = 0;
      gAlgImbProjSigC2_4_[initIter] = 0;
      gAlgImbProjSigC4_8_[initIter] = 0;
      gAlgImbProjSigC8_100_[initIter] = 0;
      gAlgImbProjSigNCF_[initIter] = 0;
      gAlgImbPerpSigNCF_[initIter] = 0;
      gAlgImbProjSigNC0_1_[initIter] = 0;
      gAlgImbProjSigNC1_2_[initIter] = 0;
      gAlgImbProjSigNC2_4_[initIter] = 0;
      gAlgImbProjSigNC4_8_[initIter] = 0;
      gAlgImbProjSigNC8_100_[initIter] = 0;

      gAlgImbProjUECF_[initIter] = 0;
      gAlgImbPerpUECF_[initIter] = 0;
      gAlgImbProjUEC0_1_[initIter] = 0;
      gAlgImbProjUEC1_2_[initIter] = 0;
      gAlgImbProjUEC2_4_[initIter] = 0;
      gAlgImbProjUEC4_8_[initIter] = 0;
      gAlgImbProjUEC8_100_[initIter] = 0;
      gAlgImbProjUENCF_[initIter] = 0;
      gAlgImbPerpUENCF_[initIter] = 0;
      gAlgImbProjUENC0_1_[initIter] = 0;
      gAlgImbProjUENC1_2_[initIter] = 0;
      gAlgImbProjUENC2_4_[initIter] = 0;
      gAlgImbProjUENC4_8_[initIter] = 0;
      gAlgImbProjUENC8_100_[initIter] = 0;

      gAlgImbHemF_[initIter] = 0;
      gAlgImbHem0_1_[initIter] = 0;
      gAlgImbHem1_2_[initIter] = 0;
      gAlgImbHem2_4_[initIter] = 0;
      gAlgImbHem4_8_[initIter] = 0;
      gAlgImbHem8_100_[initIter] = 0;
      gAlgImbHemCF_[initIter] = 0;
      gAlgImbHemC0_1_[initIter] = 0;
      gAlgImbHemC1_2_[initIter] = 0;
      gAlgImbHemC2_4_[initIter] = 0;
      gAlgImbHemC4_8_[initIter] = 0;
      gAlgImbHemC8_100_[initIter] = 0;
      gAlgImbHemNCF_[initIter] = 0;
      gAlgImbHemNC0_1_[initIter] = 0;
      gAlgImbHemNC1_2_[initIter] = 0;
      gAlgImbHemNC2_4_[initIter] = 0;
      gAlgImbHemNC4_8_[initIter] = 0;
      gAlgImbHemNC8_100_[initIter] = 0;

      gAlgImbHemSigCF_[initIter] = 0;
      gAlgImbHemSigC0_1_[initIter] = 0;
      gAlgImbHemSigC1_2_[initIter] = 0;
      gAlgImbHemSigC2_4_[initIter] = 0;
      gAlgImbHemSigC4_8_[initIter] = 0;
      gAlgImbHemSigC8_100_[initIter] = 0;
      gAlgImbHemSigNCF_[initIter] = 0;
      gAlgImbHemSigNC0_1_[initIter] = 0;
      gAlgImbHemSigNC1_2_[initIter] = 0;
      gAlgImbHemSigNC2_4_[initIter] = 0;
      gAlgImbHemSigNC4_8_[initIter] = 0;
      gAlgImbHemSigNC8_100_[initIter] = 0;

      gAlgImbHemUECF_[initIter] = 0;
      gAlgImbHemUEC0_1_[initIter] = 0;
      gAlgImbHemUEC1_2_[initIter] = 0;
      gAlgImbHemUEC2_4_[initIter] = 0;
      gAlgImbHemUEC4_8_[initIter] = 0;
      gAlgImbHemUEC8_100_[initIter] = 0;
      gAlgImbHemUENCF_[initIter] = 0;
      gAlgImbHemUENC0_1_[initIter] = 0;
      gAlgImbHemUENC1_2_[initIter] = 0;
      gAlgImbHemUENC2_4_[initIter] = 0;
      gAlgImbHemUENC4_8_[initIter] = 0;
      gAlgImbHemUENC8_100_[initIter] = 0;

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

      gAlgImbProj1CF_[initIter] = 0;
      gAlgImbProj1C0_1_[initIter] = 0;
      gAlgImbProj1C1_2_[initIter] = 0;
      gAlgImbProj1C2_4_[initIter] = 0;
      gAlgImbProj1C4_8_[initIter] = 0;
      gAlgImbProj1C8_100_[initIter] = 0;

      gAlgImbProj2CF_[initIter] = 0;
      gAlgImbProj2C0_1_[initIter] = 0;
      gAlgImbProj2C1_2_[initIter] = 0;
      gAlgImbProj2C2_4_[initIter] = 0;
      gAlgImbProj2C4_8_[initIter] = 0;
      gAlgImbProj2C8_100_[initIter] = 0;

      gAlgImbProj3CF_[initIter] = 0;
      gAlgImbProj3C0_1_[initIter] = 0;
      gAlgImbProj3C1_2_[initIter] = 0;
      gAlgImbProj3C2_4_[initIter] = 0;
      gAlgImbProj3C4_8_[initIter] = 0;
      gAlgImbProj3C8_100_[initIter] = 0;

      gAlgImbProj4CF_[initIter] = 0;
      gAlgImbProj4C0_1_[initIter] = 0;
      gAlgImbProj4C1_2_[initIter] = 0;
      gAlgImbProj4C2_4_[initIter] = 0;
      gAlgImbProj4C4_8_[initIter] = 0;
      gAlgImbProj4C8_100_[initIter] = 0;

      gAlgImbProj5CF_[initIter] = 0;
      gAlgImbProj5C0_1_[initIter] = 0;
      gAlgImbProj5C1_2_[initIter] = 0;
      gAlgImbProj5C2_4_[initIter] = 0;
      gAlgImbProj5C4_8_[initIter] = 0;
      gAlgImbProj5C8_100_[initIter] = 0;

      gAlgImbProj6CF_[initIter] = 0;
      gAlgImbProj6C0_1_[initIter] = 0;
      gAlgImbProj6C1_2_[initIter] = 0;
      gAlgImbProj6C2_4_[initIter] = 0;
      gAlgImbProj6C4_8_[initIter] = 0;
      gAlgImbProj6C8_100_[initIter] = 0;

      gAlgImbProj7CF_[initIter] = 0;
      gAlgImbProj7C0_1_[initIter] = 0;
      gAlgImbProj7C1_2_[initIter] = 0;
      gAlgImbProj7C2_4_[initIter] = 0;
      gAlgImbProj7C4_8_[initIter] = 0;
      gAlgImbProj7C8_100_[initIter] = 0;

      gAlgImbProj8CF_[initIter] = 0;
      gAlgImbProj8C0_1_[initIter] = 0;
      gAlgImbProj8C1_2_[initIter] = 0;
      gAlgImbProj8C2_4_[initIter] = 0;
      gAlgImbProj8C4_8_[initIter] = 0;
      gAlgImbProj8C8_100_[initIter] = 0;

      gAlgImbProj9CF_[initIter] = 0;
      gAlgImbProj9C0_1_[initIter] = 0;
      gAlgImbProj9C1_2_[initIter] = 0;
      gAlgImbProj9C2_4_[initIter] = 0;
      gAlgImbProj9C4_8_[initIter] = 0;
      gAlgImbProj9C8_100_[initIter] = 0;

      gAlgImbProj10CF_[initIter] = 0;
      gAlgImbProj10C0_1_[initIter] = 0;
      gAlgImbProj10C1_2_[initIter] = 0;
      gAlgImbProj10C2_4_[initIter] = 0;
      gAlgImbProj10C4_8_[initIter] = 0;
      gAlgImbProj10C8_100_[initIter] = 0;

      // hems

      gAlgImbHem1CF_[initIter] = 0;
      gAlgImbHem1C0_1_[initIter] = 0;
      gAlgImbHem1C1_2_[initIter] = 0;
      gAlgImbHem1C2_4_[initIter] = 0;
      gAlgImbHem1C4_8_[initIter] = 0;
      gAlgImbHem1C8_100_[initIter] = 0;

      gAlgImbHem2CF_[initIter] = 0;
      gAlgImbHem2C0_1_[initIter] = 0;
      gAlgImbHem2C1_2_[initIter] = 0;
      gAlgImbHem2C2_4_[initIter] = 0;
      gAlgImbHem2C4_8_[initIter] = 0;
      gAlgImbHem2C8_100_[initIter] = 0;

      gAlgImbHem3CF_[initIter] = 0;
      gAlgImbHem3C0_1_[initIter] = 0;
      gAlgImbHem3C1_2_[initIter] = 0;
      gAlgImbHem3C2_4_[initIter] = 0;
      gAlgImbHem3C4_8_[initIter] = 0;
      gAlgImbHem3C8_100_[initIter] = 0;

      gAlgImbHem4CF_[initIter] = 0;
      gAlgImbHem4C0_1_[initIter] = 0;
      gAlgImbHem4C1_2_[initIter] = 0;
      gAlgImbHem4C2_4_[initIter] = 0;
      gAlgImbHem4C4_8_[initIter] = 0;
      gAlgImbHem4C8_100_[initIter] = 0;

      gAlgImbHem5CF_[initIter] = 0;
      gAlgImbHem5C0_1_[initIter] = 0;
      gAlgImbHem5C1_2_[initIter] = 0;
      gAlgImbHem5C2_4_[initIter] = 0;
      gAlgImbHem5C4_8_[initIter] = 0;
      gAlgImbHem5C8_100_[initIter] = 0;

      gAlgImbHem6CF_[initIter] = 0;
      gAlgImbHem6C0_1_[initIter] = 0;
      gAlgImbHem6C1_2_[initIter] = 0;
      gAlgImbHem6C2_4_[initIter] = 0;
      gAlgImbHem6C4_8_[initIter] = 0;
      gAlgImbHem6C8_100_[initIter] = 0;

      gAlgImbHem7CF_[initIter] = 0;
      gAlgImbHem7C0_1_[initIter] = 0;
      gAlgImbHem7C1_2_[initIter] = 0;
      gAlgImbHem7C2_4_[initIter] = 0;
      gAlgImbHem7C4_8_[initIter] = 0;
      gAlgImbHem7C8_100_[initIter] = 0;

      gAlgImbHem8CF_[initIter] = 0;
      gAlgImbHem8C0_1_[initIter] = 0;
      gAlgImbHem8C1_2_[initIter] = 0;
      gAlgImbHem8C2_4_[initIter] = 0;
      gAlgImbHem8C4_8_[initIter] = 0;
      gAlgImbHem8C8_100_[initIter] = 0;

      gAlgImbHem9CF_[initIter] = 0;
      gAlgImbHem9C0_1_[initIter] = 0;
      gAlgImbHem9C1_2_[initIter] = 0;
      gAlgImbHem9C2_4_[initIter] = 0;
      gAlgImbHem9C4_8_[initIter] = 0;
      gAlgImbHem9C8_100_[initIter] = 0;

      gAlgImbHem10CF_[initIter] = 0;
      gAlgImbHem10C0_1_[initIter] = 0;
      gAlgImbHem10C1_2_[initIter] = 0;
      gAlgImbHem10C2_4_[initIter] = 0;
      gAlgImbHem10C4_8_[initIter] = 0;
      gAlgImbHem10C8_100_[initIter] = 0;

      //proj as

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

    }    
  }

}

#endif
