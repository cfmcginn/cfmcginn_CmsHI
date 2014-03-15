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
Float_t trkTLeadDelPhi_[MAXTRKS];
Float_t trkVsPFLeadDelPhi_[MAXTRKS];
Float_t trkVsCaloLeadDelPhi_[MAXTRKS];

Float_t trkRLeadPF_[MAXTRKS];
Float_t trkRSubLeadPF_[MAXTRKS];
Float_t trkRMinPF_[MAXTRKS];
Float_t trkPtCorrPF_[MAXTRKS];
Float_t trkPtFactPF_[MAXTRKS];

Float_t trkRLeadCalo_[MAXTRKS];
Float_t trkRSubLeadCalo_[MAXTRKS];
Float_t trkRMinCalo_[MAXTRKS];
Float_t trkPtCorrCalo_[MAXTRKS];
Float_t trkPtFactCalo_[MAXTRKS];

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

//Jet Set, Array by algorithm, according to enum above

Float_t AlgLeadJtPt_[5];
Float_t AlgLeadJtPhi_[5];
Float_t AlgLeadJtEta_[5];
Float_t AlgSubLeadJtPt_[5];
Float_t AlgSubLeadJtPhi_[5];
Float_t AlgSubLeadJtEta_[5];
Float_t AlgJtDelPhi_[5];
Float_t AlgJtAsymm_[5];
Float_t AlgLeadRefPt_[5];
Float_t AlgLeadRefEta_[5];
Float_t AlgSubLeadRefPt_[5];
Float_t AlgSubLeadRefEta_[5];

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
  trackTree_p->Branch("trkTLeadDelPhi", &trkTLeadDelPhi_, "trkTLeadDelPhi[nTrk]/F");
  trackTree_p->Branch("trkVsPFLeadDelPhi", &trkVsPFLeadDelPhi_, "trkVsPFLeadDelPhi[nTrk]/F");
  trackTree_p->Branch("trkVsCaloLeadDelPhi", &trkVsCaloLeadDelPhi_, "trkVsCaloLeadDelPhi[nTrk]/F");

  trackTree_p->Branch("trkRLeadPF", &trkRLeadPF_, "trkRLeadPF[nTrk]/F");
  trackTree_p->Branch("trkRSubLeadPF", &trkRSubLeadPF_, "trkRSubLeadPF[nTrk]/F");
  trackTree_p->Branch("trkRMinPF", &trkRMinPF_, "trkRMinPF[nTrk]/F");
  trackTree_p->Branch("trkPtCorrPF", &trkPtCorrPF_, "trkPtCorrPF[nTrk]/F");
  trackTree_p->Branch("trkPtFactPF", &trkPtFactPF_, "trkPtFactPF[nTrk]/F");

  trackTree_p->Branch("trkRLeadCalo", &trkRLeadCalo_, "trkRLeadCalo[nTrk]/F");
  trackTree_p->Branch("trkRSubLeadCalo", &trkRSubLeadCalo_, "trkRSubLeadCalo[nTrk]/F");
  trackTree_p->Branch("trkRMinCalo", &trkRMinCalo_, "trkRMinCalo[nTrk]/F");
  trackTree_p->Branch("trkPtCorrCalo", &trkPtCorrCalo_, "trkPtCorrCalo[nTrk]/F");
  trackTree_p->Branch("trkPtFactCalo", &trkPtFactCalo_, "trkPtFactCalo[nTrk]/F");

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

  jetTree_p->Branch("AlgLeadJtPt", &AlgLeadJtPt_, "AlgLeadJtPt[5]/F");
  jetTree_p->Branch("AlgLeadJtPhi", &AlgLeadJtPhi_, "AlgLeadJtPhi[5]/F");
  jetTree_p->Branch("AlgLeadJtEta", &AlgLeadJtEta_, "AlgLeadJtEta[5]/F");
  jetTree_p->Branch("AlgSubLeadJtPt", &AlgSubLeadJtPt_, "AlgSubLeadJtPt[5]/F");
  jetTree_p->Branch("AlgSubLeadJtPhi", &AlgSubLeadJtPhi_, "AlgSubLeadJtPhi[5]/F");
  jetTree_p->Branch("AlgSubLeadJtEta", &AlgSubLeadJtEta_, "AlgSubLeadJtEta[5]/F");
  jetTree_p->Branch("AlgJtDelPhi", &AlgJtDelPhi_, "AlgJtDelPhi[5]/F");
  jetTree_p->Branch("AlgJtAsymm", &AlgJtAsymm_, "AlgJtAsymm[5]/F");

  if(montecarlo){
    //refpt for jets immediately above
    jetTree_p->Branch("AlgLeadRefPt", &AlgLeadRefPt_, "AlgLeadRefPt[5]/F");
    jetTree_p->Branch("AlgLeadRefEta", &AlgLeadRefEta_, "AlgLeadRefEta[5]/F");
    jetTree_p->Branch("AlgSubLeadRefPt", &AlgSubLeadRefPt_, "AlgSubLeadRefPt[5]/F");
    jetTree_p->Branch("AlgSubLeadRefEta", &AlgSubLeadRefEta_, "AlgSubLeadRefEta[5]/F");

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
  trackTree_p->SetBranchAddress("trkTLeadDelPhi", &trkTLeadDelPhi_);
  trackTree_p->SetBranchAddress("trkVsPFLeadDelPhi", &trkVsPFLeadDelPhi_);
  trackTree_p->SetBranchAddress("trkVsCaloLeadDelPhi", &trkVsCaloLeadDelPhi_);


  trackTree_p->SetBranchAddress("trkRLeadPF", &trkRLeadPF_);
  trackTree_p->SetBranchAddress("trkRSubLeadPF", &trkRSubLeadPF_);
  trackTree_p->SetBranchAddress("trkRMinPF", &trkRMinPF_);
  trackTree_p->SetBranchAddress("trkPtCorrPF", &trkPtCorrPF_);
  trackTree_p->SetBranchAddress("trkPtFactPF", &trkPtFactPF_);

  trackTree_p->SetBranchAddress("trkRLeadCalo", &trkRLeadCalo_);
  trackTree_p->SetBranchAddress("trkRSubLeadCalo", &trkRSubLeadCalo_);
  trackTree_p->SetBranchAddress("trkRMinCalo", &trkRMinCalo_);
  trackTree_p->SetBranchAddress("trkPtCorrCalo", &trkPtCorrCalo_);
  trackTree_p->SetBranchAddress("trkPtFactCalo", &trkPtFactCalo_);

  if(montecarlo){
    trackTree_p->SetBranchAddress("trkRLeadT", &trkRLeadT_);
    trackTree_p->SetBranchAddress("trkRSubLeadT", &trkRSubLeadT_);
    trackTree_p->SetBranchAddress("trkRMinT", &trkRMinT_);
    trackTree_p->SetBranchAddress("trkPtCorrT", &trkPtCorrT_);
    trackTree_p->SetBranchAddress("trkPtFactT", &trkPtFactT_);
  }

  trackTree_p->SetBranchAddress("trkRLeadVsPF", &trkRLeadVsPF_);
  trackTree_p->SetBranchAddress("trkRSubLeadVsPF", &trkRSubLeadVsPF_);
  trackTree_p->SetBranchAddress("trkRMinVsPF", &trkRMinVsPF_);
  trackTree_p->SetBranchAddress("trkPtCorrVsPF", &trkPtCorrVsPF_);
  trackTree_p->SetBranchAddress("trkPtFactVsPF", &trkPtFactVsPF_);

  trackTree_p->SetBranchAddress("trkRLeadVsCalo", &trkRLeadVsCalo_);
  trackTree_p->SetBranchAddress("trkRSubLeadVsCalo", &trkRSubLeadVsCalo_);
  trackTree_p->SetBranchAddress("trkRMinVsCalo", &trkRMinVsCalo_);
  trackTree_p->SetBranchAddress("trkPtCorrVsCalo", &trkPtCorrVsCalo_);
  trackTree_p->SetBranchAddress("trkPtFactVsCalo", &trkPtFactVsCalo_);

  //Tracks proj. onto Alg, ordered according to enum above (corr in back 5), All, Cone, and NotCone

  trackTree_p->SetBranchAddress("rAlgImbProjF", &rAlgImbProjF_);
  trackTree_p->SetBranchAddress("rAlgImbPerpF", &rAlgImbPerpF_);
  trackTree_p->SetBranchAddress("rAlgImbProj0_1", &rAlgImbProj0_1_);
  trackTree_p->SetBranchAddress("rAlgImbProj1_2", &rAlgImbProj1_2_);
  trackTree_p->SetBranchAddress("rAlgImbProj2_4", &rAlgImbProj2_4_);
  trackTree_p->SetBranchAddress("rAlgImbProj4_8", &rAlgImbProj4_8_);
  trackTree_p->SetBranchAddress("rAlgImbProj8_100", &rAlgImbProj8_100_);
  trackTree_p->SetBranchAddress("rAlgImbProjCF", &rAlgImbProjCF_);
  trackTree_p->SetBranchAddress("rAlgImbPerpCF", &rAlgImbPerpCF_);
  trackTree_p->SetBranchAddress("rAlgImbProjC0_1", &rAlgImbProjC0_1_);
  trackTree_p->SetBranchAddress("rAlgImbProjC1_2", &rAlgImbProjC1_2_);
  trackTree_p->SetBranchAddress("rAlgImbProjC2_4", &rAlgImbProjC2_4_);
  trackTree_p->SetBranchAddress("rAlgImbProjC4_8", &rAlgImbProjC4_8_);
  trackTree_p->SetBranchAddress("rAlgImbProjC8_100", &rAlgImbProjC8_100_);
  trackTree_p->SetBranchAddress("rAlgImbProjNCF", &rAlgImbProjNCF_);
  trackTree_p->SetBranchAddress("rAlgImbPerpNCF", &rAlgImbPerpNCF_);
  trackTree_p->SetBranchAddress("rAlgImbProjNC0_1", &rAlgImbProjNC0_1_);
  trackTree_p->SetBranchAddress("rAlgImbProjNC1_2", &rAlgImbProjNC1_2_);
  trackTree_p->SetBranchAddress("rAlgImbProjNC2_4", &rAlgImbProjNC2_4_);
  trackTree_p->SetBranchAddress("rAlgImbProjNC4_8", &rAlgImbProjNC4_8_);
  trackTree_p->SetBranchAddress("rAlgImbProjNC8_100", &rAlgImbProjNC8_100_);


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

  jetTree_p->SetBranchAddress("eventSet", &eventSet_);

  jetTree_p->SetBranchAddress("AlgLeadJtPt", &AlgLeadJtPt_);
  jetTree_p->SetBranchAddress("AlgLeadJtPhi", &AlgLeadJtPhi_);
  jetTree_p->SetBranchAddress("AlgLeadJtEta", &AlgLeadJtEta_);
  jetTree_p->SetBranchAddress("AlgSubLeadJtPt", &AlgSubLeadJtPt_);
  jetTree_p->SetBranchAddress("AlgSubLeadJtPhi", &AlgSubLeadJtPhi_);
  jetTree_p->SetBranchAddress("AlgSubLeadJtEta", &AlgSubLeadJtEta_);
  jetTree_p->SetBranchAddress("AlgJtDelPhi", &AlgJtDelPhi_);
  jetTree_p->SetBranchAddress("AlgJtAsymm", &AlgJtAsymm_);

  if(montecarlo){
    //Ref pt for jet var immediately above
    jetTree_p->SetBranchAddress("AlgLeadRefPt", &AlgLeadRefPt_);
    jetTree_p->SetBranchAddress("AlgLeadRefEta", &AlgLeadRefEta_);
    jetTree_p->SetBranchAddress("AlgSubLeadRefPt", &AlgSubLeadRefPt_);
    jetTree_p->SetBranchAddress("AlgSubLeadRefEta", &AlgSubLeadRefEta_);

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

    genTree_p->SetBranchAddress("gAlgImbProjF", &gAlgImbProjF_);
    genTree_p->SetBranchAddress("gAlgImbPerpF", &gAlgImbPerpF_);
    genTree_p->SetBranchAddress("gAlgImbProj0_1", &gAlgImbProj0_1_);
    genTree_p->SetBranchAddress("gAlgImbProj1_2", &gAlgImbProj1_2_);
    genTree_p->SetBranchAddress("gAlgImbProj2_4", &gAlgImbProj2_4_);
    genTree_p->SetBranchAddress("gAlgImbProj4_8", &gAlgImbProj4_8_);
    genTree_p->SetBranchAddress("gAlgImbProj8_100", &gAlgImbProj8_100_);

    genTree_p->SetBranchAddress("gAlgImbProjCF", &gAlgImbProjCF_);
    genTree_p->SetBranchAddress("gAlgImbPerpCF", &gAlgImbPerpCF_);
    genTree_p->SetBranchAddress("gAlgImbProjC0_1", &gAlgImbProjC0_1_);
    genTree_p->SetBranchAddress("gAlgImbProjC1_2", &gAlgImbProjC1_2_);
    genTree_p->SetBranchAddress("gAlgImbProjC2_4", &gAlgImbProjC2_4_);
    genTree_p->SetBranchAddress("gAlgImbProjC4_8", &gAlgImbProjC4_8_);
    genTree_p->SetBranchAddress("gAlgImbProjC8_100", &gAlgImbProjC8_100_);

    genTree_p->SetBranchAddress("gAlgImbProjNCF", &gAlgImbProjNCF_);
    genTree_p->SetBranchAddress("gAlgImbPerpNCF", &gAlgImbPerpNCF_);
    genTree_p->SetBranchAddress("gAlgImbProjNC0_1", &gAlgImbProjNC0_1_);
    genTree_p->SetBranchAddress("gAlgImbProjNC1_2", &gAlgImbProjNC1_2_);
    genTree_p->SetBranchAddress("gAlgImbProjNC2_4", &gAlgImbProjNC2_4_);
    genTree_p->SetBranchAddress("gAlgImbProjNC4_8", &gAlgImbProjNC4_8_);
    genTree_p->SetBranchAddress("gAlgImbProjNC8_100", &gAlgImbProjNC8_100_);
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
  for(Int_t initIter = 0; initIter < 5; initIter++){
    eventSet_[initIter] = false;

    AlgLeadJtPt_[initIter] = -10;
    AlgSubLeadJtPt_[initIter] = -10;
    AlgLeadJtPhi_[initIter] = -10;
    AlgSubLeadJtPhi_[initIter] = -10;
    AlgLeadJtEta_[initIter] = -10;
    AlgSubLeadJtEta_[initIter] = -10;
    AlgJtDelPhi_[initIter] = -10;
    AlgJtAsymm_[initIter] = -10;

    if(montecarlo){
      AlgLeadRefPt_[initIter] = -10;
      AlgSubLeadRefPt_[initIter] = -10;
      AlgLeadRefEta_[initIter] = -10;
      AlgSubLeadRefEta_[initIter] = -10;
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
    }    
  }

}

#endif
