//=============================================
// Author: Chris McGinn
// 
// DiJet Analysis Class (MC)
//
//=============================================

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
//#include "/net/hisrv0001/home/cfmcginn/emDiJet/CMSSW_5_3_12_patch3/HiForestAnalysis/hiForest.h"
//#include "../../HiForestAnalysis/hiForest.h"
//Version of hiForest w/ correct track array size, and Voronoi Subtracted jets
#include "/net/hisrv0001/home/cfmcginn/emDiJet/CMSSW_5_3_12_patch3/tempHIFA/HiForestAnalysis/hiForest.h"
//#include "../../tempHIFA/HiForestAnalysis/hiForest.h"
#include "commonUtility.h"
#include "cfmDiJetSkim.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include "factorizedPtCorr.h"

const Float_t leadJtPtCut = 120.;
const Float_t subLeadJtPtCut = 50.;
const Float_t jtDelPhiCut = 7.*(TMath::Pi())/8.;
const Float_t jtEtaCut = 1.6; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed

collisionType getCType(sampleType sType);

//PF corr profs

const char* PFcorrName0_1 = "eff_pt0_1_step_cent3accept3pt3rmin2.root";
const char* PFcorrName1_3 = "eff_pt1_3_step_cent3accept3pt3rmin2.root";
const char* PFcorrName3_8 = "eff_pt3_8_step_cent3accept3pt3rmin2.root";
const char* PFcorrName8_100 = "eff_pt8_100_step_cent3accept3pt3rmin2.root";

TProfile* PFcent0_1_p;
TProfile* PFcent1_3_p;
TProfile* PFcent3_8_p;
TProfile* PFcent8_100_p;

TProfile2D* PFphiEta0_1_p;
TProfile2D* PFphiEta1_3_p;
TProfile2D* PFphiEta3_8_p;
TProfile2D* PFphiEta8_100_p;

TProfile* PFpt0_1_p;
TProfile* PFpt1_3_p;
TProfile* PFpt3_8_p;
TProfile* PFpt8_100_p;

TProfile* PFdelR0_1_p;
TProfile* PFdelR1_3_p;
TProfile* PFdelR3_8_p;
TProfile* PFdelR8_100_p;

//Calo corr profs

TProfile* Calocent0_1_p;
TProfile* Calocent1_3_p;
TProfile* Calocent3_8_p;
TProfile* Calocent8_100_p;

TProfile2D* CalophiEta0_1_p;
TProfile2D* CalophiEta1_3_p;
TProfile2D* CalophiEta3_8_p;
TProfile2D* CalophiEta8_100_p;

TProfile* Calopt0_1_p;
TProfile* Calopt1_3_p;
TProfile* Calopt3_8_p;
TProfile* Calopt8_100_p;

TProfile* CalodelR0_1_p;
TProfile* CalodelR1_3_p;
TProfile* CalodelR3_8_p;
TProfile* CalodelR8_100_p;

//T corr profs

const char* TcorrName0_1 = "eff_pt0_1_step_cent3accept3pt3rmin2_akPu3PF_dogenjet1.root";
const char* TcorrName1_3 = "eff_pt1_3_step_cent3accept3pt3rmin2_akPu3PF_dogenjet1.root";
const char* TcorrName3_8 = "eff_pt3_8_step_cent3accept3pt3rmin2_akPu3PF_dogenjet1.root";
const char* TcorrName8_100 = "eff_pt8_100_step_cent3accept3pt3rmin2_akPu3PF_dogenjet1.root";

TProfile* Tcent0_1_p;
TProfile* Tcent1_3_p;
TProfile* Tcent3_8_p;
TProfile* Tcent8_100_p;

TProfile2D* TphiEta0_1_p;
TProfile2D* TphiEta1_3_p;
TProfile2D* TphiEta3_8_p;
TProfile2D* TphiEta8_100_p;

TProfile* Tpt0_1_p;
TProfile* Tpt1_3_p;
TProfile* Tpt3_8_p;
TProfile* Tpt8_100_p;

TProfile* TdelR0_1_p;
TProfile* TdelR1_3_p;
TProfile* TdelR3_8_p;
TProfile* TdelR8_100_p;

int makeDiJetTTree(string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_CFMSKIM.root")
{
  //Define MC or Data
  bool montecarlo = false;
  if(sType == kPPMC || sType == kPAMC || sType == kHIMC)
    montecarlo = true;

  std::cout << montecarlo << std::endl;

  collisionType cType = getCType(sType);

  string buffer;
  std::vector<string> listOfFiles;
  int nLines = 0;
  ifstream inFile(fList.data());

  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
    return 1;
  }
  else{
    while(!inFile.eof()){
      inFile >> buffer;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  std::cout << "FileList Loaded" << std::endl;

  //Get PF Profs

  TFile* PFcorrFile0_1_p = new TFile(PFcorrName0_1, "READ");
  PFcent0_1_p = (TProfile*)PFcorrFile0_1_p->Get("p_eff_cent");
  PFphiEta0_1_p = (TProfile2D*)PFcorrFile0_1_p->Get("p_eff_acceptance");
  PFpt0_1_p = (TProfile*)PFcorrFile0_1_p->Get("p_eff_pt");
  PFdelR0_1_p = (TProfile*)PFcorrFile0_1_p->Get("p_eff_rmin");

  TFile* PFcorrFile1_3_p = new TFile(PFcorrName1_3, "READ");
  PFcent1_3_p = (TProfile*)PFcorrFile1_3_p->Get("p_eff_cent");
  PFphiEta1_3_p = (TProfile2D*)PFcorrFile1_3_p->Get("p_eff_acceptance");
  PFpt1_3_p = (TProfile*)PFcorrFile1_3_p->Get("p_eff_pt");
  PFdelR1_3_p = (TProfile*)PFcorrFile1_3_p->Get("p_eff_rmin");

  TFile* PFcorrFile3_8_p = new TFile(PFcorrName3_8, "READ");
  PFcent3_8_p = (TProfile*)PFcorrFile3_8_p->Get("p_eff_cent");
  PFphiEta3_8_p = (TProfile2D*)PFcorrFile3_8_p->Get("p_eff_acceptance");
  PFpt3_8_p = (TProfile*)PFcorrFile3_8_p->Get("p_eff_pt");
  PFdelR3_8_p = (TProfile*)PFcorrFile3_8_p->Get("p_eff_rmin");

  TFile* PFcorrFile8_100_p = new TFile(PFcorrName8_100, "READ");
  PFcent8_100_p = (TProfile*)PFcorrFile8_100_p->Get("p_eff_cent");
  PFphiEta8_100_p = (TProfile2D*)PFcorrFile8_100_p->Get("p_eff_acceptance");
  PFpt8_100_p = (TProfile*)PFcorrFile8_100_p->Get("p_eff_pt");
  PFdelR8_100_p = (TProfile*)PFcorrFile8_100_p->Get("p_eff_rmin");



  //Get T Profs
  TFile* TcorrFile0_1_p = new TFile(TcorrName0_1, "READ");
  Tcent0_1_p = (TProfile*)TcorrFile0_1_p->Get("p_eff_cent");
  TphiEta0_1_p = (TProfile2D*)TcorrFile0_1_p->Get("p_eff_acceptance");
  Tpt0_1_p = (TProfile*)TcorrFile0_1_p->Get("p_eff_pt");
  TdelR0_1_p = (TProfile*)TcorrFile0_1_p->Get("p_eff_rmin");

  TFile* TcorrFile1_3_p = new TFile(TcorrName1_3, "READ");
  Tcent1_3_p = (TProfile*)TcorrFile1_3_p->Get("p_eff_cent");
  TphiEta1_3_p = (TProfile2D*)TcorrFile1_3_p->Get("p_eff_acceptance");
  Tpt1_3_p = (TProfile*)TcorrFile1_3_p->Get("p_eff_pt");
  TdelR1_3_p = (TProfile*)TcorrFile1_3_p->Get("p_eff_rmin");

  TFile* TcorrFile3_8_p = new TFile(TcorrName3_8, "READ");
  Tcent3_8_p = (TProfile*)TcorrFile3_8_p->Get("p_eff_cent");
  TphiEta3_8_p = (TProfile2D*)TcorrFile3_8_p->Get("p_eff_acceptance");
  Tpt3_8_p = (TProfile*)TcorrFile3_8_p->Get("p_eff_pt");
  TdelR3_8_p = (TProfile*)TcorrFile3_8_p->Get("p_eff_rmin");

  TFile* TcorrFile8_100_p = new TFile(TcorrName8_100, "READ");
  Tcent8_100_p = (TProfile*)TcorrFile8_100_p->Get("p_eff_cent");
  TphiEta8_100_p = (TProfile2D*)TcorrFile8_100_p->Get("p_eff_acceptance");
  Tpt8_100_p = (TProfile*)TcorrFile8_100_p->Get("p_eff_pt");
  TdelR8_100_p = (TProfile*)TcorrFile8_100_p->Get("p_eff_rmin");


  TFile *outFile = new TFile(outName, "RECREATE");

  InitDiJetSkim(montecarlo);

  HiForest *c = new HiForest(listOfFiles[0].data(), "Forest", cType, montecarlo);

  c->InitTree();

  /*  if(cType == cPbPb)
    c->GetEnergyScaleTable((char*)"photonEnergyScaleTable_lowPt_v6.root");
  */

  c->LoadNoTrees();
  c->hasSkimTree = true;
  c->hasTrackTree = true;
  c->hasEvtTree = true;
  c->hasAkPu3JetTree = true;
  c->hasAkPu3CaloJetTree = true;
  c->hasAkVs3PFJetTree = true;
  c->hasAkVs3CaloJetTree = true;

  if(montecarlo)
    c->hasGenParticleTree = true;

  Long64_t nentries = c->GetEntries();

  Int_t totEv = 0;
  Int_t selectCut = 0;

  Int_t TLeadJtPtCut = 0;
  Int_t TSubLeadJtPtCut = 0;
  Int_t TDelPhiCut = 0;
  Int_t TJtEtaCut = 0;

  Int_t PFLeadJtPtCut = 0;
  Int_t PFSubLeadJtPtCut = 0;
  Int_t PFDelPhiCut = 0;
  Int_t PFJtEtaCut = 0;

  Int_t VsPFLeadJtPtCut = 0;
  Int_t VsPFSubLeadJtPtCut = 0;
  Int_t VsPFDelPhiCut = 0;
  Int_t VsPFJtEtaCut = 0;

  Int_t CaloLeadJtPtCut = 0;
  Int_t CaloSubLeadJtPtCut = 0;
  Int_t CaloDelPhiCut = 0;
  Int_t CaloJtEtaCut = 0;

  Int_t TTotTrk = 0;
  Int_t TTrkEtaCut = 0;
  Int_t TTrkPtCut = 0;
  Int_t TPurityCut = 0;
  Int_t TTrkDzCut = 0;
  Int_t TTrkDxyCut = 0;
  Int_t TTrkPtErrorCut = 0;

  Int_t TTotGen = 0;
  Int_t TGenEtaCut = 0;
  Int_t TGenPtCut = 0;
  Int_t TGenChgCut = 0;

  Int_t PFTotTrk = 0;
  Int_t PFTrkEtaCut = 0;
  Int_t PFTrkPtCut = 0;
  Int_t PFPurityCut = 0;
  Int_t PFTrkDzCut = 0;
  Int_t PFTrkDxyCut = 0;
  Int_t PFTrkPtErrorCut = 0;

  Int_t PFTotGen = 0;
  Int_t PFGenEtaCut = 0;
  Int_t PFGenPtCut = 0;
  Int_t PFGenChgCut = 0;

  Int_t CaloTotTrk = 0;
  Int_t CaloTrkEtaCut = 0;
  Int_t CaloTrkPtCut = 0;
  Int_t CaloPurityCut = 0;
  Int_t CaloTrkDzCut = 0;
  Int_t CaloTrkDxyCut = 0;
  Int_t CaloTrkPtErrorCut = 0;

  Int_t CaloTotGen = 0;
  Int_t CaloGenEtaCut = 0;
  Int_t CaloGenPtCut = 0;
  Int_t CaloGenChgCut = 0;

  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    c->GetEntry(jentry);

    totEv++;

    Bool_t TEventPass = false;
    Bool_t PFEventPass = false;
    Bool_t VsPFEventPass = false;
    Bool_t CaloEventPass = false;
    Bool_t VsCaloEventPass = false;

    if(jentry%1000 == 0)
      std::cout << jentry << std::endl;

    if(!c->selectEvent()){
      selectCut++;
      continue;
    }

    InitJetVar(montecarlo);

    //particle flow

    Int_t leadJtIndex = -1;
    Int_t subLeadJtIndex = -1;
    PFLeadJtPt_ = subLeadJtPtCut;
    PFSubLeadJtPt_ = subLeadJtPtCut;
    for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.nref; jtEntry++){
      if(c->akPu3PF.jtpt[jtEntry] > leadJtPtCut && c->akPu3PF.jtpt[jtEntry] > PFLeadJtPt_){
	subLeadJtIndex = leadJtIndex;
	PFSubLeadJtPt_ = PFLeadJtPt_;
	leadJtIndex = jtEntry;
	PFLeadJtPt_ = c->akPu3PF.jtpt[jtEntry];
      }
      else if(c->akPu3PF.jtpt[jtEntry] > PFSubLeadJtPt_){
	subLeadJtIndex = jtEntry;
	PFSubLeadJtPt_ = c->akPu3PF.jtpt[jtEntry];
      }
    }

    if(leadJtIndex < 0){
      PFLeadJtPtCut++;
      PFLeadJtPt_ = -10;
      PFSubLeadJtPt_ = -10;
    }
    else if(subLeadJtIndex < 0){
      PFSubLeadJtPtCut++;
      PFSubLeadJtPt_ = -10;
    }
    else if(getAbsDphi(c->akPu3PF.jtphi[leadJtIndex], c->akPu3PF.jtphi[subLeadJtIndex]) < jtDelPhiCut){
      PFDelPhiCut++;
    }
    else if(TMath::Abs(c->akPu3PF.jteta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akPu3PF.jteta[subLeadJtIndex]) > jtEtaCut){
      PFJtEtaCut++;
    }
    else{
      PFLeadJtPhi_ = c->akPu3PF.jtphi[leadJtIndex];
      PFSubLeadJtPhi_ = c->akPu3PF.jtphi[subLeadJtIndex];
      PFLeadJtEta_ = c->akPu3PF.jteta[leadJtIndex];
      PFSubLeadJtEta_ = c->akPu3PF.jteta[subLeadJtIndex];

      PFJtAsymm_ = (PFLeadJtPt_ - PFSubLeadJtPt_)/(PFLeadJtPt_ + PFSubLeadJtPt_);

      PFEventPass = true;

      recoPFSet_ = true;
    }


    //Vs PF jet

    leadJtIndex = -1;
    subLeadJtIndex = -1;
    VsPFLeadJtPt_ = subLeadJtPtCut;
    VsPFSubLeadJtPt_ = subLeadJtPtCut;
    for(Int_t jtEntry = 0; jtEntry < c->akVs3PF.nref; jtEntry++){
      if(c->akVs3PF.jtpt[jtEntry] > leadJtPtCut && c->akVs3PF.jtpt[jtEntry] > VsPFLeadJtPt_){
	subLeadJtIndex = leadJtIndex;
	VsPFSubLeadJtPt_ = VsPFLeadJtPt_;
	leadJtIndex = jtEntry;
	VsPFLeadJtPt_ = c->akVs3PF.jtpt[jtEntry];
      }
      else if(c->akVs3PF.jtpt[jtEntry] > VsPFSubLeadJtPt_){
	subLeadJtIndex = jtEntry;
	VsPFSubLeadJtPt_ = c->akVs3PF.jtpt[jtEntry];
      }
    }

    if(leadJtIndex < 0){
      VsPFLeadJtPtCut++;
      VsPFLeadJtPt_ = -10;
      VsPFSubLeadJtPt_ = -10;
    }
    else if(subLeadJtIndex < 0){
      VsPFSubLeadJtPtCut++;
      VsPFSubLeadJtPt_ = -10;
    }
    else if(getAbsDphi(c->akVs3PF.jtphi[leadJtIndex], c->akVs3PF.jtphi[subLeadJtIndex]) < jtDelPhiCut){
      VsPFDelPhiCut++;
    }
    else if(TMath::Abs(c->akVs3PF.jteta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akVs3PF.jteta[subLeadJtIndex]) > jtEtaCut){
      VsPFJtEtaCut++;
    }
    else{
      VsPFLeadJtPhi_ = c->akVs3PF.jtphi[leadJtIndex];
      VsPFSubLeadJtPhi_ = c->akVs3PF.jtphi[subLeadJtIndex];
      VsPFLeadJtEta_ = c->akVs3PF.jteta[leadJtIndex];
      VsPFSubLeadJtEta_ = c->akVs3PF.jteta[subLeadJtIndex];

      VsPFJtAsymm_ = (VsPFLeadJtPt_ - VsPFSubLeadJtPt_)/(VsPFLeadJtPt_ + VsPFSubLeadJtPt_);

      VsPFEventPass = true;

      recoVsPFSet_ = true;
    }


    //calo jet

    leadJtIndex = -1;
    subLeadJtIndex = -1;
    CaloLeadJtPt_ = subLeadJtPtCut;
    CaloSubLeadJtPt_ = subLeadJtPtCut;
    for(Int_t jtEntry = 0; jtEntry < c->akPu3Calo.nref; jtEntry++){
      if(c->akPu3Calo.jtpt[jtEntry] > leadJtPtCut && c->akPu3Calo.jtpt[jtEntry] > CaloLeadJtPt_){
	subLeadJtIndex = leadJtIndex;
	CaloSubLeadJtPt_ = CaloLeadJtPt_;
	leadJtIndex = jtEntry;
	CaloLeadJtPt_ = c->akPu3Calo.jtpt[jtEntry];
      }
      else if(c->akPu3Calo.jtpt[jtEntry] > CaloSubLeadJtPt_){
	subLeadJtIndex = jtEntry;
	CaloSubLeadJtPt_ = c->akPu3Calo.jtpt[jtEntry];
      }
    }

    if(leadJtIndex < 0){
      CaloLeadJtPtCut++;
      CaloLeadJtPt_ = -10;
      CaloSubLeadJtPt_ = -10;
    }
    else if(subLeadJtIndex < 0){
      CaloSubLeadJtPtCut++;
      CaloSubLeadJtPt_ = -10;
    }
    else if(getAbsDphi(c->akPu3Calo.jtphi[leadJtIndex], c->akPu3Calo.jtphi[subLeadJtIndex]) < jtDelPhiCut){
      CaloDelPhiCut++;
    }
    else if(TMath::Abs(c->akPu3Calo.jteta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akPu3Calo.jteta[subLeadJtIndex]) > jtEtaCut){
      CaloJtEtaCut++;
    }
    else{
      CaloLeadJtPhi_ = c->akPu3Calo.jtphi[leadJtIndex];
      CaloSubLeadJtPhi_ = c->akPu3Calo.jtphi[subLeadJtIndex];
      CaloLeadJtEta_ = c->akPu3Calo.jteta[leadJtIndex];
      CaloSubLeadJtEta_ = c->akPu3Calo.jteta[subLeadJtIndex];

      CaloJtAsymm_ = (CaloLeadJtPt_ - CaloSubLeadJtPt_)/(CaloLeadJtPt_ + CaloSubLeadJtPt_);

      CaloEventPass = true;

      recoCaloSet_ = true;
    }



    //Vs calo jet

    leadJtIndex = -1;
    subLeadJtIndex = -1;
    VsCaloLeadJtPt_ = subLeadJtPtCut;
    VsCaloSubLeadJtPt_ = subLeadJtPtCut;
    for(Int_t jtEntry = 0; jtEntry < c->akVs3Calo.nref; jtEntry++){
      if(c->akVs3Calo.jtpt[jtEntry] > leadJtPtCut && c->akVs3Calo.jtpt[jtEntry] > VsCaloLeadJtPt_){
	subLeadJtIndex = leadJtIndex;
	VsCaloSubLeadJtPt_ = VsCaloLeadJtPt_;
	leadJtIndex = jtEntry;
	VsCaloLeadJtPt_ = c->akVs3Calo.jtpt[jtEntry];
      }
      else if(c->akVs3Calo.jtpt[jtEntry] > VsCaloSubLeadJtPt_){
	subLeadJtIndex = jtEntry;
	VsCaloSubLeadJtPt_ = c->akVs3Calo.jtpt[jtEntry];
      }
    }

    if(leadJtIndex < 0){
      VsCaloLeadJtPtCut++;
      VsCaloLeadJtPt_ = -10;
      VsCaloSubLeadJtPt_ = -10;
    }
    else if(subLeadJtIndex < 0){
      VsCaloSubLeadJtPtCut++;
      VsCaloSubLeadJtPt_ = -10;
    }
    else if(getAbsDphi(c->akVs3Calo.jtphi[leadJtIndex], c->akVs3Calo.jtphi[subLeadJtIndex]) < jtDelPhiCut){
      VsCaloDelPhiCut++;
    }
    else if(TMath::Abs(c->akVs3Calo.jteta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akVs3Calo.jteta[subLeadJtIndex]) > jtEtaCut){
      VsCaloJtEtaCut++;
    }
    else{
      VsCaloLeadJtPhi_ = c->akVs3Calo.jtphi[leadJtIndex];
      VsCaloSubLeadJtPhi_ = c->akVs3Calo.jtphi[subLeadJtIndex];
      VsCaloLeadJtEta_ = c->akVs3Calo.jteta[leadJtIndex];
      VsCaloSubLeadJtEta_ = c->akVs3Calo.jteta[subLeadJtIndex];

      VsCaloJtAsymm_ = (VsCaloLeadJtPt_ - VsCaloSubLeadJtPt_)/(VsCaloLeadJtPt_ + VsCaloSubLeadJtPt_);

      VsCaloEventPass = true;

      recoVsCaloSet_ = true;
    }


    //truth

    leadJtIndex = -1;
    subLeadJtIndex = -1;
    TLeadJtPt_ = subLeadJtPtCut;
    TSubLeadJtPt_ = subLeadJtPtCut;
    for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.nref; jtEntry++){
      if(c->akPu3PF.refpt[jtEntry] > leadJtPtCut && c->akPu3PF.refpt[jtEntry] > TLeadJtPt_){
	subLeadJtIndex = leadJtIndex;
	TSubLeadJtPt_ = TLeadJtPt_;
	leadJtIndex = jtEntry;
	TLeadJtPt_ = c->akPu3PF.refpt[jtEntry];
      }
      else if(c->akPu3PF.refpt[jtEntry] > TSubLeadJtPt_){
	subLeadJtIndex = jtEntry;
	TSubLeadJtPt_ = c->akPu3PF.refpt[jtEntry];
      }
    }

    if(leadJtIndex < 0){
      TLeadJtPtCut++;
      TLeadJtPt_ = -10;
      TSubLeadJtPt_ = -10;
    }
    else if(subLeadJtIndex < 0){
      TSubLeadJtPtCut++;
      TSubLeadJtPt_ = -10;
    }
    else if(getAbsDphi(c->akPu3PF.refphi[leadJtIndex], c->akPu3PF.refphi[subLeadJtIndex]) < jtDelPhiCut){
      TDelPhiCut++;
    }
    else if(TMath::Abs(c->akPu3PF.refeta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akPu3PF.refeta[subLeadJtIndex]) > jtEtaCut){
      TJtEtaCut++;
    }
    else{
      TLeadJtPhi_ = c->akPu3PF.refphi[leadJtIndex];
      TSubLeadJtPhi_ = c->akPu3PF.refphi[subLeadJtIndex];
      TLeadJtEta_ = c->akPu3PF.refeta[leadJtIndex];
      TSubLeadJtEta_ = c->akPu3PF.refeta[subLeadJtIndex];

      TJtAsymm_ = (TLeadJtPt_ - TSubLeadJtPt_)/(TLeadJtPt_ + TSubLeadJtPt_);

      TEventPass = true;

      truthSet_ = true;
    }

    if(TEventPass == false && PFEventPass == false && CaloEventPass == false && VsPFEventPass == false && VsCaloEventPass == false)
      continue;

    run_ = c->evt.run;
    evt_ = c->akPu3PF.evt;
    lumi_ = c->evt.lumi;
    hiBin_ = c->evt.hiBin;

    //Iterate over tracks

    nTrk_ = 0;

    InitProjPerp(montecarlo);

    Tracks trkCollection;
    trkCollection = c->track;

    for(Int_t trkEntry = 0; trkEntry < trkCollection.nTrk; trkEntry++){
      if(TEventPass)	TTotTrk++;
      if(PFEventPass)	PFTotTrk++;
      if(CaloEventPass)	CaloTotTrk++;

      if(TMath::Abs(trkCollection.trkEta[trkEntry]) > 2.4){
	if(TEventPass)	  TTrkEtaCut++;
	if(PFEventPass)	  PFTrkEtaCut++;
	if(CaloEventPass)	  CaloTrkEtaCut++;

	continue;
      }
      
      if(trkCollection.trkPt[trkEntry] < 0.5){
	if(TEventPass)	  TTrkPtCut++;
	if(PFEventPass)	  PFTrkPtCut++;
	if(CaloEventPass)	  CaloTrkPtCut++;

	continue;
      }

      if(!trkCollection.highPurity[trkEntry]){ //Note highPuritySetWithPV seems to be wrong cut, creates diff. bet. truth and reco
	if(TEventPass)	  TPurityCut++;
	if(PFEventPass)	  PFPurityCut++;
	if(CaloEventPass)	  CaloPurityCut++;

	continue;
      }

      if(TMath::Abs(trkCollection.trkDz1[trkEntry]/trkCollection.trkDzError1[trkEntry]) > 3){
	if(TEventPass)	  TTrkDzCut++;
	if(PFEventPass)	  PFTrkDzCut++;
	if(CaloEventPass)	  CaloTrkDzCut++;

	continue;
      }

      if(TMath::Abs(trkCollection.trkDxy1[trkEntry]/trkCollection.trkDxyError1[trkEntry]) > 3){
	if(TEventPass)	  TTrkDxyCut++;
	if(PFEventPass)	  PFTrkDxyCut++;
	if(CaloEventPass)	  CaloTrkDxyCut++;

	continue;
      }

      if(trkCollection.trkPtError[trkEntry]/trkCollection.trkPt[trkEntry] > 0.1){
	if(TEventPass)	  TTrkPtErrorCut++;
	if(PFEventPass)	  PFTrkPtErrorCut++;
	if(CaloEventPass)	  CaloTrkPtErrorCut++;

	continue;
      }


      trkPt_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPhi_[nTrk_] = trkCollection.trkPhi[trkEntry];
      trkEta_[nTrk_] = trkCollection.trkEta[trkEntry];
      
      if(PFEventPass)
	trkPtPF_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));
      else
	trkPtPF_[nTrk_] = 0;

      if(CaloEventPass)
	trkPtCalo_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));
      else
	trkPtCalo_[nTrk_] = 0;

      if(montecarlo && TEventPass)
	trkPtT_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));
      else
	trkPtT_[nTrk_] = 0;


      if(PFEventPass)
	trkPFLeadDelPhi_[nTrk_] = getAbsDphi(PFLeadJtPhi_, trkCollection.trkPhi[trkEntry]);
      else
	trkPFLeadDelPhi_[nTrk_] = -10;

      if(CaloEventPass)
	trkCaloLeadDelPhi_[nTrk_] = getAbsDphi(CaloLeadJtPhi_, trkCollection.trkPhi[trkEntry]);
      else
	trkCaloLeadDelPhi_[nTrk_] = -10;


      if(PFEventPass){
	rPFImbProjF_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));
	rPFImbPerpF_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));
	if(trkCollection.trkPt[trkEntry] > 8){
	  rPFImbProjH_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));
	  rPFImbPerpH_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));

	  rPFImbProj8_100_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));
	}
	else{
	  rPFImbProjL_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));
	  rPFImbPerpL_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));

	  if(trkCollection.trkPt[trkEntry] < 1){
	    rPFImbProj0_1_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 2){
	    rPFImbProj1_2_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 4){
	    rPFImbProj2_4_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 8){
	    rPFImbProj4_8_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));
	  }

	}
      }

      if(CaloEventPass){
	rCaloImbProjF_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));
	rCaloImbPerpF_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));
	if(trkCollection.trkPt[trkEntry] > 8){
	  rCaloImbProjH_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));
	  rCaloImbPerpH_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));

	  rCaloImbProj8_100_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));
	}
	else{
	  rCaloImbProjL_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));
	  rCaloImbPerpL_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));

	  if(trkCollection.trkPt[trkEntry] < 1){
	    rCaloImbProj0_1_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 2){
	    rCaloImbProj1_2_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 4){
	    rCaloImbProj2_4_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 8){
	    rCaloImbProj4_8_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));
	  }

	}
      }

      if(montecarlo && TEventPass){
	rTImbProjF_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));
	rTImbPerpF_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));
	if(trkCollection.trkPt[trkEntry] > 8){
	  rTImbProjH_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));
	  rTImbPerpH_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));

	  rTImbProj8_100_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));
	}
	else{
	  rTImbProjL_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));
	  rTImbPerpL_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));

	  if(trkCollection.trkPt[trkEntry] < 1){
	    rTImbProj0_1_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 2){
	    rTImbProj1_2_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 4){
	    rTImbProj2_4_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 8){
	    rTImbProj4_8_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));
	  }

	}
      }

    
      nTrk_++;
      if(nTrk_ > MAXTRKS - 1){
	printf("ERROR: Trk arrays not large enough.\n");
	return(1);
      }
    }

    //Get trk rmin/corrections for 3 Jet subsets

    for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
      trkRMinPF_[trkEntry] = 6;
      trkRMinCalo_[trkEntry] = 6;
      trkRMinT_[trkEntry] = 6;

      trkPtCorrCalo_[trkEntry] = trkPt_[trkEntry];

      for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.nref; jtEntry++){
	if(c->akPu3PF.jtpt[jtEntry] < 30 || TMath::Abs(c->akPu3PF.jteta[jtEntry]) > 2.0)
	  continue;

	if(trkRMinPF_[trkEntry] > getDR(trkEta_[trkEntry], trkPhi_[trkEntry], c->akPu3PF.jteta[jtEntry], c->akPu3PF.jtphi[jtEntry]))
	  trkRMinPF_[trkEntry] = getDR(trkEta_[trkEntry], trkPhi_[trkEntry], c->akPu3PF.jteta[jtEntry], c->akPu3PF.jtphi[jtEntry]);
      }

      if(trkPt_[trkEntry] > 0.5 && trkPt_[trkEntry] < 1.0){
	trkPtCorrPF_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent0_1_p, PFphiEta0_1_p, PFpt0_1_p, PFdelR0_1_p);
	trkPtFactPF_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent0_1_p, PFphiEta0_1_p, PFpt0_1_p, PFdelR0_1_p);
      }
      else if(trkPt_[trkEntry] > 1 && trkPt_[trkEntry] < 3){
	trkPtCorrPF_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent1_3_p, PFphiEta1_3_p, PFpt1_3_p, PFdelR1_3_p);
	trkPtFactPF_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent1_3_p, PFphiEta1_3_p, PFpt1_3_p, PFdelR1_3_p);
      }
      else if(trkPt_[trkEntry] > 3 && trkPt_[trkEntry] < 8){
	trkPtCorrPF_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent3_8_p, PFphiEta3_8_p, PFpt3_8_p, PFdelR3_8_p);
	trkPtFactPF_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent3_8_p, PFphiEta3_8_p, PFpt3_8_p, PFdelR3_8_p);
      }
      else if(trkPt_[trkEntry] > 8 && trkPt_[trkEntry] < 100){
	trkPtCorrPF_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent8_100_p, PFphiEta8_100_p, PFpt8_100_p, PFdelR8_100_p);
	trkPtFactPF_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent8_100_p, PFphiEta8_100_p, PFpt8_100_p, PFdelR8_100_p);
      }
      else{
	trkPtCorrPF_[trkEntry] = trkPt_[trkEntry];
	trkPtFactPF_[trkEntry] = 1;
      }

      //Truth 

      for(Int_t refEntry = 0; refEntry < c->akPu3PF.nref; refEntry++){
	if(c->akPu3PF.refpt[refEntry] < 30 || TMath::Abs(c->akPu3PF.refeta[refEntry]) > 2.0)
	  continue;

	if(trkRMinT_[trkEntry] > getDR(trkEta_[trkEntry], trkPhi_[trkEntry], c->akPu3PF.refeta[refEntry], c->akPu3PF.refphi[refEntry]))
	  trkRMinT_[trkEntry] = getDR(trkEta_[trkEntry], trkPhi_[trkEntry], c->akPu3PF.refeta[refEntry], c->akPu3PF.refphi[refEntry]);
      }

      if(trkPt_[trkEntry] > 0.5 && trkPt_[trkEntry] < 1.0){
	trkPtCorrT_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinT_[trkEntry], Tcent0_1_p, TphiEta0_1_p, Tpt0_1_p, TdelR0_1_p);
	trkPtFactT_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinT_[trkEntry], Tcent0_1_p, TphiEta0_1_p, Tpt0_1_p, TdelR0_1_p);
      }
      else if(trkPt_[trkEntry] > 1 && trkPt_[trkEntry] < 3){
	trkPtCorrT_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinT_[trkEntry], Tcent1_3_p, TphiEta1_3_p, Tpt1_3_p, TdelR1_3_p);
	trkPtFactT_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinT_[trkEntry], Tcent1_3_p, TphiEta1_3_p, Tpt1_3_p, TdelR1_3_p);
      }
      else if(trkPt_[trkEntry] > 3 && trkPt_[trkEntry] < 8){
	trkPtCorrT_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinT_[trkEntry], Tcent3_8_p, TphiEta3_8_p, Tpt3_8_p, TdelR3_8_p);
	trkPtFactT_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinT_[trkEntry], Tcent3_8_p, TphiEta3_8_p, Tpt3_8_p, TdelR3_8_p);
      }
      else if(trkPt_[trkEntry] > 8 && trkPt_[trkEntry] < 100){
	trkPtCorrT_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinT_[trkEntry], Tcent8_100_p, TphiEta8_100_p, Tpt8_100_p, TdelR8_100_p);
	trkPtFactT_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinT_[trkEntry], Tcent8_100_p, TphiEta8_100_p, Tpt8_100_p, TdelR8_100_p);
      }
      else{
	trkPtCorrT_[trkEntry] = trkPt_[trkEntry];
	trkPtFactT_[trkEntry] = 1;
      }
    }

    //Apply corrections to appropriate subsets


    if(PFEventPass){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	rPFImbProjFCorr_ += -trkPtCorrPF_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], PFLeadJtPhi_));             
        rPFImbPerpFCorr_ += -trkPtCorrPF_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], PFLeadJtPhi_));             
        if(trkPt_[trkEntry] > 8){                                                                     
          rPFImbProjHCorr_ += -trkPtCorrPF_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], PFLeadJtPhi_));           
          rPFImbPerpHCorr_ += -trkPtCorrPF_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], PFLeadJtPhi_));           

          rPFImbProj8_100Corr_ += -trkPtCorrPF_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], PFLeadJtPhi_));           
        }                                                                                                 
        else{                                                                                             
          rPFImbProjLCorr_ += -trkPtCorrPF_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], PFLeadJtPhi_));           
          rPFImbPerpLCorr_ += -trkPtCorrPF_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], PFLeadJtPhi_));           

	  if(trkPt_[trkEntry] < 1){
	    rPFImbProj0_1Corr_ += -trkPtCorrPF_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], PFLeadJtPhi_));
	  }
	  else if(trkPt_[trkEntry] < 2){
	    rPFImbProj1_2Corr_ += -trkPtCorrPF_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], PFLeadJtPhi_));
	  }
	  else if(trkPt_[trkEntry] < 4){
	    rPFImbProj2_4Corr_ += -trkPtCorrPF_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], PFLeadJtPhi_));
	  }
	  else if(trkPt_[trkEntry] < 8){
	    rPFImbProj4_8Corr_ += -trkPtCorrPF_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], PFLeadJtPhi_));
	  }
        }                                                                                                 
      }
    }

    if(CaloEventPass){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	rCaloImbProjFCorr_ += -trkPtCorrCalo_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], CaloLeadJtPhi_));             
        rCaloImbPerpFCorr_ += -trkPtCorrCalo_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], CaloLeadJtPhi_));             
        if(trkPt_[trkEntry] > 8){                                                                     
          rCaloImbProjHCorr_ += -trkPtCorrCalo_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], CaloLeadJtPhi_));           
          rCaloImbPerpHCorr_ += -trkPtCorrCalo_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], CaloLeadJtPhi_));           

          rCaloImbProj8_100Corr_ += -trkPtCorrCalo_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], CaloLeadJtPhi_));           
        }                                                                                                 
        else{                                                                                             
          rCaloImbProjLCorr_ += -trkPtCorrCalo_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], CaloLeadJtPhi_));           
          rCaloImbPerpLCorr_ += -trkPtCorrCalo_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], CaloLeadJtPhi_));           

	  if(trkPt_[trkEntry] < 1){
	    rCaloImbProj0_1Corr_ += -trkPtCorrCalo_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], CaloLeadJtPhi_));
	  }
	  else if(trkPt_[trkEntry] < 2){
	    rCaloImbProj1_2Corr_ += -trkPtCorrCalo_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], CaloLeadJtPhi_));
	  }
	  else if(trkPt_[trkEntry] < 4){
	    rCaloImbProj2_4Corr_ += -trkPtCorrCalo_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], CaloLeadJtPhi_));
	  }
	  else if(trkPt_[trkEntry] < 8){
	    rCaloImbProj4_8Corr_ += -trkPtCorrCalo_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], CaloLeadJtPhi_));
	  }

        }                                                                                                 
      }
    }

    if(montecarlo && TEventPass){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
        rTImbProjFCorr_ += -trkPtCorrT_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], TLeadJtPhi_));            
        rTImbPerpFCorr_ += -trkPtCorrT_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], TLeadJtPhi_));            
        if(trkPt_[trkEntry] > 8){                                                                     
          rTImbProjHCorr_ += -trkPtCorrT_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], TLeadJtPhi_));          
          rTImbPerpHCorr_ += -trkPtCorrT_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], TLeadJtPhi_));          

          rTImbProj8_100Corr_ += -trkPtCorrT_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], TLeadJtPhi_));          
        }                                                                                                 
        else{                                                                                             
          rTImbProjLCorr_ += -trkPtCorrT_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], TLeadJtPhi_));          
          rTImbPerpLCorr_ += -trkPtCorrT_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], TLeadJtPhi_));          

	  if(trkPt_[trkEntry] < 1){
	    rTImbProj0_1Corr_ += -trkPtCorrT_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], TLeadJtPhi_));          
	  }
	  else if(trkPt_[trkEntry] < 2){
	    rTImbProj1_2Corr_ += -trkPtCorrT_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], TLeadJtPhi_));          
	  }
	  else if(trkPt_[trkEntry] < 4){
	    rTImbProj2_4Corr_ += -trkPtCorrT_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], TLeadJtPhi_));          
	  }
	  else if(trkPt_[trkEntry] < 8){
	    rTImbProj4_8Corr_ += -trkPtCorrT_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], TLeadJtPhi_));          
	  }

        }                                                                                                 
      }
    }

  
    if(montecarlo){
      //Iterate over truth

      nGen_ = 0;

      GenParticles genCollection;
      genCollection = c->genparticle;

      for(Int_t genEntry = 0; genEntry < genCollection.mult; genEntry++){
	if(TEventPass)	  TTotGen++;
	if(PFEventPass)	  PFTotGen++;
	if(CaloEventPass)	  CaloTotGen++;
	  
	if(genCollection.chg[genEntry] == 0){
	  if(TEventPass)	    TGenChgCut++;
	  if(PFEventPass)	    PFGenChgCut++;
	  if(CaloEventPass)	    CaloGenChgCut++;
	  
	  continue;
	}
	  
	if(TMath::Abs(genCollection.eta[genEntry]) > 2.4){
	  if(TEventPass)	    TGenEtaCut++;
	  if(PFEventPass)	    PFGenEtaCut++;
	  if(CaloEventPass)	    CaloGenEtaCut++;

	  continue;
	}
	
	if(genCollection.pt[genEntry] < 0.5){
	  if(TEventPass)	    TGenPtCut++;
	  if(PFEventPass)	    PFGenPtCut++;
	  if(CaloEventPass)	    CaloGenPtCut++;

	  continue;
	}

	genPt_[nGen_] = genCollection.pt[genEntry];
	genPhi_[nGen_] = genCollection.phi[genEntry];
	genEta_[nGen_] = genCollection.eta[genEntry];

	if(PFEventPass)
          genPtPF_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));
	else
	  genPtPF_[nGen_] = 0;

	if(CaloEventPass)
          genPtCalo_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));
	else
	  genPtCalo_[nGen_] = 0;

	if(TEventPass)
          genPtT_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	else
	  genPtT_[nGen_] = 0;

	
	if(TEventPass)
	  genLeadDelPhi_[nGen_] = getAbsDphi(TLeadJtPhi_, genCollection.phi[genEntry]);
	else
	  genLeadDelPhi_[nGen_] = -10;


	if(TEventPass){
	  gTImbProjF_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	  gTImbPerpF_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	  if(genCollection.pt[genEntry] > 8){
	    gTImbProjH_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	    gTImbPerpH_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	    gTImbProj8_100_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	  }
	  else{
	    gTImbProjL_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	    gTImbPerpL_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));

	    if(genCollection.pt[genEntry] < 1){
	      gTImbProj0_1_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	    }
	    else if(genCollection.pt[genEntry] < 2){
	      gTImbProj1_2_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	    }
	    else if(genCollection.pt[genEntry] < 4){
	      gTImbProj2_4_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	    }
	    else if(genCollection.pt[genEntry] < 8){
	      gTImbProj4_8_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	    }
	  }
	}

	if(PFEventPass){
	  gPFImbProjF_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));
	  gPFImbPerpF_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));

	  if(genCollection.pt[genEntry] > 8){
            gPFImbProjH_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));
            gPFImbPerpH_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));
            gPFImbProj8_100_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));
          }
          else{
            gPFImbProjL_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));
            gPFImbPerpL_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));

            if(genCollection.pt[genEntry] < 1){
              gPFImbProj0_1_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));
            }
            else if(genCollection.pt[genEntry] < 2){
              gPFImbProj1_2_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));
            }
            else if(genCollection.pt[genEntry] < 4){
              gPFImbProj2_4_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));
            }
            else if(genCollection.pt[genEntry] < 8){
              gPFImbProj4_8_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));
            }
          }
	}

	if(CaloEventPass){
	  gCaloImbProjF_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));
	  gCaloImbPerpF_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));

	  if(genCollection.pt[genEntry] > 8){
            gCaloImbProjH_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));
            gCaloImbPerpH_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));
            gCaloImbProj8_100_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));
          }
          else{
            gCaloImbProjL_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));
            gCaloImbPerpL_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));

            if(genCollection.pt[genEntry] < 1){
              gCaloImbProj0_1_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));
            }
            else if(genCollection.pt[genEntry] < 2){
              gCaloImbProj1_2_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));
            }
            else if(genCollection.pt[genEntry] < 4){
              gCaloImbProj2_4_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));
            }
            else if(genCollection.pt[genEntry] < 8){
              gCaloImbProj4_8_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));
            }
          }
	}
	  
	nGen_++;
	if(nGen_ > MAXGEN - 1){
	  printf("ERROR: Gen arrays not large enough.\n");
	  return(1);
	}
      }
    }

    std::cout << "Event, totTrk, trkKept: " << jentry << ", " << trkCollection.nTrk << ", " << nTrk_ << std::endl;


    jetTree_p->Fill();
    trackTree_p->Fill();

    if(montecarlo)
      genTree_p->Fill();
  }

  std::cout << "totEv: " << totEv << std::endl;
  Int_t tempTot = totEv - selectCut;
  std::cout << "selectCut: " << tempTot << std::endl;

  std::cout << std::endl;
  tempTot = tempTot - TLeadJtPtCut;
  std::cout << "TJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - TSubLeadJtPtCut;
  std::cout << "TSubLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - TDelPhiCut;
  std::cout << "TDelPhiCut: " << tempTot << std::endl;
  tempTot = tempTot - TJtEtaCut;
  std::cout << "TJtEtaCut: " << tempTot << std::endl;

  std::cout << std::endl;
  tempTot = totEv - selectCut - PFLeadJtPtCut;
  std::cout << "PFLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - PFSubLeadJtPtCut;
  std::cout << "PFSubLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - PFDelPhiCut;
  std::cout << "PFDelPhiCut: " << tempTot << std::endl;
  tempTot = tempTot - PFJtEtaCut;
  std::cout << "PFJtEtaCut: " << tempTot << std::endl;

  std::cout << std::endl;
  tempTot = totEv - selectCut - CaloLeadJtPtCut;
  std::cout << "CaloLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - CaloSubLeadJtPtCut;
  std::cout << "CaloSubLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - CaloDelPhiCut;
  std::cout << "CaloDelPhiCut: " << tempTot << std::endl;
  tempTot = tempTot - CaloJtEtaCut;
  std::cout << "CaloJtEtaCut: " << tempTot << std::endl;

  std::cout << std::endl;
  std::cout << "TTotTrk: " << TTotTrk << std::endl;
  tempTot = TTotTrk - TTrkEtaCut;
  std::cout << "TTrkEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - TTrkPtCut;
  std::cout << "TTrkPtCut: " << tempTot << std::endl;
  tempTot = tempTot - TPurityCut;
  std::cout << "TPurityCut: " << tempTot << std::endl;
  tempTot = tempTot - TTrkDzCut;
  std::cout << "TTrkDzCut: " << tempTot << std::endl;
  tempTot = tempTot - TTrkDxyCut;
  std::cout << "TTrkDxyCut: " << tempTot << std::endl;
  tempTot = tempTot - TTrkPtErrorCut;
  std::cout << "TTrkPtErrorCut: " << tempTot << std::endl;

  std::cout << std::endl;
  std::cout << "TTotGen: " << TTotGen << std::endl;
  tempTot = TTotGen - TGenChgCut;
  std::cout << "TGenChgCut: " << tempTot << std::endl;
  tempTot = tempTot - TGenEtaCut;
  std::cout << "TGenEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - TGenPtCut;
  std::cout << "TGenPtCut: " << tempTot << std::endl;
  


  std::cout << std::endl;
  std::cout << "PFTotTrk: " << PFTotTrk << std::endl;
  tempTot = PFTotTrk - PFTrkEtaCut;
  std::cout << "PFTrkEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - PFTrkPtCut;
  std::cout << "PFTrkPtCut: " << tempTot << std::endl;
  tempTot = tempTot - PFPurityCut;
  std::cout << "PFPurityCut: " << tempTot << std::endl;
  tempTot = tempTot - PFTrkDzCut;
  std::cout << "PFTrkDzCut: " << tempTot << std::endl;
  tempTot = tempTot - PFTrkDxyCut;
  std::cout << "PFTrkDxyCut: " << tempTot << std::endl;
  tempTot = tempTot - PFTrkPtErrorCut;
  std::cout << "PFTrkPtErrorCut: " << tempTot << std::endl;

  std::cout << std::endl;
  std::cout << "PFTotGen: " << PFTotGen << std::endl;
  tempTot = PFTotGen - PFGenChgCut;
  std::cout << "PFGenChgCut: " << tempTot << std::endl;
  tempTot = tempTot - PFGenEtaCut;
  std::cout << "PFGenEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - PFGenPtCut;
  std::cout << "PFGenPtCut: " << tempTot << std::endl;


  std::cout << std::endl;
  std::cout << "CaloTotTrk: " << CaloTotTrk << std::endl;
  tempTot = CaloTotTrk - CaloTrkEtaCut;
  std::cout << "CaloTrkEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - CaloTrkPtCut;
  std::cout << "CaloTrkPtCut: " << tempTot << std::endl;
  tempTot = tempTot - CaloPurityCut;
  std::cout << "CaloPurityCut: " << tempTot << std::endl;
  tempTot = tempTot - CaloTrkDzCut;
  std::cout << "CaloTrkDzCut: " << tempTot << std::endl;
  tempTot = tempTot - CaloTrkDxyCut;
  std::cout << "CaloTrkDxyCut: " << tempTot << std::endl;
  tempTot = tempTot - CaloTrkPtErrorCut;
  std::cout << "CaloTrkPtErrorCut: " << tempTot << std::endl;

  std::cout << std::endl;
  std::cout << "CaloTotGen: " << CaloTotGen << std::endl;
  tempTot = CaloTotGen - CaloGenChgCut;
  std::cout << "CaloGenChgCut: " << tempTot << std::endl;
  tempTot = tempTot - CaloGenEtaCut;
  std::cout << "CaloGenEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - CaloGenPtCut;
  std::cout << "CaloGenPtCut: " << tempTot << std::endl;

  outFile->cd();
  jetTree_p->Write();
  trackTree_p->Write();
  genTree_p->Write();
  outFile->Close();

  printf("Done.\n");
  return(0);
}


collisionType getCType(sampleType sType)
{
  switch (sType)
    {
    case kPPDATA:
    case kPPMC:
      return cPP;
    case kPADATA:
    case kPAMC:
      return cPPb;
    case kHIDATA:
    case kHIMC:
      return cPbPb;
    }
  return cPbPb; //probably a bad guess
}


int main(int argc, char *argv[])
{
  if(argc != 4)
    {
      std::cout << "Usage: sortForest <inputFile> <MCBool> <outputFile>" << std::endl;
      return 1;
    }

  int rStatus = -1;

  rStatus = makeDiJetTTree(argv[1], sampleType(atoi(argv[2])), argv[3]);

  return rStatus;
}
