//=============================================
// Author: Chris McGinn
// 
// DiJet Analysis Class (MC)
//
//=============================================

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "../../HiForestAnalysis/hiForest.h"
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

TProfile* cent0_1_p;
TProfile* cent1_3_p;
TProfile* cent3_8_p;
TProfile* cent8_100_p;
TProfile* cent100_300_p;

TProfile2D* phiEta0_1_p;
TProfile2D* phiEta1_3_p;
TProfile2D* phiEta3_8_p;
TProfile2D* phiEta8_100_p;
TProfile2D* phiEta100_300_p;

TProfile* pt0_1_p;
TProfile* pt1_3_p;
TProfile* pt3_8_p;
TProfile* pt8_100_p;
TProfile* pt100_300_p;

TProfile* delR0_1_p;
TProfile* delR1_3_p;
TProfile* delR3_8_p;
TProfile* delR8_100_p;
TProfile* delR100_300_p;

int makeDiJetTTree(string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_CFMSKIM.root", const char *corrName0_1 = "eff_pt0_1_step_cent3accept3pt3rmin2.root", const char *corrName1_3 = "eff_pt1_3_step_cent3accept3pt3rmin2.root", const char *corrName3_8 = "eff_pt3_8_step_cent3accept3pt3rmin2.root", const char *corrName8_100 = "eff_pt8_100_step_cent3accept3pt3rmin2.root", const char *corrName100_300 = "eff_pt100_300_step_cent3accept3pt3rmin2.root")
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

  TFile* corrFile0_1_p = new TFile(corrName0_1, "READ");
  cent0_1_p = (TProfile*)corrFile0_1_p->Get("p_eff_cent");
  phiEta0_1_p = (TProfile2D*)corrFile0_1_p->Get("p_eff_acceptance");
  pt0_1_p = (TProfile*)corrFile0_1_p->Get("p_eff_pt");
  delR0_1_p = (TProfile*)corrFile0_1_p->Get("p_eff_rmin");

  TFile* corrFile1_3_p = new TFile(corrName1_3, "READ");
  cent1_3_p = (TProfile*)corrFile1_3_p->Get("p_eff_cent");
  phiEta1_3_p = (TProfile2D*)corrFile1_3_p->Get("p_eff_acceptance");
  pt1_3_p = (TProfile*)corrFile1_3_p->Get("p_eff_pt");
  delR1_3_p = (TProfile*)corrFile1_3_p->Get("p_eff_rmin");

  TFile* corrFile3_8_p = new TFile(corrName3_8, "READ");
  cent3_8_p = (TProfile*)corrFile3_8_p->Get("p_eff_cent");
  phiEta3_8_p = (TProfile2D*)corrFile3_8_p->Get("p_eff_acceptance");
  pt3_8_p = (TProfile*)corrFile3_8_p->Get("p_eff_pt");
  delR3_8_p = (TProfile*)corrFile3_8_p->Get("p_eff_rmin");

  TFile* corrFile8_100_p = new TFile(corrName8_100, "READ");
  cent8_100_p = (TProfile*)corrFile8_100_p->Get("p_eff_cent");
  phiEta8_100_p = (TProfile2D*)corrFile8_100_p->Get("p_eff_acceptance");
  pt8_100_p = (TProfile*)corrFile8_100_p->Get("p_eff_pt");
  delR8_100_p = (TProfile*)corrFile8_100_p->Get("p_eff_rmin");

  TFile* corrFile100_300_p = new TFile(corrName100_300, "READ");
  cent100_300_p = (TProfile*)corrFile100_300_p->Get("p_eff_cent");
  phiEta100_300_p = (TProfile2D*)corrFile100_300_p->Get("p_eff_acceptance");
  pt100_300_p = (TProfile*)corrFile100_300_p->Get("p_eff_pt");
  delR100_300_p = (TProfile*)corrFile100_300_p->Get("p_eff_rmin");

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

  if(montecarlo)
    c->hasGenParticleTree = true;

  Long64_t nentries = c->GetEntries();

  Int_t totEv = 0;
  Int_t selectCut = 0;

  Int_t gLeadJtPtCut = 0;
  Int_t gSubLeadJtPtCut = 0;
  Int_t gDelPhiCut = 0;
  Int_t gJtEtaCut = 0;

  Int_t rPFLeadJtPtCut = 0;
  Int_t rPFSubLeadJtPtCut = 0;
  Int_t rPFDelPhiCut = 0;
  Int_t rPFJtEtaCut = 0;

  Int_t rCaloLeadJtPtCut = 0;
  Int_t rCaloSubLeadJtPtCut = 0;
  Int_t rCaloDelPhiCut = 0;
  Int_t rCaloJtEtaCut = 0;

  Int_t gTotTrk = 0;
  Int_t gTrkEtaCut = 0;
  Int_t gTrkPtCut = 0;
  Int_t gPurityCut = 0;
  Int_t gTrkDzCut = 0;
  Int_t gTrkDxyCut = 0;
  Int_t gTrkPtErrorCut = 0;

  Int_t gTotGen = 0;
  Int_t gGenEtaCut = 0;
  Int_t gGenPtCut = 0;
  Int_t gGenChgCut = 0;

  Int_t rPFTotTrk = 0;
  Int_t rPFTrkEtaCut = 0;
  Int_t rPFTrkPtCut = 0;
  Int_t rPFPurityCut = 0;
  Int_t rPFTrkDzCut = 0;
  Int_t rPFTrkDxyCut = 0;
  Int_t rPFTrkPtErrorCut = 0;

  Int_t rPFTotGen = 0;
  Int_t rPFGenEtaCut = 0;
  Int_t rPFGenPtCut = 0;
  Int_t rPFGenChgCut = 0;

  Int_t rCaloTotTrk = 0;
  Int_t rCaloTrkEtaCut = 0;
  Int_t rCaloTrkPtCut = 0;
  Int_t rCaloPurityCut = 0;
  Int_t rCaloTrkDzCut = 0;
  Int_t rCaloTrkDxyCut = 0;
  Int_t rCaloTrkPtErrorCut = 0;

  Int_t rCaloTotGen = 0;
  Int_t rCaloGenEtaCut = 0;
  Int_t rCaloGenPtCut = 0;
  Int_t rCaloGenChgCut = 0;

  defTrkCorr();

  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    c->GetEntry(jentry);

    totEv++;

    Bool_t gEventPass = false;
    Bool_t rPFEventPass = false;
    Bool_t rCaloEventPass = false;

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
    rPFLeadJtPt_ = subLeadJtPtCut;
    rPFSubLeadJtPt_ = subLeadJtPtCut;
    for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.nref; jtEntry++){
      if(c->akPu3PF.jtpt[jtEntry] > leadJtPtCut && c->akPu3PF.jtpt[jtEntry] > rPFLeadJtPt_){
	subLeadJtIndex = leadJtIndex;
	rPFSubLeadJtPt_ = rPFLeadJtPt_;
	leadJtIndex = jtEntry;
	rPFLeadJtPt_ = c->akPu3PF.jtpt[jtEntry];
      }
      else if(c->akPu3PF.jtpt[jtEntry] > rPFSubLeadJtPt_){
	subLeadJtIndex = jtEntry;
	rPFSubLeadJtPt_ = c->akPu3PF.jtpt[jtEntry];
      }
    }

    if(leadJtIndex < 0){
      rPFLeadJtPtCut++;
      rPFLeadJtPt_ = -10;
      rPFSubLeadJtPt_ = -10;
    }
    else if(subLeadJtIndex < 0){
      rPFSubLeadJtPtCut++;
      rPFSubLeadJtPt_ = -10;
    }
    else if(getAbsDphi(c->akPu3PF.jtphi[leadJtIndex], c->akPu3PF.jtphi[subLeadJtIndex]) < jtDelPhiCut){
      rPFDelPhiCut++;
    }
    else if(TMath::Abs(c->akPu3PF.jteta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akPu3PF.jteta[subLeadJtIndex]) > jtEtaCut){
      rPFJtEtaCut++;
    }
    else{
      rPFLeadJtPhi_ = c->akPu3PF.jtphi[leadJtIndex];
      rPFSubLeadJtPhi_ = c->akPu3PF.jtphi[subLeadJtIndex];
      rPFLeadJtEta_ = c->akPu3PF.jteta[leadJtIndex];
      rPFSubLeadJtEta_ = c->akPu3PF.jteta[subLeadJtIndex];

      rPFEventPass = true;

      recoPFSet_ = true;
    }

    //calo jet

    leadJtIndex = -1;
    subLeadJtIndex = -1;
    rCaloLeadJtPt_ = subLeadJtPtCut;
    rCaloSubLeadJtPt_ = subLeadJtPtCut;
    for(Int_t jtEntry = 0; jtEntry < c->akPu3Calo.nref; jtEntry++){
      if(c->akPu3Calo.jtpt[jtEntry] > leadJtPtCut && c->akPu3Calo.jtpt[jtEntry] > rCaloLeadJtPt_){
	subLeadJtIndex = leadJtIndex;
	rCaloSubLeadJtPt_ = rCaloLeadJtPt_;
	leadJtIndex = jtEntry;
	rCaloLeadJtPt_ = c->akPu3Calo.jtpt[jtEntry];
      }
      else if(c->akPu3Calo.jtpt[jtEntry] > rCaloSubLeadJtPt_){
	subLeadJtIndex = jtEntry;
	rCaloSubLeadJtPt_ = c->akPu3Calo.jtpt[jtEntry];
      }
    }

    if(leadJtIndex < 0){
      rCaloLeadJtPtCut++;
      rCaloLeadJtPt_ = -10;
      rCaloSubLeadJtPt_ = -10;
    }
    else if(subLeadJtIndex < 0){
      rCaloSubLeadJtPtCut++;
      rCaloSubLeadJtPt_ = -10;
    }
    else if(getAbsDphi(c->akPu3Calo.jtphi[leadJtIndex], c->akPu3Calo.jtphi[subLeadJtIndex]) < jtDelPhiCut){
      rCaloDelPhiCut++;
    }
    else if(TMath::Abs(c->akPu3Calo.jteta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akPu3Calo.jteta[subLeadJtIndex]) > jtEtaCut){
      rCaloJtEtaCut++;
    }
    else{
      rCaloLeadJtPhi_ = c->akPu3Calo.jtphi[leadJtIndex];
      rCaloSubLeadJtPhi_ = c->akPu3Calo.jtphi[subLeadJtIndex];
      rCaloLeadJtEta_ = c->akPu3Calo.jteta[leadJtIndex];
      rCaloSubLeadJtEta_ = c->akPu3Calo.jteta[subLeadJtIndex];

      rCaloEventPass = true;

      recoCaloSet_ = true;
    }

    //truth

    leadJtIndex = -1;
    subLeadJtIndex = -1;
    gLeadJtPt_ = subLeadJtPtCut;
    gSubLeadJtPt_ = subLeadJtPtCut;
    for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.nref; jtEntry++){
      if(c->akPu3PF.refpt[jtEntry] > leadJtPtCut && c->akPu3PF.refpt[jtEntry] > gLeadJtPt_){
	subLeadJtIndex = leadJtIndex;
	gSubLeadJtPt_ = gLeadJtPt_;
	leadJtIndex = jtEntry;
	gLeadJtPt_ = c->akPu3PF.refpt[jtEntry];
      }
      else if(c->akPu3PF.refpt[jtEntry] > gSubLeadJtPt_){
	subLeadJtIndex = jtEntry;
	gSubLeadJtPt_ = c->akPu3PF.refpt[jtEntry];
      }
    }

    if(leadJtIndex < 0){
      gLeadJtPtCut++;
      gLeadJtPt_ = -10;
      gSubLeadJtPt_ = -10;
    }
    else if(subLeadJtIndex < 0){
      gSubLeadJtPtCut++;
      gSubLeadJtPt_ = -10;
    }
    else if(getAbsDphi(c->akPu3PF.refphi[leadJtIndex], c->akPu3PF.refphi[subLeadJtIndex]) < jtDelPhiCut){
      gDelPhiCut++;
    }
    else if(TMath::Abs(c->akPu3PF.refeta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akPu3PF.refeta[subLeadJtIndex]) > jtEtaCut){
      gJtEtaCut++;
    }
    else{
      gLeadJtPhi_ = c->akPu3PF.refphi[leadJtIndex];
      gSubLeadJtPhi_ = c->akPu3PF.refphi[subLeadJtIndex];
      gLeadJtEta_ = c->akPu3PF.refeta[leadJtIndex];
      gSubLeadJtEta_ = c->akPu3PF.refeta[subLeadJtIndex];

      gRPFLeadJtPt_ = c->akPu3PF.jtpt[leadJtIndex];
      gRPFSubLeadJtPt_ = c->akPu3PF.jtpt[subLeadJtIndex];
      gRPFLeadJtPhi_ = c->akPu3PF.jtphi[leadJtIndex];
      gRPFSubLeadJtPhi_ = c->akPu3PF.jtphi[subLeadJtIndex];
      gRPFLeadJtEta_ = c->akPu3PF.jteta[leadJtIndex];
      gRPFSubLeadJtEta_ = c->akPu3PF.jteta[subLeadJtIndex];

      for(Int_t jtEntry = 0; jtEntry < c->akPu3Calo.nref; jtEntry++){
	if(getDR(c->akPu3Calo.jteta[jtEntry], gLeadJtEta_, c->akPu3Calo.jtphi[jtEntry], gLeadJtPhi_) < 0.3){
	  gRCaloLeadJtPt_ = c->akPu3Calo.jtpt[jtEntry];
	  gRCaloLeadJtPhi_ = c->akPu3Calo.jtphi[jtEntry];
	  gRCaloLeadJtEta_ = c->akPu3Calo.jteta[jtEntry];
	}

	if(getDR(c->akPu3Calo.jteta[jtEntry], gSubLeadJtEta_, c->akPu3Calo.jtphi[jtEntry], gSubLeadJtPhi_) < 0.3){
	  gRCaloSubLeadJtPt_ = c->akPu3Calo.jtpt[jtEntry];
	  gRCaloSubLeadJtPhi_ = c->akPu3Calo.jtphi[jtEntry];
	  gRCaloSubLeadJtEta_ = c->akPu3Calo.jteta[jtEntry];
	}
      }

      gEventPass = true;

      truthSet_ = true;
    }

    if(gEventPass == false && rPFEventPass == false && rCaloEventPass == false)
      continue;
    

    run_ = c->evt.run;
    evt_ = c->akPu3PF.evt;
    lumi_ = c->evt.lumi;
    hiBin_ = c->evt.hiBin;

    //Iterate over tracks

    nTrk_ = 0;

    InitProjPerp(montecarlo);

    if(montecarlo){
      for(Int_t divIter = 0; divIter < 10; divIter++){
	rDivGPt_[divIter] = 0;
      }
    }

    Tracks trkCollection;
    trkCollection = c->track;

    for(Int_t trkEntry = 0; trkEntry < trkCollection.nTrk; trkEntry++){
      if(gEventPass)	gTotTrk++;
      if(rPFEventPass)	rPFTotTrk++;
      if(rCaloEventPass)	rCaloTotTrk++;
      
      if(TMath::Abs(trkCollection.trkEta[trkEntry]) > 2.4){
	if(gEventPass)	  gTrkEtaCut++;
	if(rPFEventPass)	  rPFTrkEtaCut++;
	if(rCaloEventPass)	  rCaloTrkEtaCut++;

	continue;
      }
      
      if(trkCollection.trkPt[trkEntry] < 0.5){
	if(gEventPass)	  gTrkPtCut++;
	if(rPFEventPass)	  rPFTrkPtCut++;
	if(rCaloEventPass)	  rCaloTrkPtCut++;

	continue;
      }
      
      if(!trkCollection.highPurity[trkEntry]){ //Note highPuritySetWithPV seems to be wrong cut, creates diff. bet. truth and reco
	if(gEventPass)	  gPurityCut++;
	if(rPFEventPass)	  rPFPurityCut++;
	if(rCaloEventPass)	  rCaloPurityCut++;

	continue;
      }

      if(TMath::Abs(trkCollection.trkDz1[trkEntry]/trkCollection.trkDzError1[trkEntry]) > 3){
	if(gEventPass)	  gTrkDzCut++;
	if(rPFEventPass)	  rPFTrkDzCut++;
	if(rCaloEventPass)	  rCaloTrkDzCut++;

	continue;
      }

      if(TMath::Abs(trkCollection.trkDxy1[trkEntry]/trkCollection.trkDxyError1[trkEntry]) > 3){
	if(gEventPass)	  gTrkDxyCut++;
	if(rPFEventPass)	  rPFTrkDxyCut++;
	if(rCaloEventPass)	  rCaloTrkDxyCut++;

	continue;
      }

      if(trkCollection.trkPtError[trkEntry]/trkCollection.trkPt[trkEntry] > 0.1){
	if(gEventPass)	  gTrkPtErrorCut++;
	if(rPFEventPass)	  rPFTrkPtErrorCut++;
	if(rCaloEventPass)	  rCaloTrkPtErrorCut++;

	continue;
      }

      trkPt_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPhi_[nTrk_] = trkCollection.trkPhi[trkEntry];
      trkEta_[nTrk_] = trkCollection.trkEta[trkEntry];
      
      if(rPFEventPass)
	trkPtRPF_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));
      else
	trkPtRPF_[nTrk_] = 0;

      if(montecarlo && gEventPass)
	trkPtGRPF_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
      else
	trkPtGRPF_[nTrk_] = 0;

      if(rPFEventPass)
	trkPFLeadDelPhi_[nTrk_] = getAbsDphi(rPFLeadJtPhi_, trkCollection.trkPhi[trkEntry]);
      else
	trkPFLeadDelPhi_[nTrk_] = -10;

      if(rCaloEventPass)
	trkCaloLeadDelPhi_[nTrk_] = getAbsDphi(rCaloLeadJtPhi_, trkCollection.trkPhi[trkEntry]);
      else
	trkCaloLeadDelPhi_[nTrk_] = -10;


      if(montecarlo){
	for(Int_t divIter = 0; divIter < 10; divIter++){
	  if(trkPt_[nTrk_] > 20.)
	    break;
	  else if(2*divIter < trkPt_[nTrk_] && 2*(divIter + 1.) > trkPt_[nTrk_]){
	    rDivGPt_[divIter]++;
	    break;
	  }
	}
      }
      
      
      if(rPFEventPass){
	rPFImbProjF_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));
	rPFImbPerpF_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));
	if(trkCollection.trkPt[trkEntry] > 8){
	  rPFImbProjH_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));
	  rPFImbPerpH_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));

	  rPFImbProj8_100_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));
	}
	else{
	  rPFImbProjL_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));
	  rPFImbPerpL_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));

	  if(trkCollection.trkPt[trkEntry] < 1){
	    rPFImbProj5_1_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 2){
	    rPFImbProj1_2_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 4){
	    rPFImbProj2_4_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 8){
	    rPFImbProj4_8_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rPFLeadJtPhi_));
	  }

	}
      }

      if(montecarlo && gEventPass){
	gRPFImbProjF_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	gRPFImbPerpF_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	if(trkCollection.trkPt[trkEntry] > 8){
	  gRPFImbProjH_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	  gRPFImbPerpH_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));

	  gRPFImbProj8_100_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	}
	else{
	  gRPFImbProjL_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	  gRPFImbPerpL_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));

	  if(trkCollection.trkPt[trkEntry] < 1){
	    gRPFImbProj5_1_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 2){
	    gRPFImbProj1_2_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 4){
	    gRPFImbProj2_4_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	  }
	  else if(trkCollection.trkPt[trkEntry] < 8){
	    gRPFImbProj4_8_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	  }

	}
      }


      if(rCaloEventPass){
	rCaloImbProjF_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rCaloLeadJtPhi_));
	rCaloImbPerpF_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], rCaloLeadJtPhi_));
	if(trkCollection.trkPt[trkEntry] > 8){
	  rCaloImbProjH_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rCaloLeadJtPhi_));
	  rCaloImbPerpH_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], rCaloLeadJtPhi_));
	}
	else{
	  rCaloImbProjL_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rCaloLeadJtPhi_));
	  rCaloImbPerpL_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], rCaloLeadJtPhi_));
	}
      }

      nTrk_++;
      if(nTrk_ > MAXTRKS - 1){
	printf("ERROR: Trk arrays not large enough.\n");
	return(1);
      }
    }



    for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
      trkRMin_[trkEntry] = 6;

      for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.nref; jtEntry++){
	if(c->akPu3PF.jtpt[jtEntry] < 30 || TMath::Abs(c->akPu3PF.jteta[jtEntry]) > 2.0)
	  continue;

	if(trkRMin_[trkEntry] > getDR(trkEta_[trkEntry], trkPhi_[trkEntry], c->akPu3PF.jteta[jtEntry], c->akPu3PF.jtphi[jtEntry]))
	  trkRMin_[trkEntry] = getDR(trkEta_[trkEntry], trkPhi_[trkEntry], c->akPu3PF.jteta[jtEntry], c->akPu3PF.jtphi[jtEntry]);
      }

      if(trkPt_[trkEntry] > 0.5 && trkPt_[trkEntry] < 1.0){
	trkPtCorr_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMin_[trkEntry], cent0_1_p, phiEta0_1_p, pt0_1_p, delR0_1_p);
	trkPtFact_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMin_[trkEntry], cent0_1_p, phiEta0_1_p, pt0_1_p, delR0_1_p);
      }
      else if(trkPt_[trkEntry] > 1 && trkPt_[trkEntry] < 3){
	trkPtCorr_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMin_[trkEntry], cent1_3_p, phiEta1_3_p, pt1_3_p, delR1_3_p);
	trkPtFact_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMin_[trkEntry], cent1_3_p, phiEta1_3_p, pt1_3_p, delR1_3_p);
      }
      else if(trkPt_[trkEntry] > 3 && trkPt_[trkEntry] < 8){
	trkPtCorr_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMin_[trkEntry], cent3_8_p, phiEta3_8_p, pt3_8_p, delR3_8_p);
	trkPtFact_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMin_[trkEntry], cent3_8_p, phiEta3_8_p, pt3_8_p, delR3_8_p);
      }
      else if(trkPt_[trkEntry] > 8 && trkPt_[trkEntry] < 100){
	trkPtCorr_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMin_[trkEntry], cent8_100_p, phiEta8_100_p, pt8_100_p, delR8_100_p);
	trkPtFact_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMin_[trkEntry], cent8_100_p, phiEta8_100_p, pt8_100_p, delR8_100_p);
      }
      else if(trkPt_[trkEntry] > 100 && trkPt_[trkEntry] < 300){
	trkPtCorr_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMin_[trkEntry], cent100_300_p, phiEta100_300_p, pt100_300_p, delR100_300_p);
	trkPtFact_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMin_[trkEntry], cent100_300_p, phiEta100_300_p, pt100_300_p, delR100_300_p);
      }
      else{
	trkPtCorr_[trkEntry] = trkPt_[trkEntry];
	trkPtFact_[trkEntry] = 1;
      }
    }

    if(rPFEventPass){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	rPFImbProjFCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], rPFLeadJtPhi_));             
        rPFImbPerpFCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], rPFLeadJtPhi_));             
        if(trkPt_[trkEntry] > 8){                                                                     
          rPFImbProjHCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], rPFLeadJtPhi_));           
          rPFImbPerpHCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], rPFLeadJtPhi_));           

          rPFImbProjCorr8_100_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], rPFLeadJtPhi_));           
        }                                                                                                 
        else{                                                                                             
          rPFImbProjLCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], rPFLeadJtPhi_));           
          rPFImbPerpLCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], rPFLeadJtPhi_));           

	  if(trkPt_[trkEntry] < 1){
	    rPFImbProjCorr5_1_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], rPFLeadJtPhi_));
	  }
	  else if(trkPt_[trkEntry] < 2){
	    rPFImbProjCorr1_2_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], rPFLeadJtPhi_));
	  }
	  else if(trkPt_[trkEntry] < 4){
	    rPFImbProjCorr2_4_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], rPFLeadJtPhi_));
	  }
	  else if(trkPt_[trkEntry] < 8){
	    rPFImbProjCorr4_8_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], rPFLeadJtPhi_));
	  }

        }                                                                                                 
      }
    }

    if(montecarlo && gEventPass){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
        gRPFImbProjFCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));            
        gRPFImbPerpFCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));            
        if(trkPt_[trkEntry] > 8){                                                                     
          gRPFImbProjHCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));          
          gRPFImbPerpHCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));          

          gRPFImbProjCorr8_100_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));          
        }                                                                                                 
        else{                                                                                             
          gRPFImbProjLCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));          
          gRPFImbPerpLCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));          

	  if(trkPt_[trkEntry] < 1){
	    gRPFImbProjCorr5_1_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));          
	  }
	  else if(trkPt_[trkEntry] < 2){
	    gRPFImbProjCorr1_2_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));          
	  }
	  else if(trkPt_[trkEntry] < 4){
	    gRPFImbProjCorr2_4_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));          
	  }
	  else if(trkPt_[trkEntry] < 8){
	    gRPFImbProjCorr4_8_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));          
	  }

        }                                                                                                 
      }
    }


    Float_t rDivGPt_temp[10];

    if(montecarlo){
      //Iterate over truth

      nGen_ = 0;

      for(Int_t divIter = 0; divIter < 10; divIter++){
	rDivGPt_temp[divIter] = 0;
      }

      GenParticles genCollection;
      genCollection = c->genparticle;

      for(Int_t genEntry = 0; genEntry < genCollection.mult; genEntry++){
	if(gEventPass)	  gTotGen++;
	if(rPFEventPass)	  rPFTotGen++;
	if(rCaloEventPass)	  rCaloTotGen++;
	  
	if(genCollection.chg[genEntry] == 0){
	  if(gEventPass)	    gGenChgCut++;
	  if(rPFEventPass)	    rPFGenChgCut++;
	  if(rCaloEventPass)	    rCaloGenChgCut++;
	  
	  continue;
	}
	  
	if(TMath::Abs(genCollection.eta[genEntry]) > 2.4){
	  if(gEventPass)	    gGenEtaCut++;
	  if(rPFEventPass)	    rPFGenEtaCut++;
	  if(rCaloEventPass)	    rCaloGenEtaCut++;

	  continue;
	}
	
	if(genCollection.pt[genEntry] < 0.5){
	  if(gEventPass)	    gGenPtCut++;
	  if(rPFEventPass)	    rPFGenPtCut++;
	  if(rCaloEventPass)	    rCaloGenPtCut++;

	  continue;
	}

	genPt_[nGen_] = genCollection.pt[genEntry];
	genPhi_[nGen_] = genCollection.phi[genEntry];
	genEta_[nGen_] = genCollection.eta[genEntry];

	if(gEventPass)
          genPtJT_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	else
	  genPtJT_[nGen_] = 0;
	
	if(gEventPass)
	  genLeadDelPhi_[nGen_] = getAbsDphi(gLeadJtPhi_, genCollection.phi[genEntry]);
	else
	  genLeadDelPhi_[nGen_] = -10;

	for(Int_t divIter = 0; divIter < 10; divIter++){
	  if(genPt_[nGen_] > 20.)
	    break;
	  else if(2*divIter < genPt_[nGen_] && 2*(divIter + 1.) > genPt_[nGen_]){
	    rDivGPt_temp[divIter]++;
	    break;
	  }
	}

	if(gEventPass){
	  gImbProjF_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	  gImbPerpF_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	  if(genCollection.pt[genEntry] > 8){
	    gImbProjH_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	    gImbPerpH_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	    gImbProj8_100_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	  }
	  else{
	    gImbProjL_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	    gImbPerpL_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));

	    if(genCollection.pt[genEntry] < 1){
	      gImbProj5_1_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	    }
	    else if(genCollection.pt[genEntry] < 2){
	      gImbProj1_2_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	    }
	    else if(genCollection.pt[genEntry] < 4){
	      gImbProj2_4_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	    }
	    else if(genCollection.pt[genEntry] < 8){
	      gImbProj4_8_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
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

    if(montecarlo){
      for(Int_t divIter = 0; divIter < 10; divIter++){
	if(rDivGPt_temp[divIter] != 0)
	  rDivGPt_[divIter] = rDivGPt_[divIter]/rDivGPt_temp[divIter];
	else 
	  rDivGPt_[divIter] = -1;
      }
    }

    jetTree_p->Fill();
    trackTree_p->Fill();

    if(montecarlo)
      genTree_p->Fill();
  }

  std::cout << "totEv: " << totEv << std::endl;
  Int_t tempTot = totEv - selectCut;
  std::cout << "selectCut: " << tempTot << std::endl;

  std::cout << std::endl;
  tempTot = tempTot - gLeadJtPtCut;
  std::cout << "gLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - gSubLeadJtPtCut;
  std::cout << "gSubLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - gDelPhiCut;
  std::cout << "gDelPhiCut: " << tempTot << std::endl;
  tempTot = tempTot - gJtEtaCut;
  std::cout << "gJtEtaCut: " << tempTot << std::endl;

  std::cout << std::endl;
  tempTot = totEv - selectCut - rPFLeadJtPtCut;
  std::cout << "rPFLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - rPFSubLeadJtPtCut;
  std::cout << "rPFSubLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - rPFDelPhiCut;
  std::cout << "rPFDelPhiCut: " << tempTot << std::endl;
  tempTot = tempTot - rPFJtEtaCut;
  std::cout << "rPFJtEtaCut: " << tempTot << std::endl;

  std::cout << std::endl;
  tempTot = totEv - selectCut - rCaloLeadJtPtCut;
  std::cout << "rCaloLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - rCaloSubLeadJtPtCut;
  std::cout << "rCaloSubLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - rCaloDelPhiCut;
  std::cout << "rCaloDelPhiCut: " << tempTot << std::endl;
  tempTot = tempTot - rCaloJtEtaCut;
  std::cout << "rCaloJtEtaCut: " << tempTot << std::endl;

  std::cout << std::endl;
  std::cout << "gTotTrk: " << gTotTrk << std::endl;
  tempTot = gTotTrk - gTrkEtaCut;
  std::cout << "gTrkEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - gTrkPtCut;
  std::cout << "gTrkPtCut: " << tempTot << std::endl;
  tempTot = tempTot - gPurityCut;
  std::cout << "gPurityCut: " << tempTot << std::endl;
  tempTot = tempTot - gTrkDzCut;
  std::cout << "gTrkDzCut: " << tempTot << std::endl;
  tempTot = tempTot - gTrkDxyCut;
  std::cout << "gTrkDxyCut: " << tempTot << std::endl;
  tempTot = tempTot - gTrkPtErrorCut;
  std::cout << "gTrkPtErrorCut: " << tempTot << std::endl;

  std::cout << std::endl;
  std::cout << "gTotGen: " << gTotGen << std::endl;
  tempTot = gTotGen - gGenChgCut;
  std::cout << "gGenChgCut: " << tempTot << std::endl;
  tempTot = tempTot - gGenEtaCut;
  std::cout << "gGenEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - gGenPtCut;
  std::cout << "gGenPtCut: " << tempTot << std::endl;
  


  std::cout << std::endl;
  std::cout << "rPFTotTrk: " << rPFTotTrk << std::endl;
  tempTot = rPFTotTrk - rPFTrkEtaCut;
  std::cout << "rPFTrkEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - rPFTrkPtCut;
  std::cout << "rPFTrkPtCut: " << tempTot << std::endl;
  tempTot = tempTot - rPFPurityCut;
  std::cout << "rPFPurityCut: " << tempTot << std::endl;
  tempTot = tempTot - rPFTrkDzCut;
  std::cout << "rPFTrkDzCut: " << tempTot << std::endl;
  tempTot = tempTot - rPFTrkDxyCut;
  std::cout << "rPFTrkDxyCut: " << tempTot << std::endl;
  tempTot = tempTot - rPFTrkPtErrorCut;
  std::cout << "rPFTrkPtErrorCut: " << tempTot << std::endl;

  std::cout << std::endl;
  std::cout << "rPFTotGen: " << rPFTotGen << std::endl;
  tempTot = rPFTotGen - rPFGenChgCut;
  std::cout << "rPFGenChgCut: " << tempTot << std::endl;
  tempTot = tempTot - rPFGenEtaCut;
  std::cout << "rPFGenEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - rPFGenPtCut;
  std::cout << "rPFGenPtCut: " << tempTot << std::endl;


  std::cout << std::endl;
  std::cout << "rCaloTotTrk: " << rCaloTotTrk << std::endl;
  tempTot = rCaloTotTrk - rCaloTrkEtaCut;
  std::cout << "rCaloTrkEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - rCaloTrkPtCut;
  std::cout << "rCaloTrkPtCut: " << tempTot << std::endl;
  tempTot = tempTot - rCaloPurityCut;
  std::cout << "rCaloPurityCut: " << tempTot << std::endl;
  tempTot = tempTot - rCaloTrkDzCut;
  std::cout << "rCaloTrkDzCut: " << tempTot << std::endl;
  tempTot = tempTot - rCaloTrkDxyCut;
  std::cout << "rCaloTrkDxyCut: " << tempTot << std::endl;
  tempTot = tempTot - rCaloTrkPtErrorCut;
  std::cout << "rCaloTrkPtErrorCut: " << tempTot << std::endl;

  std::cout << std::endl;
  std::cout << "rCaloTotGen: " << rCaloTotGen << std::endl;
  tempTot = rCaloTotGen - rCaloGenChgCut;
  std::cout << "rCaloGenChgCut: " << tempTot << std::endl;
  tempTot = tempTot - rCaloGenEtaCut;
  std::cout << "rCaloGenEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - rCaloGenPtCut;
  std::cout << "rCaloGenPtCut: " << tempTot << std::endl;

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
