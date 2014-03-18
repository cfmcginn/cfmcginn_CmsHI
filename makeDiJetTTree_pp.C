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
#include "cfmDiJetSkim_pp.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include "factorizedPtCorr_pp.h"

const Float_t leadJtPtCut = 120.;
const Float_t subLeadJtPtCut = 50.;
const Float_t jtDelPhiCut = 2.*(TMath::Pi())/3.;
const Float_t jtEtaCut = 1.6; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed

collisionType getCType(sampleType sType);


void getLeadJt(Float_t& leadJtPt, Int_t& lPtCut, Float_t& subLeadJtPt, Int_t& sPtCut, Float_t& leadJtPhi, Float_t& subLeadJtPhi, Int_t& dPhiCut, Float_t& leadJtEta, Float_t& subLeadJtEta, Int_t& etaCut, Float_t& jtDelPhi, Float_t& jtAsymm, Bool_t& setPass, Float_t& refLPt, Float_t& refSLPt, Float_t& refLEta, Float_t& refSLEta, Jets jtCollection, Bool_t& jtVeto, Bool_t montecarlo = false)
{
  Int_t leadJtIndex = -1;
  Int_t subLeadJtIndex = -1;
  leadJtPt = subLeadJtPtCut;
  subLeadJtPt = subLeadJtPtCut;

  for(Int_t jtEntry = 0; jtEntry < jtCollection.nref; jtEntry++){
    if(jtCollection.jtpt[jtEntry] > leadJtPtCut && jtCollection.jtpt[jtEntry] > leadJtPt && TMath::Abs(jtCollection.jteta[jtEntry]) < jtEtaCut){
      subLeadJtIndex = leadJtIndex;
      subLeadJtPt = leadJtPt;
      leadJtIndex = jtEntry;
      leadJtPt = jtCollection.jtpt[jtEntry];
    }
    else if(jtCollection.jtpt[jtEntry] > subLeadJtPt && TMath::Abs(jtCollection.jteta[jtEntry]) < jtEtaCut){
      subLeadJtIndex = jtEntry;
      subLeadJtPt = jtCollection.jtpt[jtEntry];
    }
    else if(jtCollection.jtpt[jtEntry] > subLeadJtPtCut && TMath::Abs(jtCollection.jteta[jtEntry]) < jtEtaCut){
      jtVeto = true;
    }
  }

  if(leadJtIndex < 0){
    lPtCut++; 
    leadJtPt = -10;
    subLeadJtPt = -10;
  }
  else if(subLeadJtIndex < 0){
    sPtCut++;
    subLeadJtPt = -10;
  }
  else if(getAbsDphi(jtCollection.jtphi[leadJtIndex], jtCollection.jtphi[subLeadJtIndex]) < jtDelPhiCut){
    dPhiCut++;
  }
  else if(TMath::Abs(jtCollection.jteta[leadJtIndex]) > jtEtaCut || TMath::Abs(jtCollection.jteta[subLeadJtIndex]) > jtEtaCut){
    etaCut++;
  }
  else{
    leadJtPhi = jtCollection.jtphi[leadJtIndex];
    subLeadJtPhi = jtCollection.jtphi[subLeadJtIndex];
    leadJtEta = jtCollection.jteta[leadJtIndex];
    subLeadJtEta = jtCollection.jteta[subLeadJtIndex];

    jtDelPhi = getAbsDphi(jtCollection.jtphi[leadJtIndex], jtCollection.jtphi[subLeadJtIndex]);
    jtAsymm = (leadJtPt - subLeadJtPt)/(leadJtPt + subLeadJtPt);

    if(montecarlo){
      refLPt = jtCollection.refpt[leadJtIndex];
      refLEta = jtCollection.refeta[leadJtIndex];
      refSLPt = jtCollection.refpt[subLeadJtIndex];
      refSLEta = jtCollection.refeta[subLeadJtIndex];
    }

    setPass = true;
  }

  return;
}


void getPtProj(Float_t cutPt, Float_t inPt, Float_t phi, Float_t jtPhi, Float_t& ProjF, Float_t& PerpF, Float_t& Proj0_1, Float_t& Proj1_2, Float_t& Proj2_4, Float_t& Proj4_8, Float_t& Proj8_100)
{
  ProjF += -inPt*cos(getDPHI(phi, jtPhi));
  PerpF += -inPt*sin(getDPHI(phi, jtPhi));

  if(cutPt < 1)
    Proj0_1 += -inPt*cos(getDPHI(phi, jtPhi));
  else if(cutPt < 2)
    Proj1_2 += -inPt*cos(getDPHI(phi, jtPhi));
  else if(cutPt < 4)
    Proj2_4 += -inPt*cos(getDPHI(phi, jtPhi));
  else if(cutPt < 8)
    Proj4_8 += -inPt*cos(getDPHI(phi, jtPhi));
  else
    Proj8_100 += -inPt*cos(getDPHI(phi, jtPhi));

  return;
}


void getPtHem(Float_t cutPt, Float_t inPt, Float_t phi, Float_t jtPhi, Float_t& HemF, Float_t& Hem0_1, Float_t& Hem1_2, Float_t& Hem2_4, Float_t& Hem4_8, Float_t& Hem8_100)
{
  Float_t modInPt = inPt;
  if(cos(getDPHI(phi, jtPhi)) < 0)
    modInPt = -inPt;

  HemF += -modInPt;

  if(cutPt < 1)
    Hem0_1 += -modInPt;
  else if(cutPt < 2)
    Hem1_2 += -modInPt;
  else if(cutPt < 4)
    Hem2_4 += -modInPt;
  else if(cutPt < 8)
    Hem4_8 += -modInPt;
  else
    Hem8_100 += -modInPt;

  return;
}


Float_t getTrkRMin(Float_t phi, Float_t eta, Jets jtCollection, Bool_t isGen = false)
{
  Float_t trkRMin = 10;

  if(!isGen){
    for(Int_t jtEntry = 0; jtEntry < jtCollection.nref; jtEntry++){
      if(jtCollection.jtpt[jtEntry] < 30 || TMath::Abs(jtCollection.jteta[jtEntry]) > 2.0)
	continue;

      if(trkRMin > getDR(eta, phi, jtCollection.jteta[jtEntry], jtCollection.jtphi[jtEntry]))
	trkRMin = getDR(eta, phi, jtCollection.jteta[jtEntry], jtCollection.jtphi[jtEntry]);
    }
  }
  else if(isGen){
    for(Int_t jtEntry = 0; jtEntry < jtCollection.ngen; jtEntry++){
      if(jtCollection.genpt[jtEntry] < 30 || TMath::Abs(jtCollection.geneta[jtEntry]) > 2.0)
        continue;

      if(trkRMin > getDR(eta, phi, jtCollection.geneta[jtEntry], jtCollection.genphi[jtEntry]))
        trkRMin = getDR(eta, phi, jtCollection.geneta[jtEntry], jtCollection.genphi[jtEntry]);
    }
  }

  return trkRMin;
}


int makeDiJetTTree_pp(string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_CFMSKIM.root", Int_t num = 0)
{
  //Define MC or Data
  bool montecarlo = false;
  if(sType == kPPMC || sType == kPAMC || sType == kHIMC)
    montecarlo = true;

  std::cout << sType << std::endl;
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

  //Setup correction tables

  InitCorrFiles();
  InitCorrHists();


  TFile *outFile = new TFile(Form("%s_%d.root", outName, num), "RECREATE");

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

  c->hasAk3JetTree = true;
  c->hasAk3CaloJetTree = true;

  if(montecarlo)
    c->hasGenParticleTree = true;


  Long64_t nentries = c->GetEntries();

  std::cout << nentries << std::endl;

  Int_t totEv = 0;
  Int_t selectCut = 0;

  Int_t AlgLeadJtPtCut[5] = {0, 0, 0, 0, 0};
  Int_t AlgSubLeadJtPtCut[5] = {0, 0, 0, 0, 0};
  Int_t AlgDelPhiCut[5] = {0, 0, 0, 0, 0};
  Int_t AlgJtEtaCut[5] = {0, 0, 0, 0, 0};

  Int_t AlgTotTrk[5] = {0, 0, 0, 0, 0};
  Int_t AlgTrkEtaCut[5] = {0, 0, 0, 0, 0};
  Int_t AlgTrkPtCut[5] = {0, 0, 0, 0, 0};
  Int_t AlgPurityCut[5] = {0, 0, 0, 0, 0};
  Int_t AlgTrkDzCut[5] = {0, 0, 0, 0, 0};
  Int_t AlgTrkDxyCut[5] = {0, 0, 0, 0, 0};
  Int_t AlgTrkPtErrorCut[5] = {0, 0, 0, 0, 0};

  Int_t AlgTotGen[5] = {0, 0, 0, 0, 0};
  Int_t AlgGenEtaCut[5] = {0, 0, 0, 0, 0};
  Int_t AlgGenPtCut[5] = {0, 0, 0, 0, 0};
  Int_t AlgGenChgCut[5] = {0, 0, 0, 0, 0};


  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    c->GetEntry(jentry);

    totEv++;

    if(jentry%1000 == 0)
      std::cout << jentry << std::endl;

    if(!c->selectEvent()){
      selectCut++;
      continue;
    }

    InitJetVar(montecarlo);

    //particle flow

    Jets AlgJtCollection[4] = {c->akPu3PF, c->akPu3Calo, c->ak3PF, c->ak3Calo};
    
    for(Int_t collIter = 0; collIter < 4; collIter++){
      getLeadJt(AlgLeadJtPt_[collIter], AlgLeadJtPtCut[collIter], AlgSubLeadJtPt_[collIter], AlgSubLeadJtPtCut[collIter], AlgLeadJtPhi_[collIter], AlgSubLeadJtPhi_[collIter], AlgDelPhiCut[collIter], AlgLeadJtEta_[collIter], AlgSubLeadJtEta_[collIter], AlgJtEtaCut[collIter], AlgJtDelPhi_[collIter], AlgJtAsymm_[collIter], eventSet_[collIter], AlgLeadRefPt_[collIter], AlgSubLeadRefPt_[collIter], AlgLeadRefEta_[collIter], AlgSubLeadRefEta_[collIter], AlgJtCollection[collIter], thirdJtVeto_[collIter], montecarlo);
    }

    //truth, doesn't work w/ getLeadJt because truth doesnt get its own tree

    if(montecarlo){
      Int_t leadJtIndex = -1;
      Int_t subLeadJtIndex = -1;
      AlgLeadJtPt_[T] = subLeadJtPtCut;
      AlgSubLeadJtPt_[T] = subLeadJtPtCut;
      for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.ngen; jtEntry++){
	if(c->akPu3PF.genpt[jtEntry] > leadJtPtCut && c->akPu3PF.genpt[jtEntry] > AlgLeadJtPt_[T] && TMath::Abs(c->akPu3PF.geneta[jtEntry]) < 1.6){
	  subLeadJtIndex = leadJtIndex;
	  AlgSubLeadJtPt_[T] = AlgLeadJtPt_[T];
	  leadJtIndex = jtEntry;
	  AlgLeadJtPt_[T] = c->akPu3PF.genpt[jtEntry];
	}
	else if(c->akPu3PF.genpt[jtEntry] > AlgSubLeadJtPt_[T] && TMath::Abs(c->akPu3PF.geneta[jtEntry]) < 1.6){
	  subLeadJtIndex = jtEntry;
	  AlgSubLeadJtPt_[T] = c->akPu3PF.genpt[jtEntry];
	}
	else if(c->akPu3PF.genpt[jtEntry] > subLeadJtPtCut && TMath::Abs(c->akPu3PF.geneta[jtEntry]) < 1.6){
	  thirdJtVeto_[T] = true;
	} 
      }

      if(leadJtIndex < 0){
	AlgLeadJtPtCut[T]++;
	AlgLeadJtPt_[T] = -10;
	AlgSubLeadJtPt_[T] = -10;
      }
      else if(subLeadJtIndex < 0){
	AlgSubLeadJtPtCut[T]++;
	AlgSubLeadJtPt_[T] = -10;
      }
      else if(getAbsDphi(c->akPu3PF.genphi[leadJtIndex], c->akPu3PF.genphi[subLeadJtIndex]) < jtDelPhiCut){
	AlgDelPhiCut[T]++;
      }
      else if(TMath::Abs(c->akPu3PF.geneta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akPu3PF.geneta[subLeadJtIndex]) > jtEtaCut){
	AlgJtEtaCut[T]++;
      }
      else{
	AlgLeadJtPhi_[T] = c->akPu3PF.genphi[leadJtIndex];
	AlgSubLeadJtPhi_[T] = c->akPu3PF.genphi[subLeadJtIndex];
	AlgLeadJtEta_[T] = c->akPu3PF.geneta[leadJtIndex];
	AlgSubLeadJtEta_[T] = c->akPu3PF.geneta[subLeadJtIndex];
	
	AlgJtDelPhi_[T] = getAbsDphi(c->akPu3PF.genphi[leadJtIndex], c->akPu3PF.genphi[subLeadJtIndex]);
	AlgJtAsymm_[T] = (AlgLeadJtPt_[T] - AlgSubLeadJtPt_[T])/(AlgLeadJtPt_[T] + AlgSubLeadJtPt_[T]);
	
	eventSet_[T] = true;
      }
    }
    
    if(eventSet_[PuPF] == false && eventSet_[PuCalo] == false && eventSet_[PF] == false && eventSet_[Calo] == false && eventSet_[T] == false)
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
      for(Int_t setIter = 0; setIter < 5; setIter++){
	if(eventSet_[setIter]) AlgTotTrk[setIter] = AlgTotTrk[setIter] + 1;
      }

      if(TMath::Abs(trkCollection.trkEta[trkEntry]) > 2.4){
	for(Int_t setIter = 0; setIter < 5; setIter++){
	  if(eventSet_[setIter]) AlgTrkEtaCut[setIter] = AlgTrkEtaCut[setIter] + 1;
	}
	continue;
      }

      if(trkCollection.trkPt[trkEntry] < 0.5){
	for(Int_t setIter = 0; setIter < 5; setIter++){
	  if(eventSet_[setIter]) AlgTrkPtCut[setIter] = AlgTrkPtCut[setIter] + 1;
	}
	continue;
      }

      if(!trkCollection.highPurity[trkEntry]){
	for(Int_t setIter = 0; setIter < 5; setIter++){
	  if(eventSet_[setIter]) AlgPurityCut[setIter] = AlgPurityCut[setIter] + 1;
	}
	continue;
      }

      if(TMath::Abs(trkCollection.trkDz1[trkEntry]/trkCollection.trkDzError1[trkEntry]) > 3){
	for(Int_t setIter = 0; setIter < 5; setIter++){
	  if(eventSet_[setIter]) AlgTrkDzCut[setIter] = AlgTrkDzCut[setIter] + 1;
	}
	continue;
      }

      if(TMath::Abs(trkCollection.trkDxy1[trkEntry]/trkCollection.trkDxyError1[trkEntry]) > 3){
	for(Int_t setIter = 0; setIter < 5; setIter++){
	  if(eventSet_[setIter]) AlgTrkDxyCut[setIter] = AlgTrkDxyCut[setIter] + 1;
	}
	continue;
      }

      if(trkCollection.trkPtError[trkEntry]/trkCollection.trkPt[trkEntry] > 0.1){
	for(Int_t setIter = 0; setIter < 5; setIter++){
	  if(eventSet_[setIter]) AlgTrkPtErrorCut[setIter] = AlgTrkPtErrorCut[setIter] + 1;
	}
	continue;
      }

      /*
      if(montecarlo){
	if(trkCollection.trkFake[trkEntry] == 1)
	  continue;

	if(trkCollection.trkStatus[trkEntry] < 0)
	  continue;
      }
      */

      trkPt_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPhi_[nTrk_] = trkCollection.trkPhi[trkEntry];
      trkEta_[nTrk_] = trkCollection.trkEta[trkEntry];

      trkPtCorrPuPF_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPtCorrPuCalo_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPtCorrT_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPtCorrPF_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPtCorrCalo_[nTrk_] = trkCollection.trkPt[trkEntry];
      
      //Grab proj. Pt Spectra For Tracks in each Event Subset

      if(eventSet_[PuPF]){
	trkPtPuPF_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[PuPF]));
	trkRLeadPuPF_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgLeadJtEta_[PuPF], AlgLeadJtPhi_[PuPF]);
	trkRSubLeadPuPF_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgSubLeadJtEta_[PuPF], AlgSubLeadJtPhi_[PuPF]);
        trkPuPFLeadDelPhi_[nTrk_] = getAbsDphi(AlgLeadJtPhi_[PuPF], trkCollection.trkPhi[trkEntry]);
     }
      else{
	trkPtPuPF_[nTrk_] = 0;
	trkRLeadPuPF_[nTrk_] = -1;
	trkRSubLeadPuPF_[nTrk_] = -1;
        trkPuPFLeadDelPhi_[nTrk_] = -10;
      }

      if(eventSet_[PuCalo]){
	trkPtPuCalo_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[PuCalo]));
	trkRLeadPuCalo_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgLeadJtEta_[PuCalo], AlgLeadJtPhi_[PuCalo]);
      	trkRSubLeadPuCalo_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgSubLeadJtEta_[PuCalo], AlgSubLeadJtPhi_[PuCalo]);
        trkPuCaloLeadDelPhi_[nTrk_] = getAbsDphi(AlgLeadJtPhi_[PuCalo], trkCollection.trkPhi[trkEntry]);
      }
      else{
	trkPtPuCalo_[nTrk_] = 0;
	trkRLeadPuCalo_[nTrk_] = -1;
	trkRSubLeadPuCalo_[nTrk_] = -1;
        trkPuCaloLeadDelPhi_[nTrk_] = -10;
      }

      if(eventSet_[PF]){
	trkPtPF_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[PF]));
	trkRLeadPF_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgLeadJtEta_[PF], AlgLeadJtPhi_[PF]);
	trkRSubLeadPF_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgSubLeadJtEta_[PF], AlgSubLeadJtPhi_[PF]);
        trkPFLeadDelPhi_[nTrk_] = getAbsDphi(AlgLeadJtPhi_[PF], trkCollection.trkPhi[trkEntry]);
      }
      else{
	trkPtPF_[nTrk_] = 0;
	trkRLeadPF_[nTrk_] = -1;
	trkRSubLeadPF_[nTrk_] = -1;
        trkPFLeadDelPhi_[nTrk_] = -10;
      }

      if(eventSet_[Calo]){
	trkPtCalo_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[Calo]));
	trkRLeadCalo_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgLeadJtEta_[Calo], AlgLeadJtPhi_[Calo]);
	trkRSubLeadCalo_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgSubLeadJtEta_[Calo], AlgSubLeadJtPhi_[Calo]);
        trkCaloLeadDelPhi_[nTrk_] = getAbsDphi(AlgLeadJtPhi_[Calo], trkCollection.trkPhi[trkEntry]);
      }
      else{
	trkPtCalo_[nTrk_] = 0;
	trkRLeadCalo_[nTrk_] = -1;
	trkRSubLeadCalo_[nTrk_] = -1;
        trkCaloLeadDelPhi_[nTrk_] = -10;
      }

      if(montecarlo && eventSet_[T]){
	trkPtT_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[T]));
	trkRLeadT_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgLeadJtEta_[T], AlgLeadJtPhi_[T]);
	trkRSubLeadT_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgSubLeadJtEta_[T], AlgSubLeadJtPhi_[T]);
        trkTLeadDelPhi_[nTrk_] = getAbsDphi(AlgLeadJtPhi_[T], trkCollection.trkPhi[trkEntry]);
      }
      else{
	trkPtT_[nTrk_] = 0;
	trkRLeadT_[nTrk_] = -1;
	trkRSubLeadT_[nTrk_] = -1;
        trkTLeadDelPhi_[nTrk_] = -10;
      }


      for(Int_t jtIter = 0; jtIter < 5; jtIter++){
	if(eventSet_[jtIter]){
	  getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProjF_[jtIter], rAlgImbPerpF_[jtIter], rAlgImbProj0_1_[jtIter], rAlgImbProj1_2_[jtIter], rAlgImbProj2_4_[jtIter], rAlgImbProj4_8_[jtIter], rAlgImbProj8_100_[jtIter]);

	  getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHemF_[jtIter], rAlgImbHem0_1_[jtIter], rAlgImbHem1_2_[jtIter], rAlgImbHem2_4_[jtIter], rAlgImbHem4_8_[jtIter], rAlgImbHem8_100_[jtIter]);


	  Float_t tempLeadDelR = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgLeadJtEta_[jtIter], AlgLeadJtPhi_[jtIter]);
	  Float_t tempSubLeadDelR = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgSubLeadJtEta_[jtIter], AlgSubLeadJtPhi_[jtIter]);

	  if(tempLeadDelR > 0 && tempSubLeadDelR > 0){
	    if(tempLeadDelR < .8 || tempSubLeadDelR < .8)
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProjCF_[jtIter], rAlgImbPerpCF_[jtIter], rAlgImbProjC0_1_[jtIter], rAlgImbProjC1_2_[jtIter], rAlgImbProjC2_4_[jtIter], rAlgImbProjC4_8_[jtIter], rAlgImbProjC8_100_[jtIter]);
	    else
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProjNCF_[jtIter], rAlgImbPerpNCF_[jtIter], rAlgImbProjNC0_1_[jtIter], rAlgImbProjNC1_2_[jtIter], rAlgImbProjNC2_4_[jtIter], rAlgImbProjNC4_8_[jtIter], rAlgImbProjNC8_100_[jtIter]);


	    if(tempLeadDelR < .5 || tempSubLeadDelR < .5)
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj05CF_[jtIter], rAlgImbPerp05CF_[jtIter], rAlgImbProj05C0_1_[jtIter], rAlgImbProj05C1_2_[jtIter],rAlgImbProj05C2_4_[jtIter], rAlgImbProj05C4_8_[jtIter], rAlgImbProj05C8_100_[jtIter]);
	    else if(tempLeadDelR < 1.0 || tempSubLeadDelR < 1.0)
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj510CF_[jtIter], rAlgImbPerp510CF_[jtIter], rAlgImbProj510C0_1_[jtIter], rAlgImbProj510C1_2_[jtIter],rAlgImbProj510C2_4_[jtIter], rAlgImbProj510C4_8_[jtIter], rAlgImbProj510C8_100_[jtIter]);
	    else if(tempLeadDelR < 1.5 || tempSubLeadDelR < 1.5)
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj1015CF_[jtIter], rAlgImbPerp1015CF_[jtIter], rAlgImbProj1015C0_1_[jtIter], rAlgImbProj1015C1_2_[jtIter],rAlgImbProj1015C2_4_[jtIter], rAlgImbProj1015C4_8_[jtIter], rAlgImbProj1015C8_100_[jtIter]);
	    else if(tempLeadDelR < 2.0 || tempSubLeadDelR < 2.0)
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj1520CF_[jtIter], rAlgImbPerp1520CF_[jtIter], rAlgImbProj1520C0_1_[jtIter], rAlgImbProj1520C1_2_[jtIter],rAlgImbProj1520C2_4_[jtIter], rAlgImbProj1520C4_8_[jtIter], rAlgImbProj1520C8_100_[jtIter]);
	    else if(tempLeadDelR < 2.5 || tempSubLeadDelR < 2.5)
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj2025CF_[jtIter], rAlgImbPerp2025CF_[jtIter], rAlgImbProj2025C0_1_[jtIter], rAlgImbProj2025C1_2_[jtIter],rAlgImbProj2025C2_4_[jtIter], rAlgImbProj2025C4_8_[jtIter], rAlgImbProj2025C8_100_[jtIter]);
	    else
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj2530CF_[jtIter], rAlgImbPerp2530CF_[jtIter], rAlgImbProj2530C0_1_[jtIter], rAlgImbProj2530C1_2_[jtIter],rAlgImbProj2530C2_4_[jtIter], rAlgImbProj2530C4_8_[jtIter], rAlgImbProj2530C8_100_[jtIter]);
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
      trkRMinPuPF_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[0]);
      trkRMinPuCalo_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[1]);
      trkRMinPF_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[2]);
      trkRMinCalo_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[3]);
    }

    if(montecarlo){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	trkRMinT_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[0], true);
      }
    }
    
    for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
      Int_t ptPos = getPtBin(trkPt_[trkEntry]);
	
      trkPtFactCalo_[trkEntry] = factorizedPtCorr_pp(trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinCalo_[trkEntry], CalophiEta_p[ptPos], Calopt_p[ptPos], CalodelR_p[ptPos], 5);
      trkPtCorrCalo_[trkEntry] = trkPt_[trkEntry]*(1 - factorizedPtCorr_pp(trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinCalo_[trkEntry], FakeCalophiEta_p[ptPos], FakeCalopt_p[ptPos], FakeCalodelR_p[ptPos], 5, false))/trkPtFactCalo_[trkEntry];
      

      Float_t tempCorr[5] = {trkPtCorrPuPF_[trkEntry], trkPtCorrPuCalo_[trkEntry], trkPtCorrPF_[trkEntry], trkPtCorrCalo_[trkEntry], trkPtCorrT_[trkEntry]};
      Float_t tempLeadR[5] = {trkRLeadPuPF_[trkEntry], trkRLeadPuCalo_[trkEntry], trkRLeadPF_[trkEntry], trkRLeadCalo_[trkEntry], trkRLeadT_[trkEntry]};
      Float_t tempSubLeadR[5] = {trkRSubLeadPuPF_[trkEntry], trkRSubLeadPuCalo_[trkEntry], trkRSubLeadPF_[trkEntry], trkRSubLeadCalo_[trkEntry], trkRLeadT_[trkEntry]};
      
      for(Int_t setIter = 0; setIter < 5; setIter++){
	if(eventSet_[setIter]){
	  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProjF_[setIter + 5], rAlgImbPerpF_[setIter + 5], rAlgImbProj0_1_[setIter + 5], rAlgImbProj1_2_[setIter + 5], rAlgImbProj2_4_[setIter + 5], rAlgImbProj4_8_[setIter + 5], rAlgImbProj8_100_[setIter + 5]);
	  
	  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHemF_[setIter + 5], rAlgImbHem0_1_[setIter + 5], rAlgImbHem1_2_[setIter + 5], rAlgImbHem2_4_[setIter + 5], rAlgImbHem4_8_[setIter + 5], rAlgImbHem8_100_[setIter + 5]);
	  
	  
	  if(tempLeadR[setIter] > 0 && tempSubLeadR[setIter] > 0){
	    if(tempLeadR[setIter] < .8 || tempSubLeadR[setIter] < .8)
	      getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProjCF_[setIter + 5], rAlgImbPerpCF_[setIter + 5], rAlgImbProjC0_1_[setIter + 5], rAlgImbProjC1_2_[setIter + 5], rAlgImbProjC2_4_[setIter + 5], rAlgImbProjC4_8_[setIter + 5], rAlgImbProjC8_100_[setIter + 5]);
	    else
	      getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProjNCF_[setIter + 5], rAlgImbPerpNCF_[setIter + 5], rAlgImbProjNC0_1_[setIter + 5], rAlgImbProjNC1_2_[setIter + 5], rAlgImbProjNC2_4_[setIter + 5], rAlgImbProjNC4_8_[setIter + 5], rAlgImbProjNC8_100_[setIter + 5]);		
	    
	    
	    if(tempLeadR[setIter] < .5 || tempSubLeadR[setIter] < .5)
	      getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj05CF_[setIter + 5], rAlgImbPerp05CF_[setIter + 5], rAlgImbProj05C0_1_[setIter + 5], rAlgImbProj05C1_2_[setIter + 5],rAlgImbProj05C2_4_[setIter + 5], rAlgImbProj05C4_8_[setIter + 5], rAlgImbProj05C8_100_[setIter + 5]);
	    else if(tempLeadR[setIter] < 1.0 || tempSubLeadR[setIter] < 1.0)
	      getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj510CF_[setIter + 5], rAlgImbPerp510CF_[setIter + 5], rAlgImbProj510C0_1_[setIter + 5], rAlgImbProj510C1_2_[setIter + 5],rAlgImbProj510C2_4_[setIter + 5], rAlgImbProj510C4_8_[setIter + 5], rAlgImbProj510C8_100_[setIter + 5]);
	    else if(tempLeadR[setIter] < 1.5 || tempSubLeadR[setIter] < 1.5)
	      getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj1015CF_[setIter + 5], rAlgImbPerp1015CF_[setIter + 5], rAlgImbProj1015C0_1_[setIter + 5], rAlgImbProj1015C1_2_[setIter + 5],rAlgImbProj1015C2_4_[setIter + 5], rAlgImbProj1015C4_8_[setIter + 5], rAlgImbProj1015C8_100_[setIter + 5]);
	    else if(tempLeadR[setIter] < 2.0 || tempSubLeadR[setIter] < 2.0)
	      getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj1520CF_[setIter + 5], rAlgImbPerp1520CF_[setIter + 5], rAlgImbProj1520C0_1_[setIter + 5], rAlgImbProj1520C1_2_[setIter + 5],rAlgImbProj1520C2_4_[setIter + 5], rAlgImbProj1520C4_8_[setIter + 5], rAlgImbProj1520C8_100_[setIter + 5]);
	    else if(tempLeadR[setIter] < 2.5 || tempSubLeadR[setIter] < 2.5)
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj2025CF_[setIter + 5], rAlgImbPerp2025CF_[setIter + 5], rAlgImbProj2025C0_1_[setIter + 5], rAlgImbProj2025C1_2_[setIter + 5],rAlgImbProj2025C2_4_[setIter + 5], rAlgImbProj2025C4_8_[setIter + 5], rAlgImbProj2025C8_100_[setIter + 5]);
	    else
	      getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj2530CF_[setIter + 5], rAlgImbPerp2530CF_[setIter + 5], rAlgImbProj2530C0_1_[setIter + 5], rAlgImbProj2530C1_2_[setIter + 5],rAlgImbProj2530C2_4_[setIter + 5], rAlgImbProj2530C4_8_[setIter + 5], rAlgImbProj2530C8_100_[setIter + 5]);
	    
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
	for(Int_t setIter = 0; setIter < 5; setIter++){
	  if(eventSet_[setIter]) AlgTotGen[setIter] = AlgTotGen[setIter] + 1;
	}

	if(genCollection.chg[genEntry] == 0){
	  for(Int_t setIter = 0; setIter < 5; setIter++){
	    if(eventSet_[setIter]) AlgGenChgCut[setIter] = AlgGenChgCut[setIter] + 1;
	  }
	  continue;
	}
  
	if(TMath::Abs(genCollection.eta[genEntry]) > 2.4){
	  for(Int_t setIter = 0; setIter < 5; setIter++){
	    if(eventSet_[setIter]) AlgGenEtaCut[setIter] = AlgGenEtaCut[setIter] + 1;
	  }
	  continue;
	}

	if(genCollection.pt[genEntry] < 0.5){
	  for(Int_t setIter = 0; setIter < 5; setIter++){
	    if(eventSet_[setIter]) AlgGenPtCut[setIter] = AlgGenPtCut[setIter] + 1;
	  }
	  continue;
	}

	genPt_[nGen_] = genCollection.pt[genEntry];
	genPhi_[nGen_] = genCollection.phi[genEntry];
	genEta_[nGen_] = genCollection.eta[genEntry];

	if(eventSet_[PuPF])
          genPtPuPF_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], AlgLeadJtPhi_[PuPF]));
	else
	  genPtPuPF_[nGen_] = 0;

	if(eventSet_[PuCalo])
          genPtPuCalo_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], AlgLeadJtPhi_[PuCalo]));
	else
	  genPtPuCalo_[nGen_] = 0;

	if(eventSet_[T])
          genPtT_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], AlgLeadJtPhi_[T]));
	else
	  genPtT_[nGen_] = 0;

	if(eventSet_[PF])
          genPtPF_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], AlgLeadJtPhi_[PF]));
	else
	  genPtPF_[nGen_] = 0;

	if(eventSet_[Calo])
          genPtCalo_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], AlgLeadJtPhi_[Calo]));
	else
	  genPtCalo_[nGen_] = 0;

	
	if(eventSet_[T])
	  genLeadDelPhi_[nGen_] = getAbsDphi(AlgLeadJtPhi_[T], genCollection.phi[genEntry]);
	else
	  genLeadDelPhi_[nGen_] = -10;

       
	for(Int_t setIter = 0; setIter < 5; setIter++){
	  if(eventSet_[setIter]){
	    getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProjF_[setIter], gAlgImbPerpF_[setIter], gAlgImbProj0_1_[setIter], gAlgImbProj1_2_[setIter], gAlgImbProj2_4_[setIter], gAlgImbProj4_8_[setIter], gAlgImbProj8_100_[setIter]);

	    getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHemF_[setIter], gAlgImbHem0_1_[setIter], gAlgImbHem1_2_[setIter], gAlgImbHem2_4_[setIter], gAlgImbHem4_8_[setIter], gAlgImbHem8_100_[setIter]);


	    Float_t tempLeadDelR = getDR(genEta_[nGen_], genPhi_[nGen_], AlgLeadJtEta_[setIter], AlgLeadJtPhi_[setIter]);
	    Float_t tempSubLeadDelR = getDR(genEta_[nGen_], genPhi_[nGen_], AlgSubLeadJtEta_[setIter], AlgSubLeadJtPhi_[setIter]);

	    if(tempLeadDelR > 0 && tempSubLeadDelR > 0){
	      if(tempLeadDelR < .8 || tempSubLeadDelR < .8)
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProjCF_[setIter], gAlgImbPerpCF_[setIter], gAlgImbProjC0_1_[setIter], gAlgImbProjC1_2_[setIter], gAlgImbProjC2_4_[setIter], gAlgImbProjC4_8_[setIter], gAlgImbProjC8_100_[setIter]);
	      else
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProjNCF_[setIter], gAlgImbPerpNCF_[setIter], gAlgImbProjNC0_1_[setIter], gAlgImbProjNC1_2_[setIter], gAlgImbProjNC2_4_[setIter], gAlgImbProjNC4_8_[setIter], gAlgImbProjNC8_100_[setIter]);
	      
	      
	      if(tempLeadDelR < .5 || tempSubLeadDelR < .5)
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj05CF_[setIter], gAlgImbPerp05CF_[setIter], gAlgImbProj05C0_1_[setIter], gAlgImbProj05C1_2_[setIter],gAlgImbProj05C2_4_[setIter], gAlgImbProj05C4_8_[setIter], gAlgImbProj05C8_100_[setIter]);
	      else if(tempLeadDelR < 1.0 || tempSubLeadDelR < 1.0)
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj510CF_[setIter], gAlgImbPerp510CF_[setIter], gAlgImbProj510C0_1_[setIter], gAlgImbProj510C1_2_[setIter],gAlgImbProj510C2_4_[setIter], gAlgImbProj510C4_8_[setIter], gAlgImbProj510C8_100_[setIter]);
	      else if(tempLeadDelR < 1.5 || tempSubLeadDelR < 1.5)
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj1015CF_[setIter], gAlgImbPerp1015CF_[setIter], gAlgImbProj1015C0_1_[setIter], gAlgImbProj1015C1_2_[setIter],gAlgImbProj1015C2_4_[setIter], gAlgImbProj1015C4_8_[setIter], gAlgImbProj1015C8_100_[setIter]);
	      else if(tempLeadDelR < 2.0 || tempSubLeadDelR < 2.0)
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj1520CF_[setIter], gAlgImbPerp1520CF_[setIter], gAlgImbProj1520C0_1_[setIter], gAlgImbProj1520C1_2_[setIter],gAlgImbProj1520C2_4_[setIter], gAlgImbProj1520C4_8_[setIter], gAlgImbProj1520C8_100_[setIter]);
	      else if(tempLeadDelR < 2.5 || tempSubLeadDelR < 2.5)
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj2025CF_[setIter], gAlgImbPerp2025CF_[setIter], gAlgImbProj2025C0_1_[setIter], gAlgImbProj2025C1_2_[setIter],gAlgImbProj2025C2_4_[setIter], gAlgImbProj2025C4_8_[setIter], gAlgImbProj2025C8_100_[setIter]);
	      else
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj2530CF_[setIter], gAlgImbPerp2530CF_[setIter], gAlgImbProj2530C0_1_[setIter], gAlgImbProj2530C1_2_[setIter],gAlgImbProj2530C2_4_[setIter], gAlgImbProj2530C4_8_[setIter], gAlgImbProj2530C8_100_[setIter]);

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

    jetTree_p->Fill();
    trackTree_p->Fill();

    if(montecarlo)
      genTree_p->Fill();
    }

  std::cout << "totEv: " << totEv << std::endl;
  Int_t tempTot = totEv - selectCut;
  std::cout << "selectCut: " << tempTot << std::endl;

  for(Int_t cutIter = 0; cutIter < 5; cutIter++){
    std::cout << std::endl;
    tempTot = totEv - selectCut - AlgLeadJtPtCut[cutIter];
    std::cout << "AlgLeadJtPtCut[" << cutIter << "]: " << tempTot << std::endl;
    tempTot = tempTot - AlgSubLeadJtPtCut[cutIter];
    std::cout << "AlgSubLeadJtPtCut[" << cutIter << "]: " << tempTot << std::endl;
    tempTot = tempTot - AlgDelPhiCut[cutIter];
    std::cout << "AlgDelPhiCut[" << cutIter << "]: " << tempTot << std::endl;
    tempTot = tempTot - AlgJtEtaCut[cutIter];
    std::cout << "AlgJtEtaCut[" << cutIter << "]: " << tempTot << std::endl;
  }

  for(Int_t cutIter = 0; cutIter < 5; cutIter++){
    std::cout << std::endl;
    std::cout << "AlgTotTrk[" << cutIter << "]: " << AlgTotTrk[cutIter] << std::endl;
    tempTot = AlgTotTrk[cutIter] - AlgTrkEtaCut[cutIter];
    std::cout << "AlgTrkEtaCut[" << cutIter << "]: " << tempTot << std::endl;
    tempTot = tempTot - AlgTrkPtCut[cutIter];
    std::cout << "AlgTrkPtCut[" << cutIter << "]: " << tempTot << std::endl;
    tempTot = tempTot - AlgPurityCut[cutIter];
    std::cout << "AlgPurityCut[" << cutIter << "]: " << tempTot << std::endl;
    tempTot = tempTot - AlgTrkDzCut[cutIter];
    std::cout << "AlgTrkDzCut[" << cutIter << "]: " << tempTot << std::endl;
    tempTot = tempTot - AlgTrkDxyCut[cutIter];
    std::cout << "AlgTrkDxyCut[" << cutIter << "]: " << tempTot << std::endl;
    tempTot = tempTot - AlgTrkPtErrorCut[cutIter];
    std::cout << "AlgTrkPtError[" << cutIter << "]: " << tempTot << std::endl;

    if(montecarlo){
      std::cout << std::endl;
      std::cout << "AlgTotGen[" << cutIter << "]: " << AlgTotGen[cutIter] << std::endl;
      tempTot = AlgTotGen[cutIter] - AlgGenChgCut[cutIter];
      std::cout << "AlgGenChgCut[" << cutIter << "]: " << tempTot << std::endl;
      tempTot = tempTot - AlgGenEtaCut[cutIter];
      std::cout << "AlgGenEtaCut[" << cutIter << "]: " << tempTot << std::endl;
      tempTot = tempTot - AlgGenPtCut[cutIter];
      std::cout << "AlgGenPtCut[" << cutIter << "]: " << tempTot << std::endl;
    }
  }

  outFile->cd();
  jetTree_p->Write();
  trackTree_p->Write();

  if(montecarlo)
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
  if(argc != 5)
    {
      std::cout << "Usage: sortForest <inputFile> <MCBool> <outputFile> <#>" << std::endl;
      return 1;
    }

  int rStatus = -1;

  rStatus = makeDiJetTTree_pp(argv[1], sampleType(atoi(argv[2])), argv[3], atoi(argv[4]));

  return rStatus;
}
