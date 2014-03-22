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
const Float_t jtDelPhiCut = 0;
const Float_t jtEtaCut = 1.6; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed

collisionType getCType(sampleType sType);

const char* algType[5] = {"PuPF", "PuCalo", "VsPF", "VsCalo", "T"};

void getLeadJt(Int_t& leadJtIndex, Int_t& subLeadJtIndex, Int_t& thirdJtIndex, Int_t& fourthJtIndex, Int_t& lPtCut, Int_t& sPtCut, Int_t& dPhiCut, Int_t& etaCut, Bool_t& setPass, Jets jtCollection)
{
  leadJtIndex = -1;
  subLeadJtIndex = -1;
  thirdJtIndex = -1;
  fourthJtIndex = -1;
  Float_t leadJtPt = 30;
  Float_t subLeadJtPt = 30;
  Float_t thirdJtPt = 30;
  Float_t fourthJtPt = 30;

  for(Int_t jtEntry = 0; jtEntry < jtCollection.nref; jtEntry++){
    if(jtCollection.jtpt[jtEntry] > leadJtPtCut && jtCollection.jtpt[jtEntry] > leadJtPt && TMath::Abs(jtCollection.jteta[jtEntry]) < jtEtaCut){

      fourthJtIndex = thirdJtIndex;
      fourthJtPt = thirdJtPt;
      thirdJtIndex = subLeadJtIndex;
      thirdJtPt = subLeadJtPt;
      
      subLeadJtIndex = leadJtIndex;
      subLeadJtPt = leadJtPt;
      leadJtIndex = jtEntry;
      leadJtPt = jtCollection.jtpt[jtEntry];
    }
    else if(jtCollection.jtpt[jtEntry] > subLeadJtPtCut && jtCollection.jtpt[jtEntry] > subLeadJtPt && TMath::Abs(jtCollection.jteta[jtEntry]) < jtEtaCut){

      fourthJtIndex = thirdJtIndex;
      fourthJtPt = thirdJtPt;
      thirdJtIndex = subLeadJtIndex;
      thirdJtPt = subLeadJtPt;

      subLeadJtIndex = jtEntry;
      subLeadJtPt = jtCollection.jtpt[jtEntry];
    }
    else if(jtCollection.jtpt[jtEntry] > thirdJtPt && TMath::Abs(jtCollection.jteta[jtEntry]) < jtEtaCut){
      fourthJtIndex = thirdJtIndex;
      fourthJtPt = thirdJtPt;

      thirdJtIndex = jtEntry;
      thirdJtPt = jtCollection.jtpt[jtEntry];
    }
    else if(jtCollection.jtpt[jtEntry] > fourthJtPt && TMath::Abs(jtCollection.jteta[jtEntry]) < jtEtaCut){
      fourthJtIndex = jtEntry;
      fourthJtPt = jtCollection.jtpt[jtEntry];
    }
  }

  if(leadJtIndex < 0)
    lPtCut++; 
  else if(subLeadJtIndex < 0)
    sPtCut++;
  else if(getAbsDphi(jtCollection.jtphi[leadJtIndex], jtCollection.jtphi[subLeadJtIndex]) < jtDelPhiCut)
    dPhiCut++;
  else if(TMath::Abs(jtCollection.jteta[leadJtIndex]) > jtEtaCut || TMath::Abs(jtCollection.jteta[subLeadJtIndex]) > jtEtaCut)
    etaCut++;
  else
    setPass = true;

  return;
}


Float_t getFlippedPhi(Float_t inPhi)
{
  Float_t outPhi;

  if(TMath::Abs(inPhi) > TMath::Pi()){
    std::cout << "getFlippedPhi: inPhi is outside accepted range, return -10" << std::endl;
    return -10;
  }
  else if(inPhi > 0)
    outPhi = inPhi - TMath::Pi();
  else
    outPhi = inPhi + TMath::Pi();

  return outPhi;
}


Bool_t sameSign(Float_t num1, Float_t num2){
  if((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0)) return true;

  return false;
}


Float_t getAvePhi(Float_t inLeadPhi, Float_t inSubLeadPhi)
{
  Float_t flipPhi = getFlippedPhi(inSubLeadPhi);
  Float_t avePhi;

  if(sameSign(inLeadPhi, flipPhi) || (TMath::Abs(inLeadPhi) < TMath::Pi()/2 && TMath::Abs(flipPhi) < TMath::Pi()/2))
    avePhi = (flipPhi + inLeadPhi)/2;
  else if(TMath::Abs(inLeadPhi) > TMath::Pi()/2 && TMath::Abs(flipPhi) > TMath::Pi()/2){
    avePhi = (flipPhi + inLeadPhi)/2;
    if(avePhi > 0)
      avePhi = TMath::Pi() - avePhi;
    else
      avePhi = -TMath::Pi() - avePhi;
  }
  else{
    avePhi = 0.;
  }

  return avePhi;
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


int makeDiJetTTree(string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_CFMSKIM.root", Int_t num = 0)
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

  TFile *centHistFile_p = new TFile("centHist_eventSelect.root", "READ");
  TH1F *hist_DataOverMC_p[4];

  TFile *centHistFile_2pi3_p = new TFile("centHist_eventSelect_2pi3.root", "READ");
  TH1F *hist_DataOverMC_2pi3_p[4];

  if(montecarlo){
    for(Int_t algIter = 0; algIter < 4; algIter++){
      hist_DataOverMC_p[algIter] = (TH1F*)centHistFile_p->Get(Form("hiBin_%s_DataOverMC_h", algType[algIter]));
      hist_DataOverMC_2pi3_p[algIter] = (TH1F*)centHistFile_2pi3_p->Get(Form("hiBin_%s_DataOverMC_h", algType[algIter]));
    }
  }

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
  c->hasAkVs3PFJetTree = true;
  c->hasAkVs3CaloJetTree = true;

  if(montecarlo)
    c->hasGenParticleTree = true;

  Long64_t nentries = c->GetEntries();

  std::cout << nentries << std::endl;

  Int_t totEv = 0;
  Int_t selectCut = 0;
  Int_t vzCut = 0;

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

  std::cout << "Cuts, Lead/Sublead Pt, delphi, eta: " << leadJtPtCut << ", " << subLeadJtPtCut << ", " << jtDelPhiCut << ", " << jtEtaCut << std::endl; 

  for(Long64_t jentry = 0; jentry < 1500/*nentries*/; jentry++){
    c->GetEntry(jentry);

    totEv++;

    if(jentry%1000 == 0)
      std::cout << jentry << std::endl;

    if(!c->selectEvent()){
      selectCut++;
      continue;
    }

    if(TMath::Abs(c->evt.vz) > 15){
      vzCut++;
      continue;
    }

    InitJetVar(montecarlo);

    //particle flow

    Jets AlgJtCollection[4] = {c->akPu3PF, c->akPu3Calo, c->akVs3PF, c->akVs3Calo};

    Int_t tempJtIndex[4] = {-1, -1, -1, -1};

    for(Int_t collIter = 0; collIter < 4; collIter++){
      getLeadJt(tempJtIndex[0], tempJtIndex[1], tempJtIndex[2], tempJtIndex[3], AlgLeadJtPtCut[collIter], AlgSubLeadJtPtCut[collIter], AlgDelPhiCut[collIter], AlgJtEtaCut[collIter], eventSet_[collIter], AlgJtCollection[collIter]);
      
      if(eventSet_[collIter]){
	AlgLeadJtPt_[collIter] = AlgJtCollection[collIter].jtpt[tempJtIndex[0]];
	AlgLeadJtPhi_[collIter] = AlgJtCollection[collIter].jtphi[tempJtIndex[0]];
	AlgLeadJtEta_[collIter] = AlgJtCollection[collIter].jteta[tempJtIndex[0]];

	AlgSubLeadJtPt_[collIter] = AlgJtCollection[collIter].jtpt[tempJtIndex[1]];
	AlgSubLeadJtPhi_[collIter] = AlgJtCollection[collIter].jtphi[tempJtIndex[1]];
	AlgSubLeadJtEta_[collIter] = AlgJtCollection[collIter].jteta[tempJtIndex[1]];

	AlgJtAvePhi_[collIter] = getAvePhi(AlgLeadJtPhi_[collIter], AlgSubLeadJtPhi_[collIter]);
	AlgJtDelPhi_[collIter] = getAbsDphi(AlgLeadJtPhi_[collIter], AlgSubLeadJtPhi_[collIter]);
	AlgJtAsymm_[collIter] = (AlgLeadJtPt_[collIter] - AlgSubLeadJtPt_[collIter])/(AlgLeadJtPt_[collIter] + AlgSubLeadJtPt_[collIter]);

	if(montecarlo){
	  AlgLeadRefPt_[collIter] = AlgJtCollection[collIter].refpt[tempJtIndex[0]];
	  AlgLeadRefPhi_[collIter] = AlgJtCollection[collIter].refphi[tempJtIndex[0]];
	  AlgLeadRefEta_[collIter] = AlgJtCollection[collIter].refeta[tempJtIndex[0]];

	  AlgSubLeadRefPt_[collIter] = AlgJtCollection[collIter].refpt[tempJtIndex[1]];
	  AlgSubLeadRefPhi_[collIter] = AlgJtCollection[collIter].refphi[tempJtIndex[1]];
	  AlgSubLeadRefEta_[collIter] = AlgJtCollection[collIter].refeta[tempJtIndex[1]];


	  AlgRefDelPhi_[collIter] = getAbsDphi(AlgLeadRefPhi_[collIter], AlgSubLeadRefPhi_[collIter]);
	  AlgRefAsymm_[collIter] = (AlgLeadRefPt_[collIter] - AlgSubLeadRefPt_[collIter])/(AlgLeadRefPt_[collIter] + AlgSubLeadRefPt_[collIter]);
	}
      }

      if(tempJtIndex[2] >= 0){
	AlgThirdJtPt_[collIter] = AlgJtCollection[collIter].jtpt[tempJtIndex[2]];
        AlgThirdJtPhi_[collIter] = AlgJtCollection[collIter].jtphi[tempJtIndex[2]];
        AlgThirdJtEta_[collIter] = AlgJtCollection[collIter].jteta[tempJtIndex[2]];

	if(montecarlo){
	  AlgThirdRefPt_[collIter] = AlgJtCollection[collIter].refpt[tempJtIndex[2]];
          AlgThirdRefPhi_[collIter] = AlgJtCollection[collIter].refphi[tempJtIndex[2]];
          AlgThirdRefEta_[collIter] = AlgJtCollection[collIter].refeta[tempJtIndex[2]];
	}
      }

      if(tempJtIndex[3] >= 0){
	AlgFourthJtPt_[collIter] = AlgJtCollection[collIter].jtpt[tempJtIndex[3]];
        AlgFourthJtPhi_[collIter] = AlgJtCollection[collIter].jtphi[tempJtIndex[3]];
        AlgFourthJtEta_[collIter] = AlgJtCollection[collIter].jteta[tempJtIndex[3]];

	if(montecarlo){
	  AlgFourthRefPt_[collIter] = AlgJtCollection[collIter].refpt[tempJtIndex[3]];
          AlgFourthRefPhi_[collIter] = AlgJtCollection[collIter].refphi[tempJtIndex[3]];
          AlgFourthRefEta_[collIter] = AlgJtCollection[collIter].refeta[tempJtIndex[3]];
	}
      }

    }

    //truth, doesn't work w/ getLeadJt because truth doesnt get its own tree

    if(montecarlo){
      Int_t leadJtIndex = -1;
      Int_t subLeadJtIndex = -1;
      Int_t thirdJtIndex = -1;
      Int_t fourthJtIndex = -1;
      AlgLeadJtPt_[T] = 30;
      AlgSubLeadJtPt_[T] = 30;
      AlgThirdJtPt_[T] = 30;
      AlgFourthJtPt_[T] = 30;
      for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.ngen; jtEntry++){
	if(c->akPu3PF.genpt[jtEntry] > leadJtPtCut && c->akPu3PF.genpt[jtEntry] > AlgLeadJtPt_[T] && TMath::Abs(c->akPu3PF.geneta[jtEntry]) < 1.6){

	  fourthJtIndex = thirdJtIndex;
	  AlgFourthJtPt_[T] = AlgThirdJtPt_[T];
	  thirdJtIndex = subLeadJtIndex;
	  AlgThirdJtPt_[T] = AlgSubLeadJtPt_[T];

	  subLeadJtIndex = leadJtIndex;
	  AlgSubLeadJtPt_[T] = AlgLeadJtPt_[T];
	  leadJtIndex = jtEntry;
	  AlgLeadJtPt_[T] = c->akPu3PF.genpt[jtEntry];
	}
	else if(c->akPu3PF.genpt[jtEntry] > subLeadJtPtCut && c->akPu3PF.genpt[jtEntry] > AlgSubLeadJtPt_[T] && TMath::Abs(c->akPu3PF.geneta[jtEntry]) < 1.6){
	  fourthJtIndex = thirdJtIndex;
          AlgFourthJtPt_[T] = AlgThirdJtPt_[T];
          thirdJtIndex = subLeadJtIndex;
          AlgThirdJtPt_[T] = AlgSubLeadJtPt_[T];

	  subLeadJtIndex = jtEntry;
	  AlgSubLeadJtPt_[T] = c->akPu3PF.genpt[jtEntry];
	}
	else if(c->akPu3PF.genpt[jtEntry] > AlgThirdJtPt_[T] && TMath::Abs(c->akPu3PF.geneta[jtEntry]) < 1.6){
	  fourthJtIndex = thirdJtIndex;
          AlgFourthJtPt_[T] = AlgThirdJtPt_[T];

	  thirdJtIndex = jtEntry;
	  AlgThirdJtPt_[T] = c->akPu3PF.genpt[jtEntry];
	} 
	else if(c->akPu3PF.genpt[jtEntry] > AlgFourthJtPt_[T] && TMath::Abs(c->akPu3PF.geneta[jtEntry]) < 1.6){
	  fourthJtIndex = jtEntry;
	  AlgFourthJtPt_[T] = c->akPu3PF.genpt[jtEntry];
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

	AlgJtAvePhi_[T] = getAvePhi(AlgLeadJtPhi_[T], AlgSubLeadJtPhi_[T]);
	AlgJtDelPhi_[T] = getAbsDphi(AlgLeadJtPhi_[T], AlgSubLeadJtPhi_[T]);
	AlgJtAsymm_[T] = (AlgLeadJtPt_[T] - AlgSubLeadJtPt_[T])/(AlgLeadJtPt_[T] + AlgSubLeadJtPt_[T]);
	
	eventSet_[T] = true;

	if(thirdJtIndex >= 0){
	  AlgThirdJtPhi_[T] = c->akPu3PF.genphi[thirdJtIndex];
	  AlgThirdJtEta_[T] = c->akPu3PF.geneta[thirdJtIndex];
	}
	else
	  AlgThirdJtPt_[T] = -10;

	if(fourthJtIndex >= 0){
	  AlgFourthJtPhi_[T] = c->akPu3PF.genphi[fourthJtIndex];
	  AlgFourthJtEta_[T] = c->akPu3PF.geneta[fourthJtIndex];
	}
	else
	  AlgFourthJtPt_[T] = -10;

      }
    }
    
    if(eventSet_[PuPF] == false && eventSet_[PuCalo] == false && eventSet_[VsPF] == false && eventSet_[VsCalo] == false && eventSet_[T] == false)
      continue;

    run_ = c->evt.run;
    evt_ = c->akPu3PF.evt;
    lumi_ = c->evt.lumi;
    hiBin_ = c->evt.hiBin;

    if(montecarlo){
      for(Int_t algIter = 0; algIter < 4; algIter++){
	setWeight_[algIter] = hist_DataOverMC_p[algIter]->GetBinContent(hist_DataOverMC_p[algIter]->FindBin(hiBin_));
	setWeight_2pi3_[algIter] = hist_DataOverMC_2pi3_p[algIter]->GetBinContent(hist_DataOverMC_2pi3_p[algIter]->FindBin(hiBin_));
      }
    }

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
      trkPtCorrVsPF_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPtCorrVsCalo_[nTrk_] = trkCollection.trkPt[trkEntry];
      
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

      if(eventSet_[VsPF]){
	trkPtVsPF_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[VsPF]));
	trkRLeadVsPF_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgLeadJtEta_[VsPF], AlgLeadJtPhi_[VsPF]);
	trkRSubLeadVsPF_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgSubLeadJtEta_[VsPF], AlgSubLeadJtPhi_[VsPF]);
        trkVsPFLeadDelPhi_[nTrk_] = getAbsDphi(AlgLeadJtPhi_[VsPF], trkCollection.trkPhi[trkEntry]);
      }
      else{
	trkPtVsPF_[nTrk_] = 0;
	trkRLeadVsPF_[nTrk_] = -1;
	trkRSubLeadVsPF_[nTrk_] = -1;
        trkVsPFLeadDelPhi_[nTrk_] = -10;
      }

      if(eventSet_[VsCalo]){
	trkPtVsCalo_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[VsCalo]));
	trkRLeadVsCalo_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgLeadJtEta_[VsCalo], AlgLeadJtPhi_[VsCalo]);
	trkRSubLeadVsCalo_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgSubLeadJtEta_[VsCalo], AlgSubLeadJtPhi_[VsCalo]);
        trkVsCaloLeadDelPhi_[nTrk_] = getAbsDphi(AlgLeadJtPhi_[VsCalo], trkCollection.trkPhi[trkEntry]);
      }
      else{
	trkPtVsCalo_[nTrk_] = 0;
	trkRLeadVsCalo_[nTrk_] = -1;
	trkRSubLeadVsCalo_[nTrk_] = -1;
        trkVsCaloLeadDelPhi_[nTrk_] = -10;
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

	  Float_t tempHold = 0;

	  getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjAF_[jtIter], tempHold, rAlgImbProjA0_1_[jtIter], rAlgImbProjA1_2_[jtIter], rAlgImbProjA2_4_[jtIter], rAlgImbProjA4_8_[jtIter], rAlgImbProjA8_100_[jtIter]);

	  Float_t tempLeadDelR = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgLeadJtEta_[jtIter], AlgLeadJtPhi_[jtIter]);
	  Float_t tempSubLeadDelR = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], AlgSubLeadJtEta_[jtIter], AlgSubLeadJtPhi_[jtIter]);

	  if(tempLeadDelR > 0 && tempSubLeadDelR > 0){
	    if(tempLeadDelR < .8 || tempSubLeadDelR < .8){
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProjCF_[jtIter], rAlgImbPerpCF_[jtIter], rAlgImbProjC0_1_[jtIter], rAlgImbProjC1_2_[jtIter], rAlgImbProjC2_4_[jtIter], rAlgImbProjC4_8_[jtIter], rAlgImbProjC8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHemCF_[jtIter], rAlgImbHemC0_1_[jtIter], rAlgImbHemC1_2_[jtIter], rAlgImbHemC2_4_[jtIter], rAlgImbHemC4_8_[jtIter], rAlgImbHemC8_100_[jtIter]);

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjACF_[jtIter], tempHold, rAlgImbProjAC0_1_[jtIter], rAlgImbProjAC1_2_[jtIter], rAlgImbProjAC2_4_[jtIter], rAlgImbProjAC4_8_[jtIter], rAlgImbProjAC8_100_[jtIter]);
	    }
	    else{
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProjNCF_[jtIter], rAlgImbPerpNCF_[jtIter], rAlgImbProjNC0_1_[jtIter], rAlgImbProjNC1_2_[jtIter], rAlgImbProjNC2_4_[jtIter], rAlgImbProjNC4_8_[jtIter], rAlgImbProjNC8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHemNCF_[jtIter], rAlgImbHemNC0_1_[jtIter], rAlgImbHemNC1_2_[jtIter], rAlgImbHemNC2_4_[jtIter], rAlgImbHemNC4_8_[jtIter], rAlgImbHemNC8_100_[jtIter]);

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjANCF_[jtIter], tempHold, rAlgImbProjANC0_1_[jtIter], rAlgImbProjANC1_2_[jtIter], rAlgImbProjANC2_4_[jtIter], rAlgImbProjANC4_8_[jtIter], rAlgImbProjANC8_100_[jtIter]);
	    }

	  
	    if(tempLeadDelR < .20 || tempSubLeadDelR < .20){

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj1CF_[jtIter], tempHold, rAlgImbProj1C0_1_[jtIter], rAlgImbProj1C1_2_[jtIter], rAlgImbProj1C2_4_[jtIter], rAlgImbProj1C4_8_[jtIter], rAlgImbProj1C8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHem1CF_[jtIter], rAlgImbHem1C0_1_[jtIter], rAlgImbHem1C1_2_[jtIter], rAlgImbHem1C2_4_[jtIter], rAlgImbHem1C4_8_[jtIter], rAlgImbHem1C8_100_[jtIter]);
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA1CF_[jtIter], tempHold, rAlgImbProjA1C0_1_[jtIter], rAlgImbProjA1C1_2_[jtIter], rAlgImbProjA1C2_4_[jtIter], rAlgImbProjA1C4_8_[jtIter], rAlgImbProjA1C8_100_[jtIter]);

	    }
	    else if(tempLeadDelR < .40 || tempSubLeadDelR < .40){

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj2CF_[jtIter], tempHold, rAlgImbProj2C0_1_[jtIter], rAlgImbProj2C1_2_[jtIter], rAlgImbProj2C2_4_[jtIter], rAlgImbProj2C4_8_[jtIter], rAlgImbProj2C8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHem2CF_[jtIter], rAlgImbHem2C0_1_[jtIter], rAlgImbHem2C1_2_[jtIter], rAlgImbHem2C2_4_[jtIter], rAlgImbHem2C4_8_[jtIter], rAlgImbHem2C8_100_[jtIter]);
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA2CF_[jtIter], tempHold, rAlgImbProjA2C0_1_[jtIter], rAlgImbProjA2C1_2_[jtIter], rAlgImbProjA2C2_4_[jtIter], rAlgImbProjA2C4_8_[jtIter], rAlgImbProjA2C8_100_[jtIter]);

	    }
	    else if(tempLeadDelR < .60 || tempSubLeadDelR < .60){

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj3CF_[jtIter], tempHold, rAlgImbProj3C0_1_[jtIter], rAlgImbProj3C1_2_[jtIter], rAlgImbProj3C2_4_[jtIter], rAlgImbProj3C4_8_[jtIter], rAlgImbProj3C8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHem3CF_[jtIter], rAlgImbHem3C0_1_[jtIter], rAlgImbHem3C1_2_[jtIter], rAlgImbHem3C2_4_[jtIter], rAlgImbHem3C4_8_[jtIter], rAlgImbHem3C8_100_[jtIter]);
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA3CF_[jtIter], tempHold, rAlgImbProjA3C0_1_[jtIter], rAlgImbProjA3C1_2_[jtIter], rAlgImbProjA3C2_4_[jtIter], rAlgImbProjA3C4_8_[jtIter], rAlgImbProjA3C8_100_[jtIter]);

	    }
	    else if(tempLeadDelR < .80 || tempSubLeadDelR < .80){

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj4CF_[jtIter], tempHold, rAlgImbProj4C0_1_[jtIter], rAlgImbProj4C1_2_[jtIter], rAlgImbProj4C2_4_[jtIter], rAlgImbProj4C4_8_[jtIter], rAlgImbProj4C8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHem4CF_[jtIter], rAlgImbHem4C0_1_[jtIter], rAlgImbHem4C1_2_[jtIter], rAlgImbHem4C2_4_[jtIter], rAlgImbHem4C4_8_[jtIter], rAlgImbHem4C8_100_[jtIter]);
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA4CF_[jtIter], tempHold, rAlgImbProjA4C0_1_[jtIter], rAlgImbProjA4C1_2_[jtIter], rAlgImbProjA4C2_4_[jtIter], rAlgImbProjA4C4_8_[jtIter], rAlgImbProjA4C8_100_[jtIter]);

	    }
	    else if(tempLeadDelR < 1.0 || tempSubLeadDelR < 1.0){

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj5CF_[jtIter], tempHold, rAlgImbProj5C0_1_[jtIter], rAlgImbProj5C1_2_[jtIter], rAlgImbProj5C2_4_[jtIter], rAlgImbProj5C4_8_[jtIter], rAlgImbProj5C8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHem5CF_[jtIter], rAlgImbHem5C0_1_[jtIter], rAlgImbHem5C1_2_[jtIter], rAlgImbHem5C2_4_[jtIter], rAlgImbHem5C4_8_[jtIter], rAlgImbHem5C8_100_[jtIter]);
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA5CF_[jtIter], tempHold, rAlgImbProjA5C0_1_[jtIter], rAlgImbProjA5C1_2_[jtIter], rAlgImbProjA5C2_4_[jtIter], rAlgImbProjA5C4_8_[jtIter], rAlgImbProjA5C8_100_[jtIter]);

	    }
	    else if(tempLeadDelR < 1.2 || tempSubLeadDelR < 1.2){

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj6CF_[jtIter], tempHold, rAlgImbProj6C0_1_[jtIter], rAlgImbProj6C1_2_[jtIter], rAlgImbProj6C2_4_[jtIter], rAlgImbProj6C4_8_[jtIter], rAlgImbProj6C8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHem6CF_[jtIter], rAlgImbHem6C0_1_[jtIter], rAlgImbHem6C1_2_[jtIter], rAlgImbHem6C2_4_[jtIter], rAlgImbHem6C4_8_[jtIter], rAlgImbHem6C8_100_[jtIter]);
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA6CF_[jtIter], tempHold, rAlgImbProjA6C0_1_[jtIter], rAlgImbProjA6C1_2_[jtIter], rAlgImbProjA6C2_4_[jtIter], rAlgImbProjA6C4_8_[jtIter], rAlgImbProjA6C8_100_[jtIter]);

	    }
	    else if(tempLeadDelR < 1.4 || tempSubLeadDelR < 1.4){

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj7CF_[jtIter], tempHold, rAlgImbProj7C0_1_[jtIter], rAlgImbProj7C1_2_[jtIter], rAlgImbProj7C2_4_[jtIter], rAlgImbProj7C4_8_[jtIter], rAlgImbProj7C8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHem7CF_[jtIter], rAlgImbHem7C0_1_[jtIter], rAlgImbHem7C1_2_[jtIter], rAlgImbHem7C2_4_[jtIter], rAlgImbHem7C4_8_[jtIter], rAlgImbHem7C8_100_[jtIter]);
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA7CF_[jtIter], tempHold, rAlgImbProjA7C0_1_[jtIter], rAlgImbProjA7C1_2_[jtIter], rAlgImbProjA7C2_4_[jtIter], rAlgImbProjA7C4_8_[jtIter], rAlgImbProjA7C8_100_[jtIter]);

	    }
	    else if(tempLeadDelR < 1.6 || tempSubLeadDelR < 1.6){

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj8CF_[jtIter], tempHold, rAlgImbProj8C0_1_[jtIter], rAlgImbProj8C1_2_[jtIter], rAlgImbProj8C2_4_[jtIter], rAlgImbProj8C4_8_[jtIter], rAlgImbProj8C8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHem8CF_[jtIter], rAlgImbHem8C0_1_[jtIter], rAlgImbHem8C1_2_[jtIter], rAlgImbHem8C2_4_[jtIter], rAlgImbHem8C4_8_[jtIter], rAlgImbHem8C8_100_[jtIter]);
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA8CF_[jtIter], tempHold, rAlgImbProjA8C0_1_[jtIter], rAlgImbProjA8C1_2_[jtIter], rAlgImbProjA8C2_4_[jtIter], rAlgImbProjA8C4_8_[jtIter], rAlgImbProjA8C8_100_[jtIter]);

	    }
	    else if(tempLeadDelR < 1.8 || tempSubLeadDelR < 1.8){

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj9CF_[jtIter], tempHold, rAlgImbProj9C0_1_[jtIter], rAlgImbProj9C1_2_[jtIter], rAlgImbProj9C2_4_[jtIter], rAlgImbProj9C4_8_[jtIter], rAlgImbProj9C8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHem9CF_[jtIter], rAlgImbHem9C0_1_[jtIter], rAlgImbHem9C1_2_[jtIter], rAlgImbHem9C2_4_[jtIter], rAlgImbHem9C4_8_[jtIter], rAlgImbHem9C8_100_[jtIter]);
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA9CF_[jtIter], tempHold, rAlgImbProjA9C0_1_[jtIter], rAlgImbProjA9C1_2_[jtIter], rAlgImbProjA9C2_4_[jtIter], rAlgImbProjA9C4_8_[jtIter], rAlgImbProjA9C8_100_[jtIter]);

	    }
	    else{

	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbProj10CF_[jtIter], tempHold, rAlgImbProj10C0_1_[jtIter], rAlgImbProj10C1_2_[jtIter], rAlgImbProj10C2_4_[jtIter], rAlgImbProj10C4_8_[jtIter], rAlgImbProj10C8_100_[jtIter]);
	      getPtHem(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgLeadJtPhi_[jtIter], rAlgImbHem10CF_[jtIter], rAlgImbHem10C0_1_[jtIter], rAlgImbHem10C1_2_[jtIter], rAlgImbHem10C2_4_[jtIter], rAlgImbHem10C4_8_[jtIter], rAlgImbHem10C8_100_[jtIter]);
	      getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA10CF_[jtIter], tempHold, rAlgImbProjA10C0_1_[jtIter], rAlgImbProjA10C1_2_[jtIter], rAlgImbProjA10C2_4_[jtIter], rAlgImbProjA10C4_8_[jtIter], rAlgImbProjA10C8_100_[jtIter]);

	    }	  
	    
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
      trkRMinVsPF_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[2]);
      trkRMinVsCalo_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[3]);
    }
    
    if(montecarlo){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	trkRMinT_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[3], true);
      }
    }
    
    Int_t hiBinDiv[5] = {20, 40, 60, 100, 200};
    Int_t hiSetEff[15] = {0, 5, 10, 1, 6, 11, 2, 7, 12, 3, 8, 12, 4, 9, 12};
  
    for(Int_t hiBinIter = 0; hiBinIter < 5; hiBinIter++){
      if(hiBin_ < hiBinDiv[hiBinIter]){
	for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	  Int_t ptPos = getPtBin(trkPt_[trkEntry], hiSetEff[hiBinIter*3], hiSetEff[hiBinIter*3 + 1], hiSetEff[hiBinIter*3 + 2], 13);
	
	  trkPtFactPuPF_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPuPF_[trkEntry], PuPFcent_p[ptPos], PuPFphiEta_p[ptPos], PuPFpt_p[ptPos], PuPFdelR_p[ptPos], 3);
	  trkPtCorrPuPF_[trkEntry] = trkPt_[trkEntry]/trkPtFactPuPF_[trkEntry];

	  trkPtFactPuCalo_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPuCalo_[trkEntry], PuCalocent_p[ptPos], PuCalophiEta_p[ptPos], PuCalopt_p[ptPos], PuCalodelR_p[ptPos], 5);
	  trkPtFactPuCalo_[trkEntry] = (1 - factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPuCalo_[trkEntry], FakePuCalocent_p[ptPos], FakePuCalophiEta_p[ptPos], FakePuCalopt_p[ptPos], FakePuCalodelR_p[ptPos], 5, false))/trkPtFactPuCalo_[trkEntry];
	  trkPtCorrPuCalo_[trkEntry] = trkPt_[trkEntry]*trkPtFactPuCalo_[trkEntry];

	  trkPtFactVsCalo_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinVsCalo_[trkEntry], VsCalocent_p[ptPos], VsCalophiEta_p[ptPos], VsCalopt_p[ptPos], VsCalodelR_p[ptPos], 5);
	  trkPtFactVsCalo_[trkEntry] = (1 - factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinVsCalo_[trkEntry], FakeVsCalocent_p[ptPos], FakeVsCalophiEta_p[ptPos], FakeVsCalopt_p[ptPos], FakeVsCalodelR_p[ptPos], 5, false))/trkPtFactVsCalo_[trkEntry];
	  trkPtCorrVsCalo_[trkEntry] = trkPt_[trkEntry]*trkPtFactVsCalo_[trkEntry];

	  if(montecarlo){
	    trkPtFactT_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinT_[trkEntry], VsCalocent_p[ptPos], VsCalophiEta_p[ptPos], VsCalopt_p[ptPos], VsCalodelR_p[ptPos], 5);
	    trkPtCorrT_[trkEntry] = trkPt_[trkEntry]*(1 - factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinT_[trkEntry], FakeVsCalocent_p[ptPos], FakeVsCalophiEta_p[ptPos], FakeVsCalopt_p[ptPos], FakeVsCalodelR_p[ptPos], 5, false))/trkPtFactT_[trkEntry];
	  }


	  Float_t tempCorr[5] = {trkPtCorrPuPF_[trkEntry], trkPtCorrPuCalo_[trkEntry], trkPtCorrVsPF_[trkEntry], trkPtCorrVsCalo_[trkEntry], trkPtCorrT_[trkEntry]};
	  Float_t tempLeadR[5] = {trkRLeadPuPF_[trkEntry], trkRLeadPuCalo_[trkEntry], trkRLeadVsPF_[trkEntry], trkRLeadVsCalo_[trkEntry], trkRLeadT_[trkEntry]};
	  Float_t tempSubLeadR[5] = {trkRSubLeadPuPF_[trkEntry], trkRSubLeadPuCalo_[trkEntry], trkRSubLeadVsPF_[trkEntry], trkRSubLeadVsCalo_[trkEntry], trkRSubLeadT_[trkEntry]};


	  for(Int_t setIter = 0; setIter < 5; setIter++){
	    if(eventSet_[setIter]){
	      getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProjF_[setIter + 5], rAlgImbPerpF_[setIter + 5], rAlgImbProj0_1_[setIter + 5], rAlgImbProj1_2_[setIter + 5], rAlgImbProj2_4_[setIter + 5], rAlgImbProj4_8_[setIter + 5], rAlgImbProj8_100_[setIter + 5]);

	      getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHemF_[setIter + 5], rAlgImbHem0_1_[setIter + 5], rAlgImbHem1_2_[setIter + 5], rAlgImbHem2_4_[setIter + 5], rAlgImbHem4_8_[setIter + 5], rAlgImbHem8_100_[setIter + 5]);

	      Float_t tempHold = 0;

	      getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjAF_[setIter + 5], tempHold, rAlgImbProjA0_1_[setIter + 5], rAlgImbProjA1_2_[setIter + 5], rAlgImbProjA2_4_[setIter + 5], rAlgImbProjA4_8_[setIter + 5], rAlgImbProjA8_100_[setIter + 5]);

	      if(tempLeadR[setIter] > 0 && tempSubLeadR[setIter] > 0){
		if(tempLeadR[setIter] < .8 || tempSubLeadR[setIter] < .8){
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProjCF_[setIter + 5], rAlgImbPerpCF_[setIter + 5], rAlgImbProjC0_1_[setIter + 5], rAlgImbProjC1_2_[setIter + 5], rAlgImbProjC2_4_[setIter + 5], rAlgImbProjC4_8_[setIter + 5], rAlgImbProjC8_100_[setIter + 5]);
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHemCF_[setIter + 5], rAlgImbHemC0_1_[setIter + 5], rAlgImbHemC1_2_[setIter + 5], rAlgImbHemC2_4_[setIter + 5], rAlgImbHemC4_8_[setIter + 5], rAlgImbHemC8_100_[setIter + 5]);

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjACF_[setIter + 5], tempHold, rAlgImbProjAC0_1_[setIter + 5], rAlgImbProjAC1_2_[setIter + 5], rAlgImbProjAC2_4_[setIter + 5], rAlgImbProjAC4_8_[setIter + 5], rAlgImbProjAC8_100_[setIter + 5]);
		}
		else{
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProjNCF_[setIter + 5], rAlgImbPerpNCF_[setIter + 5], rAlgImbProjNC0_1_[setIter + 5], rAlgImbProjNC1_2_[setIter + 5], rAlgImbProjNC2_4_[setIter + 5], rAlgImbProjNC4_8_[setIter + 5], rAlgImbProjNC8_100_[setIter + 5]);		
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHemNCF_[setIter + 5], rAlgImbHemNC0_1_[setIter + 5], rAlgImbHemNC1_2_[setIter + 5], rAlgImbHemNC2_4_[setIter + 5], rAlgImbHemNC4_8_[setIter + 5], rAlgImbHemNC8_100_[setIter + 5]);

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjANCF_[setIter + 5], tempHold, rAlgImbProjANC0_1_[setIter + 5], rAlgImbProjANC1_2_[setIter + 5], rAlgImbProjANC2_4_[setIter + 5], rAlgImbProjANC4_8_[setIter + 5], rAlgImbProjANC8_100_[setIter + 5]);		
		}

		if(tempLeadR[setIter] < .20 || tempSubLeadR[setIter] < .20){

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj1CF_[setIter + 5], tempHold, rAlgImbProj1C0_1_[setIter + 5], rAlgImbProj1C1_2_[setIter + 5], rAlgImbProj1C2_4_[setIter + 5], rAlgImbProj1C4_8_[setIter + 5], rAlgImbProj1C8_100_[setIter + 5]);
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHem1CF_[setIter + 5], rAlgImbHem1C0_1_[setIter + 5], rAlgImbHem1C1_2_[setIter + 5], rAlgImbHem1C2_4_[setIter + 5], rAlgImbHem1C4_8_[setIter + 5], rAlgImbHem1C8_100_[setIter + 5]);
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA1CF_[setIter + 5], tempHold, rAlgImbProjA1C0_1_[setIter + 5], rAlgImbProjA1C1_2_[setIter + 5], rAlgImbProjA1C2_4_[setIter + 5], rAlgImbProjA1C4_8_[setIter + 5], rAlgImbProjA1C8_100_[setIter + 5]);

		}
		else if(tempLeadR[setIter] < .40 || tempSubLeadR[setIter] < .40){

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj2CF_[setIter + 5], tempHold, rAlgImbProj2C0_1_[setIter + 5], rAlgImbProj2C1_2_[setIter + 5], rAlgImbProj2C2_4_[setIter + 5], rAlgImbProj2C4_8_[setIter + 5], rAlgImbProj2C8_100_[setIter + 5]);
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHem2CF_[setIter + 5], rAlgImbHem2C0_1_[setIter + 5], rAlgImbHem2C1_2_[setIter + 5], rAlgImbHem2C2_4_[setIter + 5], rAlgImbHem2C4_8_[setIter + 5], rAlgImbHem2C8_100_[setIter + 5]);
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA2CF_[setIter + 5], tempHold, rAlgImbProjA2C0_1_[setIter + 5], rAlgImbProjA2C1_2_[setIter + 5], rAlgImbProjA2C2_4_[setIter + 5], rAlgImbProjA2C4_8_[setIter + 5], rAlgImbProjA2C8_100_[setIter + 5]);

		}
		else if(tempLeadR[setIter] < .60 || tempSubLeadR[setIter] < .60){

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj3CF_[setIter + 5], tempHold, rAlgImbProj3C0_1_[setIter + 5], rAlgImbProj3C1_2_[setIter + 5], rAlgImbProj3C2_4_[setIter + 5], rAlgImbProj3C4_8_[setIter + 5], rAlgImbProj3C8_100_[setIter + 5]);
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHem3CF_[setIter + 5], rAlgImbHem3C0_1_[setIter + 5], rAlgImbHem3C1_2_[setIter + 5], rAlgImbHem3C2_4_[setIter + 5], rAlgImbHem3C4_8_[setIter + 5], rAlgImbHem3C8_100_[setIter + 5]);
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA3CF_[setIter + 5], tempHold, rAlgImbProjA3C0_1_[setIter + 5], rAlgImbProjA3C1_2_[setIter + 5], rAlgImbProjA3C2_4_[setIter + 5], rAlgImbProjA3C4_8_[setIter + 5], rAlgImbProjA3C8_100_[setIter + 5]);

		}
		else if(tempLeadR[setIter] < .80 || tempSubLeadR[setIter] < .80){

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj4CF_[setIter + 5], tempHold, rAlgImbProj4C0_1_[setIter + 5], rAlgImbProj4C1_2_[setIter + 5], rAlgImbProj4C2_4_[setIter + 5], rAlgImbProj4C4_8_[setIter + 5], rAlgImbProj4C8_100_[setIter + 5]);
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHem4CF_[setIter + 5], rAlgImbHem4C0_1_[setIter + 5], rAlgImbHem4C1_2_[setIter + 5], rAlgImbHem4C2_4_[setIter + 5], rAlgImbHem4C4_8_[setIter + 5], rAlgImbHem4C8_100_[setIter + 5]);
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA4CF_[setIter + 5], tempHold, rAlgImbProjA4C0_1_[setIter + 5], rAlgImbProjA4C1_2_[setIter + 5], rAlgImbProjA4C2_4_[setIter + 5], rAlgImbProjA4C4_8_[setIter + 5], rAlgImbProjA4C8_100_[setIter + 5]);

		}
		else if(tempLeadR[setIter] < 1.0 || tempSubLeadR[setIter] < 1.0){

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj5CF_[setIter + 5], tempHold, rAlgImbProj5C0_1_[setIter + 5], rAlgImbProj5C1_2_[setIter + 5], rAlgImbProj5C2_4_[setIter + 5], rAlgImbProj5C4_8_[setIter + 5], rAlgImbProj5C8_100_[setIter + 5]);
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHem5CF_[setIter + 5], rAlgImbHem5C0_1_[setIter + 5], rAlgImbHem5C1_2_[setIter + 5], rAlgImbHem5C2_4_[setIter + 5], rAlgImbHem5C4_8_[setIter + 5], rAlgImbHem5C8_100_[setIter + 5]);
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA5CF_[setIter + 5], tempHold, rAlgImbProjA5C0_1_[setIter + 5], rAlgImbProjA5C1_2_[setIter + 5], rAlgImbProjA5C2_4_[setIter + 5], rAlgImbProjA5C4_8_[setIter + 5], rAlgImbProjA5C8_100_[setIter + 5]);

		}
		else if(tempLeadR[setIter] < 1.2 || tempSubLeadR[setIter] < 1.2){

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj6CF_[setIter + 5], tempHold, rAlgImbProj6C0_1_[setIter + 5], rAlgImbProj6C1_2_[setIter + 5], rAlgImbProj6C2_4_[setIter + 5], rAlgImbProj6C4_8_[setIter + 5], rAlgImbProj6C8_100_[setIter + 5]);
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHem6CF_[setIter + 5], rAlgImbHem6C0_1_[setIter + 5], rAlgImbHem6C1_2_[setIter + 5], rAlgImbHem6C2_4_[setIter + 5], rAlgImbHem6C4_8_[setIter + 5], rAlgImbHem6C8_100_[setIter + 5]);
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA6CF_[setIter + 5], tempHold, rAlgImbProjA6C0_1_[setIter + 5], rAlgImbProjA6C1_2_[setIter + 5], rAlgImbProjA6C2_4_[setIter + 5], rAlgImbProjA6C4_8_[setIter + 5], rAlgImbProjA6C8_100_[setIter + 5]);

		}
		else if(tempLeadR[setIter] < 1.4 || tempSubLeadR[setIter] < 1.4){

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj7CF_[setIter + 5], tempHold, rAlgImbProj7C0_1_[setIter + 5], rAlgImbProj7C1_2_[setIter + 5], rAlgImbProj7C2_4_[setIter + 5], rAlgImbProj7C4_8_[setIter + 5], rAlgImbProj7C8_100_[setIter + 5]);
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHem7CF_[setIter + 5], rAlgImbHem7C0_1_[setIter + 5], rAlgImbHem7C1_2_[setIter + 5], rAlgImbHem7C2_4_[setIter + 5], rAlgImbHem7C4_8_[setIter + 5], rAlgImbHem7C8_100_[setIter + 5]);
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA7CF_[setIter + 5], tempHold, rAlgImbProjA7C0_1_[setIter + 5], rAlgImbProjA7C1_2_[setIter + 5], rAlgImbProjA7C2_4_[setIter + 5], rAlgImbProjA7C4_8_[setIter + 5], rAlgImbProjA7C8_100_[setIter + 5]);

		}
		else if(tempLeadR[setIter] < 1.6 || tempSubLeadR[setIter] < 1.6){

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj8CF_[setIter + 5], tempHold, rAlgImbProj8C0_1_[setIter + 5], rAlgImbProj8C1_2_[setIter + 5], rAlgImbProj8C2_4_[setIter + 5], rAlgImbProj8C4_8_[setIter + 5], rAlgImbProj8C8_100_[setIter + 5]);
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHem8CF_[setIter + 5], rAlgImbHem8C0_1_[setIter + 5], rAlgImbHem8C1_2_[setIter + 5], rAlgImbHem8C2_4_[setIter + 5], rAlgImbHem8C4_8_[setIter + 5], rAlgImbHem8C8_100_[setIter + 5]);
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA8CF_[setIter + 5], tempHold, rAlgImbProjA8C0_1_[setIter + 5], rAlgImbProjA8C1_2_[setIter + 5], rAlgImbProjA8C2_4_[setIter + 5], rAlgImbProjA8C4_8_[setIter + 5], rAlgImbProjA8C8_100_[setIter + 5]);

		}
		else if(tempLeadR[setIter] < 1.8 || tempSubLeadR[setIter] < 1.8){

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj9CF_[setIter + 5], tempHold, rAlgImbProj9C0_1_[setIter + 5], rAlgImbProj9C1_2_[setIter + 5], rAlgImbProj9C2_4_[setIter + 5], rAlgImbProj9C4_8_[setIter + 5], rAlgImbProj9C8_100_[setIter + 5]);
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHem9CF_[setIter + 5], rAlgImbHem9C0_1_[setIter + 5], rAlgImbHem9C1_2_[setIter + 5], rAlgImbHem9C2_4_[setIter + 5], rAlgImbHem9C4_8_[setIter + 5], rAlgImbHem9C8_100_[setIter + 5]);
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA9CF_[setIter + 5], tempHold, rAlgImbProjA9C0_1_[setIter + 5], rAlgImbProjA9C1_2_[setIter + 5], rAlgImbProjA9C2_4_[setIter + 5], rAlgImbProjA9C4_8_[setIter + 5], rAlgImbProjA9C8_100_[setIter + 5]);
		  
		}
		else{

		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbProj10CF_[setIter + 5], tempHold, rAlgImbProj10C0_1_[setIter + 5], rAlgImbProj10C1_2_[setIter + 5], rAlgImbProj10C2_4_[setIter + 5], rAlgImbProj10C4_8_[setIter + 5], rAlgImbProj10C8_100_[setIter + 5]);
		  getPtHem(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgLeadJtPhi_[setIter], rAlgImbHem10CF_[setIter + 5], rAlgImbHem10C0_1_[setIter + 5], rAlgImbHem10C1_2_[setIter + 5], rAlgImbHem10C2_4_[setIter + 5], rAlgImbHem10C4_8_[setIter + 5], rAlgImbHem10C8_100_[setIter + 5]);
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA10CF_[setIter + 5], tempHold, rAlgImbProjA10C0_1_[setIter + 5], rAlgImbProjA10C1_2_[setIter + 5], rAlgImbProjA10C2_4_[setIter + 5], rAlgImbProjA10C4_8_[setIter + 5], rAlgImbProjA10C8_100_[setIter + 5]);
		  
		}
		    
	      }

	    }
	  }
	  
	}
	break;
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

	if(eventSet_[VsPF])
          genPtVsPF_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], AlgLeadJtPhi_[VsPF]));
	else
	  genPtVsPF_[nGen_] = 0;

	if(eventSet_[VsCalo])
          genPtVsCalo_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], AlgLeadJtPhi_[VsCalo]));
	else
	  genPtVsCalo_[nGen_] = 0;

	
	if(eventSet_[T])
	  genLeadDelPhi_[nGen_] = getAbsDphi(AlgLeadJtPhi_[T], genCollection.phi[genEntry]);
	else
	  genLeadDelPhi_[nGen_] = -10;
      
       
	for(Int_t setIter = 0; setIter < 5; setIter++){
	  if(eventSet_[setIter]){
	    getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProjF_[setIter], gAlgImbPerpF_[setIter], gAlgImbProj0_1_[setIter], gAlgImbProj1_2_[setIter], gAlgImbProj2_4_[setIter], gAlgImbProj4_8_[setIter], gAlgImbProj8_100_[setIter]);

	    getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHemF_[setIter], gAlgImbHem0_1_[setIter], gAlgImbHem1_2_[setIter], gAlgImbHem2_4_[setIter], gAlgImbHem4_8_[setIter], gAlgImbHem8_100_[setIter]);

	    Float_t tempHold = 0;
	    
	    getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjAF_[setIter], tempHold, gAlgImbProjA0_1_[setIter], gAlgImbProjA1_2_[setIter], gAlgImbProjA2_4_[setIter], gAlgImbProjA4_8_[setIter], gAlgImbProjA8_100_[setIter]);

	  
	    Float_t tempLeadDelR = getDR(genEta_[nGen_], genPhi_[nGen_], AlgLeadJtEta_[setIter], AlgLeadJtPhi_[setIter]);
	    Float_t tempSubLeadDelR = getDR(genEta_[nGen_], genPhi_[nGen_], AlgSubLeadJtEta_[setIter], AlgSubLeadJtPhi_[setIter]);

	    if(tempLeadDelR > 0 && tempSubLeadDelR > 0){
	      if(tempLeadDelR < .8 || tempSubLeadDelR < .8){
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProjCF_[setIter], gAlgImbPerpCF_[setIter], gAlgImbProjC0_1_[setIter], gAlgImbProjC1_2_[setIter], gAlgImbProjC2_4_[setIter], gAlgImbProjC4_8_[setIter], gAlgImbProjC8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHemCF_[setIter], gAlgImbHemC0_1_[setIter], gAlgImbHemC1_2_[setIter], gAlgImbHemC2_4_[setIter], gAlgImbHemC4_8_[setIter], gAlgImbHemC8_100_[setIter]);

		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjACF_[setIter], tempHold, gAlgImbProjAC0_1_[setIter], gAlgImbProjAC1_2_[setIter], gAlgImbProjAC2_4_[setIter], gAlgImbProjAC4_8_[setIter], gAlgImbProjAC8_100_[setIter]);

		if(genCollection.sube[genEntry] == 0){
		  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProjSigCF_[setIter], gAlgImbPerpSigCF_[setIter], gAlgImbProjSigC0_1_[setIter], gAlgImbProjSigC1_2_[setIter], gAlgImbProjSigC2_4_[setIter], gAlgImbProjSigC4_8_[setIter], gAlgImbProjSigC8_100_[setIter]);
		  getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHemSigCF_[setIter], gAlgImbHemSigC0_1_[setIter], gAlgImbHemSigC1_2_[setIter], gAlgImbHemSigC2_4_[setIter], gAlgImbHemSigC4_8_[setIter], gAlgImbHemSigC8_100_[setIter]);

		  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjSigACF_[setIter], tempHold, gAlgImbProjSigAC0_1_[setIter], gAlgImbProjSigAC1_2_[setIter], gAlgImbProjSigAC2_4_[setIter], gAlgImbProjSigAC4_8_[setIter], gAlgImbProjSigAC8_100_[setIter]);
		}
		else{
		  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProjUECF_[setIter], gAlgImbPerpUECF_[setIter], gAlgImbProjUEC0_1_[setIter], gAlgImbProjUEC1_2_[setIter], gAlgImbProjUEC2_4_[setIter], gAlgImbProjUEC4_8_[setIter], gAlgImbProjUEC8_100_[setIter]);
		  getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHemUECF_[setIter], gAlgImbHemUEC0_1_[setIter], gAlgImbHemUEC1_2_[setIter], gAlgImbHemUEC2_4_[setIter], gAlgImbHemUEC4_8_[setIter], gAlgImbHemUEC8_100_[setIter]);

		  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjUEACF_[setIter], tempHold, gAlgImbProjUEAC0_1_[setIter], gAlgImbProjUEAC1_2_[setIter], gAlgImbProjUEAC2_4_[setIter], gAlgImbProjUEAC4_8_[setIter], gAlgImbProjUEAC8_100_[setIter]);
		}

	      }
	      else{
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProjNCF_[setIter], gAlgImbPerpNCF_[setIter], gAlgImbProjNC0_1_[setIter], gAlgImbProjNC1_2_[setIter], gAlgImbProjNC2_4_[setIter], gAlgImbProjNC4_8_[setIter], gAlgImbProjNC8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHemNCF_[setIter], gAlgImbHemNC0_1_[setIter], gAlgImbHemNC1_2_[setIter], gAlgImbHemNC2_4_[setIter], gAlgImbHemNC4_8_[setIter], gAlgImbHemNC8_100_[setIter]);
	
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjANCF_[setIter], tempHold, gAlgImbProjANC0_1_[setIter], gAlgImbProjANC1_2_[setIter], gAlgImbProjANC2_4_[setIter], gAlgImbProjANC4_8_[setIter], gAlgImbProjANC8_100_[setIter]);

		if(genCollection.sube[genEntry] == 0){
		  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProjSigNCF_[setIter], gAlgImbPerpSigNCF_[setIter], gAlgImbProjSigNC0_1_[setIter], gAlgImbProjSigNC1_2_[setIter], gAlgImbProjSigNC2_4_[setIter], gAlgImbProjSigNC4_8_[setIter], gAlgImbProjSigNC8_100_[setIter]);
		  getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHemSigNCF_[setIter], gAlgImbHemSigNC0_1_[setIter], gAlgImbHemSigNC1_2_[setIter], gAlgImbHemSigNC2_4_[setIter], gAlgImbHemSigNC4_8_[setIter], gAlgImbHemSigNC8_100_[setIter]);

		  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjSigANCF_[setIter], tempHold, gAlgImbProjSigANC0_1_[setIter], gAlgImbProjSigANC1_2_[setIter], gAlgImbProjSigANC2_4_[setIter], gAlgImbProjSigANC4_8_[setIter], gAlgImbProjSigANC8_100_[setIter]);
		}
		else{
		  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProjUENCF_[setIter], gAlgImbPerpUENCF_[setIter], gAlgImbProjUENC0_1_[setIter], gAlgImbProjUENC1_2_[setIter], gAlgImbProjUENC2_4_[setIter], gAlgImbProjUENC4_8_[setIter], gAlgImbProjUENC8_100_[setIter]);
		  getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHemUENCF_[setIter], gAlgImbHemUENC0_1_[setIter], gAlgImbHemUENC1_2_[setIter], gAlgImbHemUENC2_4_[setIter], gAlgImbHemUENC4_8_[setIter], gAlgImbHemUENC8_100_[setIter]);

		  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjUEANCF_[setIter], tempHold, gAlgImbProjUEANC0_1_[setIter], gAlgImbProjUEANC1_2_[setIter], gAlgImbProjUEANC2_4_[setIter], gAlgImbProjUEANC4_8_[setIter], gAlgImbProjUEANC8_100_[setIter]);
		}

	      }	      
	      

	      if(tempLeadDelR < .20 || tempSubLeadDelR < .20){

		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj1CF_[setIter], tempHold, gAlgImbProj1C0_1_[setIter], gAlgImbProj1C1_2_[setIter], gAlgImbProj1C2_4_[setIter], gAlgImbProj1C4_8_[setIter], gAlgImbProj1C8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHem1CF_[setIter], gAlgImbHem1C0_1_[setIter], gAlgImbHem1C1_2_[setIter], gAlgImbHem1C2_4_[setIter], gAlgImbHem1C4_8_[setIter], gAlgImbHem1C8_100_[setIter]);
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA1CF_[setIter], tempHold, gAlgImbProjA1C0_1_[setIter], gAlgImbProjA1C1_2_[setIter], gAlgImbProjA1C2_4_[setIter], gAlgImbProjA1C4_8_[setIter], gAlgImbProjA1C8_100_[setIter]);

	      }
	      else if(tempLeadDelR < .40 || tempSubLeadDelR < .40){

		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj2CF_[setIter], tempHold, gAlgImbProj2C0_1_[setIter], gAlgImbProj2C1_2_[setIter], gAlgImbProj2C2_4_[setIter], gAlgImbProj2C4_8_[setIter], gAlgImbProj2C8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHem2CF_[setIter], gAlgImbHem2C0_1_[setIter], gAlgImbHem2C1_2_[setIter], gAlgImbHem2C2_4_[setIter], gAlgImbHem2C4_8_[setIter], gAlgImbHem2C8_100_[setIter]);
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA2CF_[setIter], tempHold, gAlgImbProjA2C0_1_[setIter], gAlgImbProjA2C1_2_[setIter], gAlgImbProjA2C2_4_[setIter], gAlgImbProjA2C4_8_[setIter], gAlgImbProjA2C8_100_[setIter]);

	      }
	      else if(tempLeadDelR < .60 || tempSubLeadDelR < .60){

		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj3CF_[setIter], tempHold, gAlgImbProj3C0_1_[setIter], gAlgImbProj3C1_2_[setIter], gAlgImbProj3C2_4_[setIter], gAlgImbProj3C4_8_[setIter], gAlgImbProj3C8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHem3CF_[setIter], gAlgImbHem3C0_1_[setIter], gAlgImbHem3C1_2_[setIter], gAlgImbHem3C2_4_[setIter], gAlgImbHem3C4_8_[setIter], gAlgImbHem3C8_100_[setIter]);
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA3CF_[setIter], tempHold, gAlgImbProjA3C0_1_[setIter], gAlgImbProjA3C1_2_[setIter], gAlgImbProjA3C2_4_[setIter], gAlgImbProjA3C4_8_[setIter], gAlgImbProjA3C8_100_[setIter]);

	      }
	      else if(tempLeadDelR < .80 || tempSubLeadDelR < .80){

		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj4CF_[setIter], tempHold, gAlgImbProj4C0_1_[setIter], gAlgImbProj4C1_2_[setIter], gAlgImbProj4C2_4_[setIter], gAlgImbProj4C4_8_[setIter], gAlgImbProj4C8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHem4CF_[setIter], gAlgImbHem4C0_1_[setIter], gAlgImbHem4C1_2_[setIter], gAlgImbHem4C2_4_[setIter], gAlgImbHem4C4_8_[setIter], gAlgImbHem4C8_100_[setIter]);
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA4CF_[setIter], tempHold, gAlgImbProjA4C0_1_[setIter], gAlgImbProjA4C1_2_[setIter], gAlgImbProjA4C2_4_[setIter], gAlgImbProjA4C4_8_[setIter], gAlgImbProjA4C8_100_[setIter]);

	      }
	      else if(tempLeadDelR < 1.0 || tempSubLeadDelR < 1.0){

		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj5CF_[setIter], tempHold, gAlgImbProj5C0_1_[setIter], gAlgImbProj5C1_2_[setIter], gAlgImbProj5C2_4_[setIter], gAlgImbProj5C4_8_[setIter], gAlgImbProj5C8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHem5CF_[setIter], gAlgImbHem5C0_1_[setIter], gAlgImbHem5C1_2_[setIter], gAlgImbHem5C2_4_[setIter], gAlgImbHem5C4_8_[setIter], gAlgImbHem5C8_100_[setIter]);
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA5CF_[setIter], tempHold, gAlgImbProjA5C0_1_[setIter], gAlgImbProjA5C1_2_[setIter], gAlgImbProjA5C2_4_[setIter], gAlgImbProjA5C4_8_[setIter], gAlgImbProjA5C8_100_[setIter]);

	      }
	      else if(tempLeadDelR < 1.2 || tempSubLeadDelR < 1.2){

		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj6CF_[setIter], tempHold, gAlgImbProj6C0_1_[setIter], gAlgImbProj6C1_2_[setIter], gAlgImbProj6C2_4_[setIter], gAlgImbProj6C4_8_[setIter], gAlgImbProj6C8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHem6CF_[setIter], gAlgImbHem6C0_1_[setIter], gAlgImbHem6C1_2_[setIter], gAlgImbHem6C2_4_[setIter], gAlgImbHem6C4_8_[setIter], gAlgImbHem6C8_100_[setIter]);
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA6CF_[setIter], tempHold, gAlgImbProjA6C0_1_[setIter], gAlgImbProjA6C1_2_[setIter], gAlgImbProjA6C2_4_[setIter], gAlgImbProjA6C4_8_[setIter], gAlgImbProjA6C8_100_[setIter]);

	      }
	      else if(tempLeadDelR < 1.4 || tempSubLeadDelR < 1.4){

		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj7CF_[setIter], tempHold, gAlgImbProj7C0_1_[setIter], gAlgImbProj7C1_2_[setIter], gAlgImbProj7C2_4_[setIter], gAlgImbProj7C4_8_[setIter], gAlgImbProj7C8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHem7CF_[setIter], gAlgImbHem7C0_1_[setIter], gAlgImbHem7C1_2_[setIter], gAlgImbHem7C2_4_[setIter], gAlgImbHem7C4_8_[setIter], gAlgImbHem7C8_100_[setIter]);
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA7CF_[setIter], tempHold, gAlgImbProjA7C0_1_[setIter], gAlgImbProjA7C1_2_[setIter], gAlgImbProjA7C2_4_[setIter], gAlgImbProjA7C4_8_[setIter], gAlgImbProjA7C8_100_[setIter]);

	      }
	      else if(tempLeadDelR < 1.6 || tempSubLeadDelR < 1.6){

		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj8CF_[setIter], tempHold, gAlgImbProj8C0_1_[setIter], gAlgImbProj8C1_2_[setIter], gAlgImbProj8C2_4_[setIter], gAlgImbProj8C4_8_[setIter], gAlgImbProj8C8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHem8CF_[setIter], gAlgImbHem8C0_1_[setIter], gAlgImbHem8C1_2_[setIter], gAlgImbHem8C2_4_[setIter], gAlgImbHem8C4_8_[setIter], gAlgImbHem8C8_100_[setIter]);
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA8CF_[setIter], tempHold, gAlgImbProjA8C0_1_[setIter], gAlgImbProjA8C1_2_[setIter], gAlgImbProjA8C2_4_[setIter], gAlgImbProjA8C4_8_[setIter], gAlgImbProjA8C8_100_[setIter]);

	      }
	      else if(tempLeadDelR < 1.8 || tempSubLeadDelR < 1.8){

		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj9CF_[setIter], tempHold, gAlgImbProj9C0_1_[setIter], gAlgImbProj9C1_2_[setIter], gAlgImbProj9C2_4_[setIter], gAlgImbProj9C4_8_[setIter], gAlgImbProj9C8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHem9CF_[setIter], gAlgImbHem9C0_1_[setIter], gAlgImbHem9C1_2_[setIter], gAlgImbHem9C2_4_[setIter], gAlgImbHem9C4_8_[setIter], gAlgImbHem9C8_100_[setIter]);
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA9CF_[setIter], tempHold, gAlgImbProjA9C0_1_[setIter], gAlgImbProjA9C1_2_[setIter], gAlgImbProjA9C2_4_[setIter], gAlgImbProjA9C4_8_[setIter], gAlgImbProjA9C8_100_[setIter]);

	      }
	      else{

		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbProj10CF_[setIter], tempHold, gAlgImbProj10C0_1_[setIter], gAlgImbProj10C1_2_[setIter], gAlgImbProj10C2_4_[setIter], gAlgImbProj10C4_8_[setIter], gAlgImbProj10C8_100_[setIter]);
		getPtHem(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgLeadJtPhi_[setIter], gAlgImbHem10CF_[setIter], gAlgImbHem10C0_1_[setIter], gAlgImbHem10C1_2_[setIter], gAlgImbHem10C2_4_[setIter], gAlgImbHem10C4_8_[setIter], gAlgImbHem10C8_100_[setIter]);
		getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA10CF_[setIter], tempHold, gAlgImbProjA10C0_1_[setIter], gAlgImbProjA10C1_2_[setIter], gAlgImbProjA10C2_4_[setIter], gAlgImbProjA10C4_8_[setIter], gAlgImbProjA10C8_100_[setIter]);

	      }

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
  tempTot = tempTot - vzCut;
  std::cout << "vzCut: " << tempTot << std::endl;

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

  delete c;
  delete outFile;
  delete centHistFile_p;

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

  rStatus = makeDiJetTTree(argv[1], sampleType(atoi(argv[2])), argv[3], atoi(argv[4]));

  return rStatus;
}
