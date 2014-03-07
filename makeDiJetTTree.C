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
const Float_t jtDelPhiCut = 2.*(TMath::Pi())/3.;
const Float_t jtEtaCut = 1.6; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed

collisionType getCType(sampleType sType);

void getLeadJt(Float_t& leadJtPt, Int_t& lPtCut, Float_t& subLeadJtPt, Int_t& sPtCut, Float_t& leadJtPhi, Float_t& subLeadJtPhi, Int_t& dPhiCut, Float_t& leadJtEta, Float_t& subLeadJtEta, Int_t& etaCut, Float_t& jtDelPhi, Float_t& jtAsymm, Bool_t& setPass, Jets jtCollection)
{
  Int_t leadJtIndex = -1;
  Int_t subLeadJtIndex = -1;
  leadJtPt = subLeadJtPtCut;
  subLeadJtPt = subLeadJtPtCut;

  for(Int_t jtEntry = 0; jtEntry < jtCollection.nref; jtEntry++){
    if(jtCollection.jtpt[jtEntry] > leadJtPtCut && jtCollection.jtpt[jtEntry] > leadJtPt){
      subLeadJtIndex = leadJtIndex;
      subLeadJtPt = leadJtPt;
      leadJtIndex = jtEntry;
      leadJtPt = jtCollection.jtpt[jtEntry];
    }
    else if(jtCollection.jtpt[jtEntry] > subLeadJtPt){
      subLeadJtIndex = jtEntry;
      subLeadJtPt = jtCollection.jtpt[jtEntry];
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


int makeDiJetTTree(string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_CFMSKIM.root")
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

  std::cout << nentries << std::endl;

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

  Int_t VsCaloLeadJtPtCut = 0;
  Int_t VsCaloSubLeadJtPtCut = 0;
  Int_t VsCaloDelPhiCut = 0;
  Int_t VsCaloJtEtaCut = 0;

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

  for(Long64_t jentry = 0; jentry < 20000 /*nentries*/; jentry++){
    c->GetEntry(jentry);

    totEv++;

    if(jentry%1000 == 0)
      std::cout << jentry << std::endl;

    if(jentry > 32800)
      std::cout << jentry << std::endl;

    if(!c->selectEvent()){
      selectCut++;
      continue;
    }

    InitJetVar(montecarlo);

    //particle flow

    Jets PFJtCollection = c->akPu3PF;
    
    getLeadJt(PFLeadJtPt_, PFLeadJtPtCut, PFSubLeadJtPt_, PFSubLeadJtPtCut, PFLeadJtPhi_, PFSubLeadJtPhi_, PFDelPhiCut, PFLeadJtEta_, PFSubLeadJtEta_, PFJtEtaCut, PFJtDelPhi_, PFJtAsymm_, recoPFSet_, PFJtCollection);
    
    //Vs PF jet

    Jets VsPFJtCollection = c->akVs3PF;
    
    getLeadJt(VsPFLeadJtPt_, VsPFLeadJtPtCut, VsPFSubLeadJtPt_, VsPFSubLeadJtPtCut, VsPFLeadJtPhi_, VsPFSubLeadJtPhi_, VsPFDelPhiCut, VsPFLeadJtEta_, VsPFSubLeadJtEta_, VsPFJtEtaCut, VsPFJtDelPhi_, VsPFJtAsymm_, recoVsPFSet_, VsPFJtCollection);
    
    //calo jet

    Jets CaloJtCollection = c->akPu3Calo;
    
    getLeadJt(CaloLeadJtPt_, CaloLeadJtPtCut, CaloSubLeadJtPt_, CaloSubLeadJtPtCut, CaloLeadJtPhi_, CaloSubLeadJtPhi_, CaloDelPhiCut, CaloLeadJtEta_, CaloSubLeadJtEta_, CaloJtEtaCut, CaloJtDelPhi_, CaloJtAsymm_, recoCaloSet_, CaloJtCollection);
    
    //Vs calo jet

    Jets VsCaloJtCollection = c->akVs3Calo;
    
    getLeadJt(VsCaloLeadJtPt_, VsCaloLeadJtPtCut, VsCaloSubLeadJtPt_, VsCaloSubLeadJtPtCut, VsCaloLeadJtPhi_, VsCaloSubLeadJtPhi_, VsCaloDelPhiCut, VsCaloLeadJtEta_, VsCaloSubLeadJtEta_, VsCaloJtEtaCut, VsCaloJtDelPhi_, VsCaloJtAsymm_, recoVsCaloSet_, VsCaloJtCollection);
    
    //truth, doesn't work w/ getLeadJt because truth doesnt get its own tree

    if(montecarlo){

      Int_t leadJtIndex = -1;
      Int_t subLeadJtIndex = -1;
      TLeadJtPt_ = subLeadJtPtCut;
      TSubLeadJtPt_ = subLeadJtPtCut;
      for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.ngen; jtEntry++){
	if(c->akPu3PF.genpt[jtEntry] > leadJtPtCut && c->akPu3PF.genpt[jtEntry] > TLeadJtPt_){
	  subLeadJtIndex = leadJtIndex;
	  TSubLeadJtPt_ = TLeadJtPt_;
	  leadJtIndex = jtEntry;
	  TLeadJtPt_ = c->akPu3PF.genpt[jtEntry];
	}
	else if(c->akPu3PF.genpt[jtEntry] > TSubLeadJtPt_){
	  subLeadJtIndex = jtEntry;
	  TSubLeadJtPt_ = c->akPu3PF.genpt[jtEntry];
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
      else if(getAbsDphi(c->akPu3PF.genphi[leadJtIndex], c->akPu3PF.genphi[subLeadJtIndex]) < jtDelPhiCut){
	TDelPhiCut++;
      }
      else if(TMath::Abs(c->akPu3PF.geneta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akPu3PF.geneta[subLeadJtIndex]) > jtEtaCut){
	TJtEtaCut++;
      }
      else{
	TLeadJtPhi_ = c->akPu3PF.genphi[leadJtIndex];
	TSubLeadJtPhi_ = c->akPu3PF.genphi[subLeadJtIndex];
	TLeadJtEta_ = c->akPu3PF.geneta[leadJtIndex];
	TSubLeadJtEta_ = c->akPu3PF.geneta[subLeadJtIndex];
	
	TJtDelPhi_ = getAbsDphi(c->akPu3PF.genphi[leadJtIndex], c->akPu3PF.genphi[subLeadJtIndex]);
	TJtAsymm_ = (TLeadJtPt_ - TSubLeadJtPt_)/(TLeadJtPt_ + TSubLeadJtPt_);
	
	truthSet_ = true;
      }
    }
    
    if(truthSet_ == false && recoPFSet_ == false && recoCaloSet_ == false && recoVsPFSet_ == false && recoVsCaloSet_ == false)
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
      if(truthSet_)	TTotTrk++;
      if(recoPFSet_)	PFTotTrk++;
      if(recoCaloSet_)	CaloTotTrk++;

      if(TMath::Abs(trkCollection.trkEta[trkEntry]) > 2.4){
	if(truthSet_)	  TTrkEtaCut++;
	if(recoPFSet_)	  PFTrkEtaCut++;
	if(recoCaloSet_)	  CaloTrkEtaCut++;

	continue;
      }
      
      if(trkCollection.trkPt[trkEntry] < 0.5){
	if(truthSet_)	  TTrkPtCut++;
	if(recoPFSet_)	  PFTrkPtCut++;
	if(recoCaloSet_)	  CaloTrkPtCut++;

	continue;
      }

      if(!trkCollection.highPurity[trkEntry]){ //Note highPuritySetWithPV seems to be wrong cut, creates diff. bet. truth and reco
	if(truthSet_)	  TPurityCut++;
	if(recoPFSet_)	  PFPurityCut++;
	if(recoCaloSet_)	  CaloPurityCut++;

	continue;
      }

      if(TMath::Abs(trkCollection.trkDz1[trkEntry]/trkCollection.trkDzError1[trkEntry]) > 3){
	if(truthSet_)	  TTrkDzCut++;
	if(recoPFSet_)	  PFTrkDzCut++;
	if(recoCaloSet_)	  CaloTrkDzCut++;

	continue;
      }

      if(TMath::Abs(trkCollection.trkDxy1[trkEntry]/trkCollection.trkDxyError1[trkEntry]) > 3){
	if(truthSet_)	  TTrkDxyCut++;
	if(recoPFSet_)	  PFTrkDxyCut++;
	if(recoCaloSet_)	  CaloTrkDxyCut++;

	continue;
      }

      if(trkCollection.trkPtError[trkEntry]/trkCollection.trkPt[trkEntry] > 0.1){
	if(truthSet_)	  TTrkPtErrorCut++;
	if(recoPFSet_)	  PFTrkPtErrorCut++;
	if(recoCaloSet_)	  CaloTrkPtErrorCut++;

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

      trkPtCorrPF_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPtCorrCalo_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPtCorrT_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPtCorrVsPF_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPtCorrVsCalo_[nTrk_] = trkCollection.trkPt[trkEntry];
      
      //Grab proj. Pt Spectra For Tracks in each Event Subset

      if(recoPFSet_){
	trkPtPF_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], PFLeadJtPhi_));
	trkRLeadPF_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], PFLeadJtEta_, PFLeadJtPhi_);
	trkRSubLeadPF_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], PFSubLeadJtEta_, PFSubLeadJtPhi_);
        trkPFLeadDelPhi_[nTrk_] = getAbsDphi(PFLeadJtPhi_, trkCollection.trkPhi[trkEntry]);
      }
      else{
	trkPtPF_[nTrk_] = 0;
	trkRLeadPF_[nTrk_] = -1;
	trkRSubLeadPF_[nTrk_] = -1;
        trkPFLeadDelPhi_[nTrk_] = -10;
      }

      if(recoCaloSet_){
	trkPtCalo_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_));
	trkRLeadCalo_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], CaloLeadJtEta_, CaloLeadJtPhi_);
      	trkRSubLeadCalo_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], CaloSubLeadJtEta_, CaloSubLeadJtPhi_);
        trkCaloLeadDelPhi_[nTrk_] = getAbsDphi(CaloLeadJtPhi_, trkCollection.trkPhi[trkEntry]);
      }
      else{
	trkPtCalo_[nTrk_] = 0;
	trkRLeadCalo_[nTrk_] = -1;
	trkRSubLeadCalo_[nTrk_] = -1;
        trkCaloLeadDelPhi_[nTrk_] = -10;
      }

      if(recoVsPFSet_){
	trkPtVsPF_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], VsPFLeadJtPhi_));
	trkRLeadVsPF_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], VsPFLeadJtEta_, VsPFLeadJtPhi_);
	trkRSubLeadVsPF_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], VsPFSubLeadJtEta_, VsPFSubLeadJtPhi_);
        trkVsPFLeadDelPhi_[nTrk_] = getAbsDphi(VsPFLeadJtPhi_, trkCollection.trkPhi[trkEntry]);
      }
      else{
	trkPtVsPF_[nTrk_] = 0;
	trkRLeadVsPF_[nTrk_] = -1;
	trkRSubLeadVsPF_[nTrk_] = -1;
        trkVsPFLeadDelPhi_[nTrk_] = -10;
      }

      if(recoVsCaloSet_){
	trkPtVsCalo_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], VsCaloLeadJtPhi_));
	trkRLeadVsCalo_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], VsCaloLeadJtEta_, VsCaloLeadJtPhi_);
	trkRSubLeadVsCalo_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], VsCaloSubLeadJtEta_, VsCaloSubLeadJtPhi_);
        trkVsCaloLeadDelPhi_[nTrk_] = getAbsDphi(VsCaloLeadJtPhi_, trkCollection.trkPhi[trkEntry]);
      }
      else{
	trkPtVsCalo_[nTrk_] = 0;
	trkRLeadVsCalo_[nTrk_] = -1;
	trkRSubLeadVsCalo_[nTrk_] = -1;
        trkVsCaloLeadDelPhi_[nTrk_] = -10;
      }

      if(montecarlo && truthSet_){
	trkPtT_[nTrk_] = -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], TLeadJtPhi_));
	trkRLeadT_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], TLeadJtEta_, TLeadJtPhi_);
	trkRSubLeadT_[nTrk_] = getDR(trkEta_[nTrk_], trkPhi_[nTrk_], TSubLeadJtEta_, TSubLeadJtPhi_);
        trkTLeadDelPhi_[nTrk_] = getAbsDphi(TLeadJtPhi_, trkCollection.trkPhi[trkEntry]);
      }
      else{
	trkPtT_[nTrk_] = 0;
	trkRLeadT_[nTrk_] = -1;
	trkRSubLeadT_[nTrk_] = -1;
        trkTLeadDelPhi_[nTrk_] = -10;
      }


      if(recoPFSet_){
	getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], PFLeadJtPhi_, rPFImbProjF_, rPFImbPerpF_, rPFImbProj0_1_, rPFImbProj1_2_, rPFImbProj2_4_, rPFImbProj4_8_, rPFImbProj8_100_);

	if(trkRLeadPF_[nTrk_] > 0 && trkRSubLeadPF_[nTrk_] > 0){
	  if(trkRLeadPF_[nTrk_] < .8 || trkRSubLeadPF_[nTrk_] < .8)
	    getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], PFLeadJtPhi_, rPFImbProjCF_, rPFImbPerpCF_, rPFImbProjC0_1_, rPFImbProjC1_2_, rPFImbProjC2_4_, rPFImbProjC4_8_, rPFImbProjC8_100_);
	  else
	    getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], PFLeadJtPhi_, rPFImbProjNCF_, rPFImbPerpNCF_, rPFImbProjNC0_1_, rPFImbProjNC1_2_, rPFImbProjNC2_4_, rPFImbProjNC4_8_, rPFImbProjNC8_100_);
	}
      }

      if(recoCaloSet_){
	getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_, rCaloImbProjF_, rCaloImbPerpF_, rCaloImbProj0_1_, rCaloImbProj1_2_, rCaloImbProj2_4_, rCaloImbProj4_8_, rCaloImbProj8_100_);

	if(trkRLeadCalo_[nTrk_] > 0 && trkRSubLeadCalo_[nTrk_] > 0){
	  if(trkRLeadCalo_[nTrk_] < .8 || trkRSubLeadCalo_[nTrk_] < .8)
	    getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_, rCaloImbProjCF_, rCaloImbPerpCF_, rCaloImbProjC0_1_, rCaloImbProjC1_2_, rCaloImbProjC2_4_, rCaloImbProjC4_8_, rCaloImbProjC8_100_);
	  else
	    getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], CaloLeadJtPhi_, rCaloImbProjNCF_, rCaloImbPerpNCF_, rCaloImbProjNC0_1_, rCaloImbProjNC1_2_, rCaloImbProjNC2_4_, rCaloImbProjNC4_8_, rCaloImbProjNC8_100_);
	}
      }

      if(recoVsPFSet_){
	getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], VsPFLeadJtPhi_, rVsPFImbProjF_, rVsPFImbPerpF_, rVsPFImbProj0_1_, rVsPFImbProj1_2_, rVsPFImbProj2_4_, rVsPFImbProj4_8_, rVsPFImbProj8_100_);

	if(trkRLeadVsPF_[nTrk_] > 0 && trkRSubLeadVsPF_[nTrk_] > 0){
	  if(trkRLeadVsPF_[nTrk_] < .8 || trkRSubLeadVsPF_[nTrk_] < .8)
	    getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], VsPFLeadJtPhi_, rVsPFImbProjCF_, rVsPFImbPerpCF_, rVsPFImbProjC0_1_, rVsPFImbProjC1_2_, rVsPFImbProjC2_4_, rVsPFImbProjC4_8_, rVsPFImbProjC8_100_);
	  else
	    getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], VsPFLeadJtPhi_, rVsPFImbProjNCF_, rVsPFImbPerpNCF_, rVsPFImbProjNC0_1_, rVsPFImbProjNC1_2_, rVsPFImbProjNC2_4_, rVsPFImbProjNC4_8_, rVsPFImbProjNC8_100_);
	}
      }

      if(recoVsCaloSet_){
	getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], VsCaloLeadJtPhi_, rVsCaloImbProjF_, rVsCaloImbPerpF_, rVsCaloImbProj0_1_, rVsCaloImbProj1_2_, rVsCaloImbProj2_4_, rVsCaloImbProj4_8_, rVsCaloImbProj8_100_);

	if(trkRLeadVsCalo_[nTrk_] > 0 && trkRSubLeadVsCalo_[nTrk_] > 0){
	  if(trkRLeadVsCalo_[nTrk_] < .8 || trkRSubLeadVsCalo_[nTrk_] < .8)
	    getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], VsCaloLeadJtPhi_, rVsCaloImbProjCF_, rVsCaloImbPerpCF_, rVsCaloImbProjC0_1_, rVsCaloImbProjC1_2_, rVsCaloImbProjC2_4_, rVsCaloImbProjC4_8_, rVsCaloImbProjC8_100_);
	  else
	    getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], VsCaloLeadJtPhi_, rVsCaloImbProjNCF_, rVsCaloImbPerpNCF_, rVsCaloImbProjNC0_1_, rVsCaloImbProjNC1_2_, rVsCaloImbProjNC2_4_, rVsCaloImbProjNC4_8_, rVsCaloImbProjNC8_100_);
	}
      }

      if(montecarlo && truthSet_){
	getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], TLeadJtPhi_, rTImbProjF_, rTImbPerpF_, rTImbProj0_1_, rTImbProj1_2_, rTImbProj2_4_, rTImbProj4_8_, rTImbProj8_100_);

	if(trkRLeadT_[nTrk_] > 0 && trkRSubLeadT_[nTrk_] > 0){
	  if(trkRLeadT_[nTrk_] < .8 || trkRSubLeadT_[nTrk_] < .8)
	    getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], TLeadJtPhi_, rTImbProjCF_, rTImbPerpCF_, rTImbProjC0_1_, rTImbProjC1_2_, rTImbProjC2_4_, rTImbProjC4_8_, rTImbProjC8_100_);
	  else
	    getPtProj(trkCollection.trkPt[trkEntry], trkCollection.trkPt[trkEntry], trkCollection.trkPhi[trkEntry], TLeadJtPhi_, rTImbProjNCF_, rTImbPerpNCF_, rTImbProjNC0_1_, rTImbProjNC1_2_, rTImbProjNC2_4_, rTImbProjNC4_8_, rTImbProjNC8_100_);
	}
      }
    
      nTrk_++;
      if(nTrk_ > MAXTRKS - 1){
	printf("ERROR: Trk arrays not large enough.\n");
	return(1);
      }
    }

    //Get trk rmin/corrections for 3 Jet subsets

    //Test hiBin first, less touches over trk for loop

    
    for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
      trkRMinPF_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], PFJtCollection);
      trkRMinCalo_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], CaloJtCollection);
      trkRMinVsPF_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], VsPFJtCollection);
      trkRMinVsCalo_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], VsCaloJtCollection);
    }

    if(montecarlo){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	trkRMinT_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], PFJtCollection, true);
      }
    }

    
    Int_t hiBinDiv[5] = {20, 40, 60, 100, 200};
    Int_t hiSetEff[15] = {0, 5, 10, 1, 6, 11, 2, 7, 12, 3, 8, 12, 4, 9, 12};
    Int_t hiSetFake[15] = {0, 5, 10, 1, 6, 11, 2, 7, 12, 3, 8, 12, 4, 9, 13};
  

    for(Int_t hiBinIter = 0; hiBinIter < 5; hiBinIter++){
      if(hiBin_ < hiBinDiv[hiBinIter]){
	for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	  Int_t ptPosEff = getPtBin(trkPt_[trkEntry], hiSetEff[hiBinIter*3], hiSetEff[hiBinIter*3 + 1], hiSetEff[hiBinIter*3 + 2], 13);
	  Int_t ptPosFake = getPtBin(trkPt_[trkEntry], hiSetFake[hiBinIter*3], hiSetFake[hiBinIter*3 + 1], hiSetFake[hiBinIter*3 + 2], 14);
	
	  trkPtCorrPF_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent_p[ptPosEff], PFphiEta_p[ptPosEff], PFpt_p[ptPosEff], PFdelR_p[ptPosEff], 3);
	  trkPtFactPF_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent_p[ptPosEff], PFphiEta_p[ptPosEff], PFpt_p[ptPosEff], PFdelR_p[ptPosEff], 3);

	  trkPtCorrCalo_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinCalo_[trkEntry], Calocent_p[ptPosEff], CalophiEta_p[ptPosEff], Calopt_p[ptPosEff], CalodelR_p[ptPosEff], 5);
	  trkPtFactCalo_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinCalo_[trkEntry], Calocent_p[ptPosEff], CalophiEta_p[ptPosEff], Calopt_p[ptPosEff], CalodelR_p[ptPosEff], 5);


	  trkPtFactVsCalo_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinVsCalo_[trkEntry], VsCalocent_p[ptPosEff], VsCalophiEta_p[ptPosEff], VsCalopt_p[ptPosEff], VsCalodelR_p[ptPosEff], 5);
	  trkPtCorrVsCalo_[trkEntry] = trkPt_[trkEntry]*(1 - factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinVsCalo_[trkEntry], FakeVsCalocent_p[ptPosFake], FakeVsCalophiEta_p[ptPosFake], FakeVsCalopt_p[ptPosFake], FakeVsCalodelR_p[ptPosFake], 5, false))/trkPtFactVsCalo_[trkEntry];

	}
	break;
      }
    }
    

    //Apply corrections to appropriate subsets

    if(recoPFSet_){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	getPtProj(trkPt_[trkEntry], trkPtCorrPF_[trkEntry], trkPhi_[trkEntry], PFLeadJtPhi_, rPFImbProjFCorr_, rPFImbPerpFCorr_, rPFImbProj0_1Corr_, rPFImbProj1_2Corr_, rPFImbProj2_4Corr_, rPFImbProj4_8Corr_, rPFImbProj8_100Corr_);

	//fix here

	if(trkRLeadPF_[trkEntry] > 0 && trkRSubLeadPF_[trkEntry] > 0){
	  if(trkRLeadPF_[trkEntry] < .8 || trkRSubLeadPF_[trkEntry] < .8)
	    getPtProj(trkPt_[trkEntry], trkPtCorrPF_[trkEntry], trkPhi_[trkEntry], PFLeadJtPhi_, rPFImbProjCFCorr_, rPFImbPerpCFCorr_, rPFImbProjC0_1Corr_, rPFImbProjC1_2Corr_, rPFImbProjC2_4Corr_, rPFImbProjC4_8Corr_, rPFImbProjC8_100Corr_);
	  else
	    getPtProj(trkPt_[trkEntry], trkPtCorrPF_[trkEntry], trkPhi_[trkEntry], PFLeadJtPhi_, rPFImbProjNCFCorr_, rPFImbPerpNCFCorr_, rPFImbProjNC0_1Corr_, rPFImbProjNC1_2Corr_, rPFImbProjNC2_4Corr_, rPFImbProjNC4_8Corr_, rPFImbProjNC8_100Corr_);
	}
      }
    }
    
    if(recoCaloSet_){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	getPtProj(trkPt_[trkEntry], trkPtCorrCalo_[trkEntry], trkPhi_[trkEntry], CaloLeadJtPhi_, rCaloImbProjFCorr_, rCaloImbPerpFCorr_, rCaloImbProj0_1Corr_, rCaloImbProj1_2Corr_, rCaloImbProj2_4Corr_, rCaloImbProj4_8Corr_, rCaloImbProj8_100Corr_);

	if(trkRLeadCalo_[trkEntry] > 0 && trkRSubLeadCalo_[trkEntry] > 0){
	  if(trkRLeadCalo_[trkEntry] < .8 || trkRSubLeadCalo_[trkEntry] < .8)
	    getPtProj(trkPt_[trkEntry], trkPtCorrCalo_[trkEntry], trkPhi_[trkEntry], CaloLeadJtPhi_, rCaloImbProjCFCorr_, rCaloImbPerpCFCorr_, rCaloImbProjC0_1Corr_, rCaloImbProjC1_2Corr_, rCaloImbProjC2_4Corr_, rCaloImbProjC4_8Corr_, rCaloImbProjC8_100Corr_);
	  else
	    getPtProj(trkPt_[trkEntry], trkPtCorrCalo_[trkEntry], trkPhi_[trkEntry], CaloLeadJtPhi_, rCaloImbProjNCFCorr_, rCaloImbPerpNCFCorr_, rCaloImbProjNC0_1Corr_, rCaloImbProjNC1_2Corr_, rCaloImbProjNC2_4Corr_, rCaloImbProjNC4_8Corr_, rCaloImbProjNC8_100Corr_);
	}
      }
    }    

    if(recoVsPFSet_){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	getPtProj(trkPt_[trkEntry], trkPtCorrVsPF_[trkEntry], trkPhi_[trkEntry], VsPFLeadJtPhi_, rVsPFImbProjFCorr_, rVsPFImbPerpFCorr_, rVsPFImbProj0_1Corr_, rVsPFImbProj1_2Corr_, rVsPFImbProj2_4Corr_, rVsPFImbProj4_8Corr_, rVsPFImbProj8_100Corr_);

	if(trkRLeadVsPF_[trkEntry] > 0 && trkRSubLeadVsPF_[trkEntry] > 0){
	  if(trkRLeadVsPF_[trkEntry] < .8 || trkRSubLeadVsPF_[trkEntry] < .8)
	    getPtProj(trkPt_[trkEntry], trkPtCorrVsPF_[trkEntry], trkPhi_[trkEntry], VsPFLeadJtPhi_, rVsPFImbProjCFCorr_, rVsPFImbPerpCFCorr_, rVsPFImbProjC0_1Corr_, rVsPFImbProjC1_2Corr_, rVsPFImbProjC2_4Corr_, rVsPFImbProjC4_8Corr_, rVsPFImbProjC8_100Corr_);
	  else
	    getPtProj(trkPt_[trkEntry], trkPtCorrVsPF_[trkEntry], trkPhi_[trkEntry], VsPFLeadJtPhi_, rVsPFImbProjNCFCorr_, rVsPFImbPerpNCFCorr_, rVsPFImbProjNC0_1Corr_, rVsPFImbProjNC1_2Corr_, rVsPFImbProjNC2_4Corr_, rVsPFImbProjNC4_8Corr_, rVsPFImbProjNC8_100Corr_);
	}
      }
    }

    if(recoVsCaloSet_){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	getPtProj(trkPt_[trkEntry], trkPtCorrVsCalo_[trkEntry], trkPhi_[trkEntry], VsCaloLeadJtPhi_, rVsCaloImbProjFCorr_, rVsCaloImbPerpFCorr_, rVsCaloImbProj0_1Corr_, rVsCaloImbProj1_2Corr_, rVsCaloImbProj2_4Corr_, rVsCaloImbProj4_8Corr_, rVsCaloImbProj8_100Corr_);

	if(trkRLeadVsCalo_[trkEntry] > 0 && trkRSubLeadVsCalo_[trkEntry] > 0){
	  if(trkRLeadVsCalo_[trkEntry] < .8 || trkRSubLeadVsCalo_[trkEntry] < .8)
	    getPtProj(trkPt_[trkEntry], trkPtCorrVsCalo_[trkEntry], trkPhi_[trkEntry], VsCaloLeadJtPhi_, rVsCaloImbProjCFCorr_, rVsCaloImbPerpCFCorr_, rVsCaloImbProjC0_1Corr_, rVsCaloImbProjC1_2Corr_, rVsCaloImbProjC2_4Corr_, rVsCaloImbProjC4_8Corr_, rVsCaloImbProjC8_100Corr_);
	  else
	    getPtProj(trkPt_[trkEntry], trkPtCorrVsCalo_[trkEntry], trkPhi_[trkEntry], VsCaloLeadJtPhi_, rVsCaloImbProjNCFCorr_, rVsCaloImbPerpNCFCorr_, rVsCaloImbProjNC0_1Corr_, rVsCaloImbProjNC1_2Corr_, rVsCaloImbProjNC2_4Corr_, rVsCaloImbProjNC4_8Corr_, rVsCaloImbProjNC8_100Corr_);
	}
      }
    }

    if(montecarlo && truthSet_){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	getPtProj(trkPt_[trkEntry], trkPtCorrT_[trkEntry], trkPhi_[trkEntry], TLeadJtPhi_, rTImbProjFCorr_, rTImbPerpFCorr_, rTImbProj0_1Corr_, rTImbProj1_2Corr_, rTImbProj2_4Corr_, rTImbProj4_8Corr_, rTImbProj8_100Corr_);  

	if(trkRLeadT_[trkEntry] > 0 && trkRSubLeadT_[trkEntry] > 0){
	  if(trkRLeadT_[trkEntry] < .8 || trkRSubLeadT_[trkEntry] < .8)
	    getPtProj(trkPt_[trkEntry], trkPtCorrT_[trkEntry], trkPhi_[trkEntry], TLeadJtPhi_, rTImbProjCFCorr_, rTImbPerpCFCorr_, rTImbProjC0_1Corr_, rTImbProjC1_2Corr_, rTImbProjC2_4Corr_, rTImbProjC4_8Corr_, rTImbProjC8_100Corr_);
	  else
	    getPtProj(trkPt_[trkEntry], trkPtCorrT_[trkEntry], trkPhi_[trkEntry], TLeadJtPhi_, rTImbProjNCFCorr_, rTImbPerpNCFCorr_, rTImbProjNC0_1Corr_, rTImbProjNC1_2Corr_, rTImbProjNC2_4Corr_, rTImbProjNC4_8Corr_, rTImbProjNC8_100Corr_);
	}
      }
    }

    if(montecarlo){
      //Iterate over truth

      nGen_ = 0;

      GenParticles genCollection;
      genCollection = c->genparticle;

      for(Int_t genEntry = 0; genEntry < genCollection.mult; genEntry++){
	if(truthSet_)	  TTotGen++;
	if(recoPFSet_)	  PFTotGen++;
	if(recoCaloSet_)	  CaloTotGen++;
	  
	if(genCollection.chg[genEntry] == 0){
	  if(truthSet_)	    TGenChgCut++;
	  if(recoPFSet_)	    PFGenChgCut++;
	  if(recoCaloSet_)	    CaloGenChgCut++;
	  
	  continue;
	}
	  
	if(TMath::Abs(genCollection.eta[genEntry]) > 2.4){
	  if(truthSet_)	    TGenEtaCut++;
	  if(recoPFSet_)	    PFGenEtaCut++;
	  if(recoCaloSet_)	    CaloGenEtaCut++;

	  continue;
	}
	
	if(genCollection.pt[genEntry] < 0.5){
	  if(truthSet_)	    TGenPtCut++;
	  if(recoPFSet_)	    PFGenPtCut++;
	  if(recoCaloSet_)	    CaloGenPtCut++;

	  continue;
	}

	genPt_[nGen_] = genCollection.pt[genEntry];
	genPhi_[nGen_] = genCollection.phi[genEntry];
	genEta_[nGen_] = genCollection.eta[genEntry];

	if(recoPFSet_)
          genPtPF_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], PFLeadJtPhi_));
	else
	  genPtPF_[nGen_] = 0;

	if(recoCaloSet_)
          genPtCalo_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], CaloLeadJtPhi_));
	else
	  genPtCalo_[nGen_] = 0;

	if(truthSet_)
          genPtT_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], TLeadJtPhi_));
	else
	  genPtT_[nGen_] = 0;

	if(recoVsPFSet_)
          genPtVsPF_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], VsPFLeadJtPhi_));
	else
	  genPtVsPF_[nGen_] = 0;

	if(recoVsCaloSet_)
          genPtVsCalo_[nGen_] = -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], VsCaloLeadJtPhi_));
	else
	  genPtVsCalo_[nGen_] = 0;

	
	if(truthSet_)
	  genLeadDelPhi_[nGen_] = getAbsDphi(TLeadJtPhi_, genCollection.phi[genEntry]);
	else
	  genLeadDelPhi_[nGen_] = -10;


	if(truthSet_)
	  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], TLeadJtPhi_, gTImbProjF_, gTImbPerpF_, gTImbProj0_1_, gTImbProj1_2_, gTImbProj2_4_, gTImbProj4_8_, gTImbProj8_100_);
	
	if(recoPFSet_)
	  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], PFLeadJtPhi_, gPFImbProjF_, gPFImbPerpF_, gPFImbProj0_1_, gPFImbProj1_2_, gPFImbProj2_4_, gPFImbProj4_8_, gPFImbProj8_100_);

	if(recoCaloSet_)
	  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], CaloLeadJtPhi_, gCaloImbProjF_, gCaloImbPerpF_, gCaloImbProj0_1_, gCaloImbProj1_2_, gCaloImbProj2_4_, gCaloImbProj4_8_, gCaloImbProj8_100_);

	if(recoVsPFSet_)
	  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], VsPFLeadJtPhi_, gVsPFImbProjF_, gVsPFImbPerpF_, gVsPFImbProj0_1_, gVsPFImbProj1_2_, gVsPFImbProj2_4_, gVsPFImbProj4_8_, gVsPFImbProj8_100_);

	if(recoVsCaloSet_)
	  getPtProj(genCollection.pt[genEntry], genCollection.pt[genEntry], genCollection.phi[genEntry], VsCaloLeadJtPhi_, gVsCaloImbProjF_, gVsCaloImbPerpF_, gVsCaloImbProj0_1_, gVsCaloImbProj1_2_, gVsCaloImbProj2_4_, gVsCaloImbProj4_8_, gVsCaloImbProj8_100_);
	  

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
  if(argc != 4)
    {
      std::cout << "Usage: sortForest <inputFile> <MCBool> <outputFile>" << std::endl;
      return 1;
    }

  int rStatus = -1;

  rStatus = makeDiJetTTree(argv[1], sampleType(atoi(argv[2])), argv[3]);

  return rStatus;
}
