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

const Float_t leadJtPtCut = 120.;
const Float_t subLeadJtPtCut = 50.;
const Float_t jtDelPhiCut = 7.*(TMath::Pi())/8.;
const Float_t jtEtaCut = 2.4; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed

collisionType getCType(sampleType sType);

int makeDiJetTTree(string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_CFMSKIM.root")
{
  //Define MC or Data
  bool montecarlo = false;
  if(sType == kPPMC || sType == kPAMC || sType == kHIMC)
    montecarlo = true;

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
  c->hasGenParticleTree = true;

  Long64_t nentries = c->GetEntries();

  Int_t totEv = 0;
  Int_t selectCut = 0;

  Int_t gLeadJtPtCut = 0;
  Int_t gSubLeadJtPtCut = 0;
  Int_t gDelPhiCut = 0;
  Int_t gJtEtaCut = 0;

  Int_t rLeadJtPtCut = 0;
  Int_t rSubLeadJtPtCut = 0;
  Int_t rDelPhiCut = 0;
  Int_t rJtEtaCut = 0;

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

  Int_t rTotTrk = 0;
  Int_t rTrkEtaCut = 0;
  Int_t rTrkPtCut = 0;
  Int_t rPurityCut = 0;
  Int_t rTrkDzCut = 0;
  Int_t rTrkDxyCut = 0;
  Int_t rTrkPtErrorCut = 0;

  Int_t rTotGen = 0;
  Int_t rGenEtaCut = 0;
  Int_t rGenPtCut = 0;
  Int_t rGenChgCut = 0;

  defTrkCorr();

  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    c->GetEntry(jentry);

    totEv++;

    Bool_t gEventPass = false;
    Bool_t rEventPass = false;

    if(jentry%10000 == 0)
      std::cout << jentry << std::endl;

    if(!c->selectEvent()){
      selectCut++;
      continue;
    }

    InitJetVar(montecarlo);

    Int_t leadJtIndex = -1;
    Int_t subLeadJtIndex = -1;
    rLeadJtPt_ = subLeadJtPtCut;
    rSubLeadJtPt_ = subLeadJtPtCut;
    for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.nref; jtEntry++){
      if(c->akPu3PF.jtpt[jtEntry] > leadJtPtCut && c->akPu3PF.jtpt[jtEntry] > rLeadJtPt_){
	subLeadJtIndex = leadJtIndex;
	rSubLeadJtPt_ = rLeadJtPt_;
	leadJtIndex = jtEntry;
	rLeadJtPt_ = c->akPu3PF.jtpt[jtEntry];
      }
      else if(c->akPu3PF.jtpt[jtEntry] > rSubLeadJtPt_){
	subLeadJtIndex = jtEntry;
	rSubLeadJtPt_ = c->akPu3PF.jtpt[jtEntry];
      }
    }

    if(leadJtIndex < 0){
      rLeadJtPtCut++;
      rLeadJtPt_ = -10;
      rSubLeadJtPt_ = -10;
    }
    else if(subLeadJtIndex < 0){
      rSubLeadJtPtCut++;
      rSubLeadJtPt_ = -10;
    }
    else if(getAbsDphi(c->akPu3PF.jtphi[leadJtIndex], c->akPu3PF.jtphi[subLeadJtIndex]) < jtDelPhiCut){
      rDelPhiCut++;
    }
    else if(TMath::Abs(c->akPu3PF.jteta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akPu3PF.jteta[subLeadJtIndex]) > jtEtaCut){
      rJtEtaCut++;
    }
    else{
      rLeadJtPhi_ = c->akPu3PF.jtphi[leadJtIndex];
      rSubLeadJtPhi_ = c->akPu3PF.jtphi[subLeadJtIndex];
      rLeadJtEta_ = c->akPu3PF.jteta[leadJtIndex];
      rSubLeadJtEta_ = c->akPu3PF.jteta[subLeadJtIndex];

      rEventPass = true;
    }

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

      gRLeadJtPt_ = c->akPu3PF.jtpt[leadJtIndex];
      gRSubLeadJtPt_ = c->akPu3PF.jtpt[subLeadJtIndex];
      gRLeadJtPhi_ = c->akPu3PF.jtphi[leadJtIndex];
      gRSubLeadJtPhi_ = c->akPu3PF.jtphi[subLeadJtIndex];
      gRLeadJtEta_ = c->akPu3PF.jteta[leadJtIndex];
      gRSubLeadJtEta_ = c->akPu3PF.jteta[subLeadJtIndex];

      gEventPass = true;
    }

    if(gEventPass == false && rEventPass == false)
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
      if(rEventPass)	rTotTrk++;
      
      if(TMath::Abs(trkCollection.trkEta[trkEntry]) > 2.4){
	if(gEventPass)	  gTrkEtaCut++;
	if(rEventPass)	  rTrkEtaCut++;

	continue;
      }

      if(trkCollection.trkPt[trkEntry] < 0.9){
	if(gEventPass)	  gTrkPtCut++;
	if(rEventPass)	  rTrkPtCut++;

	continue;
      }
      
      if(!trkCollection.highPurity[trkEntry]){ //Note highPuritySetWithPV seems to be wrong cut, creates diff. bet. truth and reco
	if(gEventPass)	  gPurityCut++;
	if(rEventPass)	  rPurityCut++;

	continue;
      }

      if(TMath::Abs(trkCollection.trkDz1[trkEntry]/trkCollection.trkDzError1[trkEntry]) > 3){
	if(gEventPass)	  gTrkDzCut++;
	if(rEventPass)	  rTrkDzCut++;

	continue;
      }

      if(TMath::Abs(trkCollection.trkDxy1[trkEntry]/trkCollection.trkDxyError1[trkEntry]) > 3){
	if(gEventPass)	  gTrkDxyCut++;
	if(rEventPass)	  rTrkDxyCut++;

	continue;
      }

      if(trkCollection.trkPtError[trkEntry]/trkCollection.trkPt[trkEntry] > 0.1){
	if(gEventPass)	  gTrkPtErrorCut++;
	if(rEventPass)	  rTrkPtErrorCut++;

	continue;
      }

      trkPt_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPhi_[nTrk_] = trkCollection.trkPhi[trkEntry];
      trkEta_[nTrk_] = trkCollection.trkEta[trkEntry];
      
      if(rLeadJtPhi_ > -9)
	trkLeadDelPhi_[nTrk_] = getAbsDphi(rLeadJtPhi_, trkCollection.trkPhi[trkEntry]);
      else
	trkLeadDelPhi_[nTrk_] = -10;

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
      
    
      if(rLeadJtPhi_ > -9){
	rImbProjF_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rLeadJtPhi_));
	rImbPerpF_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], rLeadJtPhi_));
	if(trkCollection.trkPt[trkEntry] > 8){
	  rImbProjH_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rLeadJtPhi_));
	  rImbPerpH_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], rLeadJtPhi_));
	}
	else{
	  rImbProjL_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], rLeadJtPhi_));
	  rImbPerpL_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], rLeadJtPhi_));
	}
      }

      if(montecarlo && gLeadJtPhi_ > -9){
	gRImbProjF_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	gRImbPerpF_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	if(trkCollection.trkPt[trkEntry] > 8){
	  gRImbProjH_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	  gRImbPerpH_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	}
	else{
	  gRImbProjL_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	  gRImbPerpL_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	}
      }

      nTrk_++;
      if(nTrk_ > MAXTRKS - 1){
	printf("ERROR: Trk arrays not large enough.\n");
	return(1);
      }
    }

    for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
      if(trkPt_[trkEntry] < 4.9){
	for(Int_t ptIter = 0; ptIter < 40; ptIter++){
  
	  if(trkPt_[trkEntry] >= ((Float_t)ptIter)/10 + 0.9 && trkPt_[trkEntry] < ((Float_t)ptIter)/10 + 1.){
	    trkPtCorr_[trkEntry] = trkPt_[trkEntry]/trkCorrLow_[ptIter];
	    trkPtFact_[trkEntry] = trkCorrLow_[ptIter];
	    break;
	  }

	}
      }
      else{
	trkPtCorr_[trkEntry] = trkPt_[trkEntry];
	trkPtFact_[trkEntry] = 1;
      }
    
      if(rLeadJtPhi_ > -9){
	rImbProjFCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], rLeadJtPhi_));
	rImbPerpFCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], rLeadJtPhi_));
	if(trkPtCorr_[trkEntry] > 8){
	  rImbProjHCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], rLeadJtPhi_));
	  rImbPerpHCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], rLeadJtPhi_));
	}
	else{
	  rImbProjLCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], rLeadJtPhi_));
	  rImbPerpLCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], rLeadJtPhi_));
	}
      }

      if(montecarlo && gLeadJtPhi_ > -9){
	gRImbProjFCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));
	gRImbPerpFCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));
	if(trkPtCorr_[trkEntry] > 8){
	  gRImbProjHCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));
	  gRImbPerpHCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));
	}
	else{
	  gRImbProjLCorr_ += -trkPtCorr_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));
	  gRImbPerpLCorr_ += -trkPtCorr_[trkEntry]*sin(getDPHI(trkPhi_[trkEntry], gLeadJtPhi_));
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
	if(rEventPass)	  rTotGen++;
	  
	if(genCollection.chg[genEntry] == 0){
	  if(gEventPass)	    gGenChgCut++;
	  if(rEventPass)	    rGenChgCut++;
	  
	  continue;
	}
	  
	if(TMath::Abs(genCollection.eta[genEntry]) > 2.4){
	  if(gEventPass)	    gGenEtaCut++;
	  if(rEventPass)	    rGenEtaCut++;

	  continue;
	}
	
	if(genCollection.pt[genEntry] < 0.5){
	  if(gEventPass)	    gGenPtCut++;
	  if(rEventPass)	    rGenPtCut++;

	  continue;
	}
	
	genPt_[nGen_] = genCollection.pt[genEntry];
	genPhi_[nGen_] = genCollection.phi[genEntry];
	genEta_[nGen_] = genCollection.eta[genEntry];
	
	if(gLeadJtPhi_ > -9)
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

	if(gLeadJtPhi_ > -9){
	  gImbProjF_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	  gImbPerpF_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	  if(genCollection.pt[genEntry] > 8){
	    gImbProjH_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	    gImbPerpH_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	  }
	  else{
	    gImbProjL_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	    gImbPerpL_ += -genCollection.pt[genEntry]*sin(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
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
  tempTot = totEv - selectCut - rLeadJtPtCut;
  std::cout << "rLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - rSubLeadJtPtCut;
  std::cout << "rSubLeadJtPtCut: " << tempTot << std::endl;
  tempTot = tempTot - rDelPhiCut;
  std::cout << "rDelPhiCut: " << tempTot << std::endl;
  tempTot = tempTot - rJtEtaCut;
  std::cout << "rJtEtaCut: " << tempTot << std::endl;

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
  std::cout << "rTotTrk: " << rTotTrk << std::endl;
  tempTot = rTotTrk - rTrkEtaCut;
  std::cout << "rTrkEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - rTrkPtCut;
  std::cout << "rTrkPtCut: " << tempTot << std::endl;
  tempTot = tempTot - rPurityCut;
  std::cout << "rPurityCut: " << tempTot << std::endl;
  tempTot = tempTot - rTrkDzCut;
  std::cout << "rTrkDzCut: " << tempTot << std::endl;
  tempTot = tempTot - rTrkDxyCut;
  std::cout << "rTrkDxyCut: " << tempTot << std::endl;
  tempTot = tempTot - rTrkPtErrorCut;
  std::cout << "rTrkPtErrorCut: " << tempTot << std::endl;

  std::cout << std::endl;
  std::cout << "rTotGen: " << rTotGen << std::endl;
  tempTot = rTotGen - rGenChgCut;
  std::cout << "rGenChgCut: " << tempTot << std::endl;
  tempTot = tempTot - rGenEtaCut;
  std::cout << "rGenEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - rGenPtCut;
  std::cout << "rGenPtCut: " << tempTot << std::endl;



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
