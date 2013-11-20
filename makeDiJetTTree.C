//=============================================
// Author: Chris McGinn
// 
// DiJet Analysis Class (MC)
//
//=============================================

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "../HiForestAnalysis/hiForest.h"
#include "../gammaJetAnalysis/commonUtility.h"
#include "cfmDiJetSkim.h"

const Float_t leadJtPtCut = 120.;
const Float_t subLeadJtPtCut = 50.;
const Float_t jtDelPhiCut = 7.*(TMath::Pi())/8.;
const Float_t jtEtaCut = 10.;

collisionType getCType(sampleType sType);

int makeDiJetTTree(const char* inName, sampleType sType, const char *outName)
{
  bool montecarlo = false;
  if(sType == kPPMC || sType == kPAMC || sType == kHIMC)
    montecarlo = true;

  collisionType cType = getCType(sType);

  TFile *outFile = new TFile(outName, "RECREATE");

  InitDiJetSkim(1);

  HiForest *c = new HiForest(inName, "Forest", cType, montecarlo);

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

  Int_t totTrk = 0;
  Int_t purityCut = 0;
  Int_t trkEtaCut = 0;
  Int_t trkPtCut = 0;

  Int_t totGen = 0;
  Int_t genEtaCut = 0;
  Int_t genPtCut = 0;
  Int_t genChgCut = 0;

  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    c->GetEntry(jentry);

    totEv++;

    if(jentry%10000 == 0)
      std::cout << jentry << std::endl;

    if(!c->selectEvent()){
      selectCut++;
      continue;
    }


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

    if(leadJtIndex < 0)
      rLeadJtPtCut++;
    else if(subLeadJtIndex < 0)
      rSubLeadJtPtCut++;
    else if(getAbsDphi(c->akPu3PF.jtphi[leadJtIndex], c->akPu3PF.jtphi[subLeadJtIndex]) < jtDelPhiCut)
      rDelPhiCut++;
    else if(TMath::Abs(c->akPu3PF.jteta[leadJtIndex]) > jtEtaCut || TMath::Abs(c->akPu3PF.jteta[subLeadJtIndex]) > jtEtaCut)
      rJtEtaCut++;


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
      continue;
    }
    else if(subLeadJtIndex < 0){
      gSubLeadJtPtCut++;
      continue;
    }

    gLeadJtPhi_ = c->akPu3PF.refphi[leadJtIndex];
    gSubLeadJtPhi_ = c->akPu3PF.refphi[subLeadJtIndex];

    if(getAbsDphi(gLeadJtPhi_, gSubLeadJtPhi_) < jtDelPhiCut){
      gDelPhiCut++;
      continue;
    }

    gLeadJtEta_ = c->akPu3PF.refeta[leadJtIndex];
    gSubLeadJtEta_ = c->akPu3PF.refeta[subLeadJtIndex];

    if(TMath::Abs(gLeadJtEta_) > jtEtaCut || TMath::Abs(gSubLeadJtEta_) > jtEtaCut){
      gJtEtaCut++;
      continue;
    }

    rLeadJtPt_ = c->akPu3PF.jtpt[leadJtIndex];
    rSubLeadJtPt_ = c->akPu3PF.jtpt[subLeadJtIndex];
    rLeadJtPhi_ = c->akPu3PF.jtphi[leadJtIndex];
    rSubLeadJtPhi_ = c->akPu3PF.jtphi[subLeadJtIndex];
    rLeadJtEta_ = c->akPu3PF.jteta[leadJtIndex];
    rSubLeadJtEta_ = c->akPu3PF.jteta[subLeadJtIndex];

    run_ = c->evt.run;
    evt_ = c->akPu3PF.evt;
    lumi_ = c->evt.lumi;
    hiBin_ = c->evt.hiBin;

    //Iterate over tracks

    nTrk_ = 0;
    rImbProjF_ = 0;
    rImbProjH_ = 0;
    rImbProjL_ = 0;

    rImbPerpF_ = 0;
    rImbPerpH_ = 0;
    rImbPerpL_ = 0;

    if(montecarlo){
      for(Int_t divIter = 0; divIter < 10; divIter++){
	rDivGPt_[divIter] = 0;
      }
    }

    Tracks trkCollection;
    trkCollection = c->track;

    for(Int_t trkEntry = 0; trkEntry < trkCollection.nTrk; trkEntry++){
      totTrk++;

      if(TMath::Abs(trkCollection.trkEta[trkEntry]) > 2.4){
	trkEtaCut++;
        continue;
      }

      if(trkCollection.trkPt[trkEntry] < 0.5){
	trkPtCut++;
        continue;
      }

      if(!trkCollection.highPurity[trkEntry]){ //Note highPuritySetWithPV seems to be wrong cut, creates diff. bet. truth and reco
	purityCut++;
        continue;
      }

      trkPt_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPhi_[nTrk_] = trkCollection.trkPhi[trkEntry];
      trkEta_[nTrk_] = trkCollection.trkEta[trkEntry];
      trkLeadDelPhi_[nTrk_] = getAbsDphi(rLeadJtPhi_, trkCollection.trkPhi[trkEntry]);

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

      rImbProjF_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
      rImbPerpF_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
      if(trkCollection.trkPt[trkEntry] > 8){
	rImbProjH_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	rImbPerpH_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
      }
      else{
	rImbProjL_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
	rImbPerpL_ += -trkCollection.trkPt[trkEntry]*sin(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
      }
      nTrk_++;
      if(nTrk_ > MAXTRKS - 1){
        printf("ERROR: Trk arrays not large enough.\n");
        return(1);
      }
    }

    Float_t rDivGPt_temp[10];

    if(montecarlo){
      //Iterate over truth

      nGen_ = 0;
      gImbProjF_ = 0;
      gImbProjH_ = 0;
      gImbProjL_ = 0;

      gImbPerpF_ = 0;
      gImbPerpH_ = 0;
      gImbPerpL_ = 0;

      for(Int_t divIter = 0; divIter < 10; divIter++){
	rDivGPt_temp[divIter] = 0;
      }

      GenParticles genCollection;
      genCollection = c->genparticle;

      for(Int_t genEntry = 0; genEntry < genCollection.mult; genEntry++){
	totGen++;

	if(genCollection.chg[genEntry] == 0){
	  genChgCut++;
	  continue;
	}

	if(TMath::Abs(genCollection.eta[genEntry]) > 2.4){
	  genEtaCut++;
	  continue;
	}
	
	if(genCollection.pt[genEntry] < 0.5){
	  genPtCut++;
	  continue;
	}
		
	genPt_[nGen_] = genCollection.pt[genEntry];
	genPhi_[nGen_] = genCollection.phi[genEntry];
	genEta_[nGen_] = genCollection.eta[genEntry];
	genLeadDelPhi_[nGen_] = getAbsDphi(gLeadJtPhi_, genCollection.phi[genEntry]);

	for(Int_t divIter = 0; divIter < 10; divIter++){
	  if(genPt_[nGen_] > 20.)
	    break;
	  else if(2*divIter < genPt_[nGen_] && 2*(divIter + 1.) > genPt_[nGen_]){
	    rDivGPt_temp[divIter]++;
	    break;
	  }
	}
	
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
  std::cout << "totTrk: " << totTrk << std::endl;
  tempTot = totTrk - trkEtaCut;
  std::cout << "trkEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - trkPtCut;
  std::cout << "trkPtCut: " << tempTot << std::endl;
  tempTot = tempTot - purityCut;
  std::cout << "purityCut: " << tempTot << std::endl;

  std::cout << std::endl;
  std::cout << "totGen: " << totGen << std::endl;
  tempTot = totGen - genChgCut;
  std::cout << "genChgCut: " << tempTot << std::endl;
  tempTot = tempTot - genEtaCut;
  std::cout << "genEtaCut: " << tempTot << std::endl;
  tempTot = tempTot - genPtCut;
  std::cout << "genPtCut: " << tempTot << std::endl;
  
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
