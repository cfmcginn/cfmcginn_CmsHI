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

  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    c->GetEntry(jentry);

    if(jentry%10000 == 0)
      std::cout << jentry << std::endl;

    if(!c->selectEvent())
      continue;

    if(c->akPu3PF.nref < 2)
      continue;

    Int_t leadJtIndex = -1;
    Int_t subLeadJtIndex = -1;
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

    if(leadJtIndex < 0 || subLeadJtIndex < 0)
      continue;

    gLeadJtPhi_ = c->akPu3PF.refphi[leadJtIndex];
    gSubLeadJtPhi_ = c->akPu3PF.refphi[subLeadJtIndex];

    if(getAbsDphi(gLeadJtPhi_, gSubLeadJtPhi_) < 2*(TMath::Pi())/3)
      continue;

    gLeadJtEta_ = c->akPu3PF.refeta[leadJtIndex];
    gSubLeadJtEta_ = c->akPu3PF.refeta[subLeadJtIndex];

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

    Tracks trkCollection;
    trkCollection = c->track;

    for(Int_t trkEntry = 0; trkEntry < trkCollection.nTrk; trkEntry++){
      if(!trkCollection.highPuritySetWithPV[trkEntry])
        continue;

      if(TMath::Abs(trkCollection.trkEta[trkEntry]) > 2.4)
        continue;

      if(trkCollection.trkPt[trkEntry] < 0.5)
        continue;

      trkPt_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPhi_[nTrk_] = trkCollection.trkPhi[trkEntry];
      trkEta_[nTrk_] = trkCollection.trkEta[trkEntry];

      rImbProjF_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
      if(trkCollection.trkPt[trkEntry] > 8)
	rImbProjH_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));
      else
	rImbProjL_ += -trkCollection.trkPt[trkEntry]*cos(getDPHI(trkCollection.trkPhi[trkEntry], gLeadJtPhi_));

      nTrk_++;
      if(nTrk_ > MAXTRKS - 1){
        printf("ERROR: Trk arrays not large enough.\n");
        return(1);
      }
    }


    if(montecarlo){
      //Iterate over truth

      nGen_ = 0;
      gImbProjF_ = 0;
      gImbProjH_ = 0;
      gImbProjL_ = 0;

      GenParticles genCollection;
      genCollection = c->genparticle;

      for(Int_t genEntry = 0; genEntry < genCollection.mult; genEntry++){
	if(TMath::Abs(genCollection.eta[genEntry]) > 2.4)
	  continue;
	
	if(genCollection.pt[genEntry] < 0.5)
	  continue;
	
	if(genCollection.chg[genEntry] == 0)
	  continue;
	
	genPt_[nGen_] = genCollection.pt[genEntry];
	genPhi_[nGen_] = genCollection.phi[genEntry];
	genEta_[nGen_] = genCollection.eta[genEntry];
	
	gImbProjF_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	if(genCollection.pt[genEntry] > 8)
	  gImbProjH_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));
	else
	  gImbProjL_ += -genCollection.pt[genEntry]*cos(getDPHI(genCollection.phi[genEntry], gLeadJtPhi_));

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
