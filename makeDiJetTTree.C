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

  Long64_t nentries = c->GetEntries();

  for(Long64_t jentry = 0; jentry < 1000 /*nentries*/; jentry++){
    c->GetEntry(jentry);

    if(jentry%10000 == 0)
      std::cout << jentry << std::endl;

    if(!c->selectEvent())
      continue;

    if(c->akPu3PF.nref < 2)
      continue;

    Int_t leadJtIndex = -1;
    Int_t subLeadJtIndex = -1;
    leadJtPt_ = subLeadJtPtCut;
    subLeadJtPt_ = subLeadJtPtCut;
    for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.nref; jtEntry++){
      if(c->akPu3PF.refpt[jtEntry] > leadJtPtCut && c->akPu3PF.refpt[jtEntry] > leadJtPt_){
	subLeadJtIndex = leadJtIndex;
	subLeadJtPt_ = leadJtPt_;
	leadJtIndex = jtEntry;
	leadJtPt_ = c->akPu3PF.refpt[jtEntry];
      }
      else if(c->akPu3PF.refpt[jtEntry] > subLeadJtPt_){
	subLeadJtIndex = jtEntry;
	subLeadJtPt_ = c->akPu3PF.refpt[jtEntry];
      }
    }

    if(leadJtIndex < 0 || subLeadJtIndex < 0)
      continue;

    leadJtPhi_ = c->akPu3PF.refphi[leadJtIndex];
    subLeadJtPhi_ = c->akPu3PF.refphi[subLeadJtIndex];

    if(getAbsDphi(leadJtPhi_, subLeadJtPhi_) < 2*(TMath::Pi())/3)
      continue;

    leadJtEta_ = c->akPu3PF.refeta[leadJtIndex];
    subLeadJtEta_ = c->akPu3PF.refeta[subLeadJtIndex];

    run_ = c->evt.run;
    evt_ = c->akPu3PF.evt;
    lumi_ = c->evt.lumi;
    hiBin_ = c->evt.hiBin;

    //    jetTree_p->Fill();
  }

  
  outFile->cd();
  //  jetTree_p->Write();
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
