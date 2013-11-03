#include "cfmDiJetSkim.h"
#include "../gammaJetAnalysis/commonUtility.h"
#include "TDatime.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TCut.h"

TFile* inFile_p = 0;
TFile* outFile_p = 0;

TTree* inTree_p = 0;


void makeAsymmHist(TTree* getTree_p, const char* outName, Int_t binLow, Int_t binHi)
{
  TH1F* asymmHist_p;

  TCut centCut = Form("hiBin >= %d && hiBin <= %d", binLow, binHi);
  getTree_p->Draw("(leadJtPt - subLeadJtPt)/(leadJtPt + subLeadJtPt) >> asymmHist", centCut);
  asymmHist_p = (TH1F*)inFile_p->Get("asymmHist");

  outFile_p = new TFile(outName, "UPDATE");
  asymmHist_p->Write();
  outFile_p->Close();
  delete outFile_p;
}


void cfmDiJetHist(const char* inName, const char* outName)
{
  inFile_p = new TFile(inName, "READ");
  inTree_p = (TTree*)inFile_p->Get("jetTree");

  makeAsymmHist(inTree_p, outName, 0, 10);

  inFile_p->Close();
  delete inFile_p;
}


