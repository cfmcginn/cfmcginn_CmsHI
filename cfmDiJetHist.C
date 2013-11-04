#include "cfmDiJetSkim.h"
#include "../gammaJetAnalysis/commonUtility.h"
#include "TDatime.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TCut.h"

TFile* inFile_p = 0;
TFile* outFile_p = 0;

TTree* inTree_p = 0;


void niceTH1(TH1F* uglyTH1, float max , float min, float ndivX, float ndivY)
{
  handsomeTH1N(uglyTH1);
  uglyTH1->SetMaximum(max);
  uglyTH1->SetMinimum(min);
  uglyTH1->SetNdivisions(ndivX);
  uglyTH1->SetNdivisions(ndivY, "Y");
}


void makeAsymmHist(TTree* getTree_p, const char* outName, const char* gorr, Int_t nBins, Int_t histLow, Int_t histHi, Int_t centLow, Int_t centHi)
{
  inFile_p->cd();

  const char* title = Form("%sAsymm_%d%d", gorr, (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5));

  TCanvas* asymmCanvas_p = new TCanvas(Form("%s_c", title), Form("%s_c", title), 1);
  TH1F* asymmHist_p;

  TString name = Form("(%sLeadJtPt - %sSubLeadJtPt)/(%sLeadJtPt + %sSubLeadJtPt) >> %s_h(%d, %d, %d)", gorr, gorr, gorr, gorr, title, nBins, histLow, histHi);
  TCut centCut = Form("hiBin >= %d && hiBin <= %d && %sLeadJtPt > 120 && %sSubLeadJtPt > 50", centLow, centHi, gorr, gorr);
  getTree_p->Draw(name, centCut);
  asymmHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  TLatex* cent_p;

  if(*gorr == 103){
    niceTH1(asymmHist_p, .5, 0., 405, 505);
    cent_p = new TLatex(.6, .3, Form("Truth, %d-%d", (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5)));
  }
  else if(*gorr == 114){
    niceTH1(asymmHist_p, .4, 0., 405, 504);
    cent_p = new TLatex(.6, .3, Form("Reco, %d-%d", (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5)));
  }
  else{
    std::cout << "Error: Neither 'g' or 'r' input." << std::endl;
    return;
  }

  asymmHist_p->SetYTitle("Event Fraction");
  asymmHist_p->SetXTitle("A_{J} = (p_{T,1} - p_{T,2})/(p_{T,1} + p_{T,2})");
  cent_p->SetNDC();
  cent_p->Draw();

  outFile_p = new TFile(outName, "UPDATE");
  asymmHist_p->Write();
  asymmCanvas_p->Write();
  outFile_p->Close();

  delete outFile_p;
  delete cent_p;
  delete asymmCanvas_p;
}


void makeAsymmPanel(const char* fileName, const char* gorr)
{

}


void cfmDiJetHist(const char* inName, const char* outName)
{
  inFile_p = new TFile(inName, "READ");
  inTree_p = (TTree*)inFile_p->Get("jetTree");

  makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 0, 39);
  makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 0, 3);
  makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 4, 7);
  makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 8, 11);
  makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 12, 19);
  makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 20, 39);

  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 0, 39);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 0, 3);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 4, 7);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 8, 11);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 12, 19);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 20, 39);
  
  inFile_p->Close();
  delete inFile_p;
}


