#include <iostream>
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCut.h"
#include "commonUtility.h"

TFile* inFilePtHat_p[5];

void pthatMergePlots(const char* outName)
{
  TH1::SetDefaultSumw2();

  TFile* outFile_p = new TFile(outName, "UPDATE");

  const char* fileName[5];
  fileName[0] = "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet30_FOREST_Track8_Jet24_FixedPtHatJES_v0/0.root";
  fileName[1] = "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet50_FOREST_Track8_Jet24_FixedPtHatJES_v0/0.root";
  fileName[2] = "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHatJES_v0/0.root";
  fileName[3] = "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHatJES_v0/0.root";
  fileName[4] = "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHatJES_v0/0.root";

  TChain* ptHatChain_p = new TChain("akVs3CaloJetAnalyzer/t");
  TChain* vzChain_p = new TChain("hiEvtAnalyzer/HiTree");

  for(unsigned int iter = 0; iter < 5; iter++){
    ptHatChain_p->Add(fileName[iter]);
    vzChain_p->Add(fileName[iter]);
  }

  ptHatChain_p->AddFriend(vzChain_p);

  ptHatChain_p->Draw("pthat >> ptHatFull_h(100, 0, 500)");

  TH1F* getPtHatFull_p = (TH1F*)outFile_p->Get("ptHatFull_h");

  TCanvas* ptHatFull_p = new TCanvas("ptHatFull_p", "ptHatFull_p", 1);
  gPad->SetLogy();
  handsomeTH1(getPtHatFull_p);
  getPtHatFull_p->DrawCopy("E1 HIST");

  ptHatFull_p->Write();
  claverCanvasSaving(ptHatFull_p, "pngDir/ptHatFull_c", "png");

  ptHatChain_p->Draw("vz >> vzNoWeight");
  TH1F* getVzNoWeight_p = (TH1F*)outFile_p->Get("vzNoWeight");

  std::cout << "Vz Mean, Error, no weight: " << getVzNoWeight_p->GetMean() << ", " << getVzNoWeight_p->GetMeanError() << std::endl;

  Float_t crossSections[6] = {.01075, .001025, .00009865, .00003069, .00001129, 0.000000000};
  Int_t ptHatCuts[6] = {30, 50, 80, 100, 120, 100000};
  Float_t numEntries[5];
  Float_t weight[5];

  TH1F* getHatHist_p[5];
  TH1F* getVzHist_p[5];

  TH1F* getPtHatFullWeight_p = new TH1F("getPtHatFullWeight_p", "getPtHatFullWeight_p", 100, 0, 500);
  TH1F* getVzFullWeight_p = new TH1F("getVzFullWeight_p", "getVzFullWeight_p", 100, -100, 100);

  for(unsigned int hatIter = 0; hatIter < 5; hatIter++){
    TCut ptHatCut = Form("pthat > %d && pthat < %d", ptHatCuts[hatIter], ptHatCuts[hatIter+1]);
    std::cout << ptHatCut << std::endl;
    ptHatChain_p->Draw(Form("pthat >> %dhat_h(100, 0, 500)", hatIter), ptHatCut);

    getHatHist_p[hatIter] = (TH1F*)outFile_p->Get(Form("%dhat_h", hatIter));

    numEntries[hatIter] = getHatHist_p[hatIter]->GetEntries();

    weight[hatIter] = (crossSections[hatIter] - crossSections[hatIter + 1])/numEntries[hatIter];

    std::cout << weight[hatIter] << std::endl;

    getPtHatFullWeight_p->Add(getHatHist_p[hatIter], weight[hatIter]);

    ptHatChain_p->Draw(Form("vz >> %dvz_h(100, -100, 100)", hatIter), ptHatCut);
    getVzHist_p[hatIter] = (TH1F*)outFile_p->Get(Form("%dvz_h", hatIter));
    getVzFullWeight_p->Add(getVzHist_p[hatIter], weight[hatIter]);

  }

  TCanvas* ptHatFullWeight_p = new TCanvas("ptHatFullWeight_p", "ptHatFullWeight_p", 1);
  gPad->SetLogy();
  handsomeTH1(getPtHatFullWeight_p);
  getPtHatFullWeight_p->DrawCopy("E1 HIST");

  ptHatFullWeight_p->Write();
  claverCanvasSaving(ptHatFullWeight_p, "pngDir/ptHatFullWeight_c", "png");

  std::cout << "Vz Mean, Error, weighted: " << getVzFullWeight_p->GetMean() << ", " << getVzFullWeight_p->GetMeanError() << std::endl;

  getVzFullWeight_p->Write();

  delete ptHatFullWeight_p;
  delete getVzFullWeight_p;
  delete getPtHatFullWeight_p;
  delete ptHatFull_p;
  delete vzChain_p;
  delete ptHatChain_p;

  outFile_p->Close();

  delete outFile_p;

  return;
}
