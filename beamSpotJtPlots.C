#include <iostream>
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "commonUtility.h"

TFile* inFileGood_p = 0;
TFile* inFileBad_p = 0;

void makeJtPlots(const char* jetAlg, const char* outName = "spotPlots.root", Bool_t isPrivate = true)
{
  const char* privStand;

  if(isPrivate)
    privStand = "Private";
  else
    privStand = "Standard";

  TTree* inTreeGood_p = 0;
  TTree* inTreeBad_p = 0;

  inTreeGood_p = (TTree*)inFileGood_p->Get(Form("%sJetAnalyzer/t", jetAlg));
  inTreeGood_p->AddFriend("skimanalysis/HltTree");
  inTreeGood_p->AddFriend("hltanalysis/HltTree");

  inTreeBad_p = (TTree*)inFileBad_p->Get(Form("%sJetAnalyzer/t", jetAlg));
  inTreeBad_p->AddFriend("skimanalysis/HltTree");
  inTreeBad_p->AddFriend("hltanalysis/HltTree");

  inFileGood_p->cd();

  inTreeGood_p->Draw("jtpt >> jtPtGood_h(38, 10, 200)", "jtpt < 200 && pcollisionEventSelection && pHBHENoiseFilter && HLT_HIJet80_v1");
  TH1F* getHistPtGood_p = (TH1F*)inFileGood_p->Get("jtPtGood_h");
  inTreeGood_p->Draw("jteta >> jtEtaGood_h(40, -2.0, 2.0)", "TMath::Abs(jteta) < 2.0 && pcollisionEventSelection && pHBHENoiseFilter && HLT_HIJet80_v1");
  TH1F* getHistEtaGood_p = (TH1F*)inFileGood_p->Get("jtEtaGood_h");
  inTreeGood_p->Draw("jtphi >> jtPhiGood_h(40, -3.15, 3.15)", "pcollisionEventSelection && pHBHENoiseFilter && HLT_HIJet80_v1");
  TH1F* getHistPhiGood_p = (TH1F*)inFileGood_p->Get("jtPhiGood_h");

  getHistPtGood_p->Scale(1/getHistPtGood_p->GetEntries());
  getHistPtGood_p->SetMarkerStyle(20);
  getHistPtGood_p->SetYTitle(Form("Event Fraction (%s)", privStand));
  getHistPtGood_p->GetYaxis()->CenterTitle();
  getHistPtGood_p->SetLabelSize(0);
  getHistPtGood_p->SetMaximum(.6);
  getHistPtGood_p->SetMinimum(TMath::Power(10, -5));

  getHistEtaGood_p->Scale(1/getHistEtaGood_p->GetEntries());
  getHistEtaGood_p->SetMarkerStyle(20);
  getHistEtaGood_p->SetYTitle(Form("Event Fraction (%s)", privStand));
  getHistEtaGood_p->GetYaxis()->CenterTitle();
  getHistEtaGood_p->SetLabelSize(0);
  getHistEtaGood_p->SetMaximum(.045);
  getHistEtaGood_p->SetMinimum(0);

  getHistPhiGood_p->Scale(1/getHistPhiGood_p->GetEntries());
  getHistPhiGood_p->SetMarkerStyle(20);
  getHistPhiGood_p->SetYTitle(Form("Event Fraction (%s)", privStand));
  getHistPhiGood_p->GetYaxis()->CenterTitle();
  getHistPhiGood_p->SetLabelSize(0);
  getHistPhiGood_p->SetMaximum(.035);
  getHistPhiGood_p->SetMinimum(0);

  inFileBad_p->cd();

  inTreeBad_p->Draw("jtpt >> jtPtBad_h(38, 10, 200)", "jtpt < 200 && pcollisionEventSelection && pHBHENoiseFilter && HLT_HIJet80_v1");
  TH1F* getHistPtBad_p = (TH1F*)inFileBad_p->Get("jtPtBad_h");
  inTreeBad_p->Draw("jteta >> jtEtaBad_h(40, -2.0, 2.0)", "TMath::Abs(jteta) < 2.0 && pcollisionEventSelection && pHBHENoiseFilter && HLT_HIJet80_v1");
  TH1F* getHistEtaBad_p = (TH1F*)inFileBad_p->Get("jtEtaBad_h");
  inTreeBad_p->Draw("jtphi >> jtPhiBad_h(40, -3.15, 3.15)", "pcollisionEventSelection && pHBHENoiseFilter && HLT_HIJet80_v1");
  TH1F* getHistPhiBad_p = (TH1F*)inFileBad_p->Get("jtPhiBad_h");

  getHistPtBad_p->Scale(1/getHistPtBad_p->GetEntries());
  getHistPtBad_p->SetMarkerStyle(25);
  getHistPtBad_p->SetMarkerColor(1);
  getHistPtBad_p->SetYTitle(Form("Event Fraction (%s)", privStand));
  getHistPtBad_p->GetYaxis()->CenterTitle();
  getHistPtBad_p->SetLabelSize(0);
  getHistPtBad_p->SetMaximum(.6);
  getHistPtBad_p->SetMinimum(TMath::Power(10, -5));

  getHistEtaBad_p->Scale(1/getHistEtaBad_p->GetEntries());
  getHistEtaBad_p->SetMarkerStyle(25);
  getHistEtaBad_p->SetMarkerColor(1);
  getHistEtaBad_p->SetYTitle(Form("Event Fraction (%s)", privStand));
  getHistEtaBad_p->GetYaxis()->CenterTitle();
  getHistEtaBad_p->SetLabelSize(0);
  getHistEtaBad_p->SetMaximum(.045);
  getHistEtaBad_p->SetMinimum(0);

  getHistPhiBad_p->Scale(1/getHistPhiBad_p->GetEntries());
  getHistPhiBad_p->SetMarkerStyle(25);
  getHistPhiBad_p->SetMarkerColor(1);
  getHistPhiBad_p->SetYTitle(Form("Event Fraction (%s)", privStand));
  getHistPhiBad_p->GetYaxis()->CenterTitle();
  getHistPhiBad_p->SetLabelSize(0);
  getHistPhiBad_p->SetMaximum(.035);
  getHistPhiBad_p->SetMinimum(0);

  TLegend* leg = new TLegend(.5, .70, .7, .85, privStand);

  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(.05);
  leg->SetBorderSize(0);

  leg->AddEntry(getHistPtGood_p, "Good Beamspot", "P");
  leg->AddEntry(getHistPtBad_p, "Bad Beamspot", "P");

  TCanvas* jtPtPlots_p = new TCanvas(Form("%sPtPlots_%s_c", jetAlg, privStand), Form("%sPtPlots_%s_c", jetAlg, privStand), 700, 1000);
  jtPtPlots_p->Divide(1, 2, 0, 0);
  jtPtPlots_p->cd(1);
  gPad->SetLogy();
  getHistPtGood_p->DrawCopy("E1");
  getHistPtBad_p->DrawCopy("SAME E1");
  leg->Draw("SAME");

  jtPtPlots_p->cd(2);
  TH1F* divHistPt_p = new TH1F("divHistPt_p", "divHistPt_p", 38, 10, 200);
  divHistPt_p->Divide(getHistPtBad_p, getHistPtGood_p);
  divHistPt_p->SetXTitle(Form("%s Jet p_{T} (GeV/c)", jetAlg));
  divHistPt_p->SetYTitle(Form("Bad/Good (%s)", privStand));
  divHistPt_p->GetYaxis()->CenterTitle();
  divHistPt_p->SetMaximum(1.5);
  divHistPt_p->SetMinimum(.5);
  divHistPt_p->DrawCopy("E1");

  TLatex* cutLabel_p = new TLatex();
  cutLabel_p->SetNDC();
  cutLabel_p->DrawLatex(.25, .90, "pCollisionEventSelection");
  cutLabel_p->DrawLatex(.25, .85, "pHBHENoiseFilter");
  cutLabel_p->DrawLatex(.25, .80, "HLT_HIJet80_v1");

  TLegend* leg2 = new TLegend(.2, .75, .4, .90, privStand);

  leg2->SetFillColor(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(.05);
  leg2->SetBorderSize(0);

  leg2->AddEntry(getHistEtaGood_p, "Good Beamspot", "P");
  leg2->AddEntry(getHistEtaBad_p, "Bad Beamspot", "P");

  TCanvas* jtEtaPlots_p = new TCanvas(Form("%sEtaPlots_%s_c", jetAlg, privStand), Form("%sEtaPlots_%s_c", jetAlg, privStand), 700, 1000);
  jtEtaPlots_p->Divide(1, 2, 0, 0);
  jtEtaPlots_p->cd(1);
  getHistEtaGood_p->DrawCopy("E1");
  getHistEtaBad_p->DrawCopy("SAME E1");
  leg2->Draw("SAME");
  
  jtEtaPlots_p->cd(2);
  TH1F* divHistEta_p = new TH1F("divHistEta_p", "divHistEta_p", 40, -2.0, 2.0);
  divHistEta_p->Divide(getHistEtaBad_p, getHistEtaGood_p);
  divHistEta_p->SetXTitle(Form("%s Jet #eta", jetAlg));
  divHistEta_p->SetYTitle(Form("Bad/Good (%s)", privStand));
  divHistEta_p->GetYaxis()->CenterTitle();
  divHistEta_p->SetMaximum(1.5);
  divHistEta_p->SetMinimum(.5);
  divHistEta_p->DrawCopy("E1");

  cutLabel_p->DrawLatex(.25, .90, "pCollisionEventSelection");
  cutLabel_p->DrawLatex(.25, .85, "pHBHENoiseFilter");
  cutLabel_p->DrawLatex(.25, .80, "HLT_HIJet80_v1");


  TLegend* leg3 = new TLegend(.2, .35, .4, .50, privStand);

  leg3->SetFillColor(0);
  leg3->SetTextFont(42);
  leg3->SetTextSize(.05);
  leg3->SetBorderSize(0);

  leg3->AddEntry(getHistPhiGood_p, "Good Beamspot", "P");
  leg3->AddEntry(getHistPhiBad_p, "Bad Beamspot", "P");

  TCanvas* jtPhiPlots_p = new TCanvas(Form("%sPhiPlots_%s_c", jetAlg, privStand), Form("%sPhiPlots_%s_c", jetAlg, privStand), 700, 1000);
  jtPhiPlots_p->Divide(1, 2, 0, 0);
  jtPhiPlots_p->cd(1);
  getHistPhiGood_p->DrawCopy("E1");
  getHistPhiBad_p->DrawCopy("SAME E1");
  leg3->Draw("SAME");
  
  jtPhiPlots_p->cd(2);
  TH1F* divHistPhi_p = new TH1F("divHistPhi_p", "divHistPhi_p", 40, -3.15, 3.15);
  divHistPhi_p->Divide(getHistPhiBad_p, getHistPhiGood_p);
  divHistPhi_p->SetXTitle(Form("%s Jet #phi", jetAlg));
  divHistPhi_p->SetYTitle(Form("Bad/Good (%s)", privStand));
  divHistPhi_p->GetYaxis()->CenterTitle();
  divHistPhi_p->SetMaximum(1.5);
  divHistPhi_p->SetMinimum(.5);
  divHistPhi_p->DrawCopy("E1");

  cutLabel_p->DrawLatex(.25, .90, "pCollisionEventSelection");
  cutLabel_p->DrawLatex(.25, .85, "pHBHENoiseFilter");
  cutLabel_p->DrawLatex(.25, .80, "HLT_HIJet80_v1");


  TFile* outFile_p = new TFile(outName, "UPDATE");
  jtPtPlots_p->Write();
  claverCanvasSaving(jtPtPlots_p, Form("pngDir/%sPtPlots_%s_c", jetAlg, privStand), "png");
  jtEtaPlots_p->Write();
  claverCanvasSaving(jtEtaPlots_p, Form("pngDir/%sEtaPlots_%s_c", jetAlg, privStand), "png");
  jtPhiPlots_p->Write();
  claverCanvasSaving(jtPhiPlots_p, Form("pngDir/%sPhiPlots_%s_c", jetAlg, privStand), "png");
  outFile_p->Close();

  delete outFile_p;
  delete cutLabel_p;
  delete divHistPt_p;
  delete divHistEta_p;
  delete divHistPhi_p;
  delete leg;
  delete leg2;
  delete leg3;
  delete jtPhiPlots_p;
  delete jtEtaPlots_p;
  delete jtPtPlots_p;

  return;
}



void beamSpotJtPlots(const char* inGoodName, const char* inBadName, const char* outName = "spotPlots.root", Bool_t isPrivate = true)
{
  TH1::SetDefaultSumw2();

  inFileGood_p = new TFile(inGoodName, "READ");
  inFileBad_p = new TFile(inBadName, "READ");

  makeJtPlots("akVs3Calo", outName, isPrivate);
  makeJtPlots("akVs3PF", outName, isPrivate);

  inFileGood_p->Close();
  inFileBad_p->Close();

  delete inFileBad_p;
  delete inFileGood_p;
  return;
}
