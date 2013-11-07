#include "TTree.h"
#include "../gammaJetAnalysis/commonUtility.h"
#include "TDatime.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TProfile.h"

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


void niceTProf(TProfile* uglyTProf, float max , float min, float ndivX, float ndivY)
{
  uglyTProf->GetYaxis()->SetTitleOffset(1.25);
  uglyTProf->GetXaxis()->CenterTitle();
  uglyTProf->GetYaxis()->CenterTitle();
  uglyTProf->SetMaximum(max);
  uglyTProf->SetMinimum(min);
  uglyTProf->SetNdivisions(ndivX);
  uglyTProf->SetNdivisions(ndivY, "Y");
}


void makePtSpectra(TTree* getTree_p, const char* outName, const char* gorr, const char* genTrk, Int_t nBins, Float_t histLow, Float_t histHi, Int_t centLow, Int_t centHi)
{
  inFile_p->cd();

  const char* title = Form("%sPt_%d%d", gorr, (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5));

  TCanvas* ptCanvas_p = new TCanvas(Form("%s_c", title), Form("%s_c", title), 1);
  TH1F* ptHist_p;

  TString name = Form("%sPt >> %s_h(%d, %f, %f)", genTrk, title, nBins, histLow, histHi);
  TCut centCut = Form("hiBin >= %d && hiBin <= %d", centLow, centHi);
  getTree_p->Draw(name, centCut);
  ptHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  TLatex* cent_p;

  if(*gorr == 103)
    cent_p = new TLatex(.6, .3, Form("Truth, %d-%d%%", (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5)));
  else if(*gorr == 114)
    cent_p = new TLatex(.6, .3, Form("Reco, %d-%d%%", (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5)));
  else{
    std::cout << "Error: Neither 'g' or 'r' input." << std::endl;
    return;
  }

  ptHist_p->SetYTitle("Multiplicity");
  ptHist_p->SetXTitle("p_{T} (GeV/c)");
  cent_p->SetNDC();
  cent_p->Draw();

  ptCanvas_p->SetLogy();

  outFile_p = new TFile(outName, "UPDATE");
  ptHist_p->Write();
  ptCanvas_p->Write();
  outFile_p->Close();

  delete outFile_p;
  delete cent_p;
  delete ptCanvas_p;
}


void makeDivHist(const char* fileName, Int_t centLow, Int_t centHi)
{
  TFile* divFile_p = new TFile(fileName, "UPDATE");

  TH1F* numHist_p = (TH1F*)divFile_p->Get(Form("rPt_%d%d_h", centLow, centHi));
  TH1F* denomHist_p = (TH1F*)divFile_p->Get(Form("gPt_%d%d_h", centLow, centHi));

  TCanvas* divCanvas_p = new TCanvas(Form("rDivGPt_%d%d_c", centLow, centHi), Form("rDivGPt_%d%d_c", centLow, centHi), 1);

  numHist_p->Divide(denomHist_p);
  numHist_p->SetYTitle("Mult Fraction");
  numHist_p->Draw();

  divCanvas_p->Write();
  divFile_p->Close();

  delete divCanvas_p;
  delete divFile_p;
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


void addHistToPanel(TFile* file_p, TH1F* hist_p, TCanvas* canv_p, const char* gorr, Int_t centLow, Int_t centHi, Int_t pos)
{
  hist_p = (TH1F*)file_p->Get(Form("%sAsymm_%d%d_h", gorr, centLow, centHi));
  canv_p->cd(pos);
  hist_p->Draw();

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->DrawLatex(.6, .3, Form("%d-%d%%", centLow, centHi));
  
  if(*gorr == 103 && pos == 1)
    label_p->DrawLatex(.4, .85, Form("DiJet Asymmetry, Truth"));
  else if(*gorr == 114 && pos == 1)
    label_p->DrawLatex(.4, .85, Form("DiJet Asymmetry, Reco"));

  delete label_p;
}


void makeAsymmPanel(const char* fileName, const char* gorr)
{
  TFile* panelFile_p = new TFile(fileName, "UPDATE");
  TH1F* asymmHist_p;

  TCanvas* asymmPanel_p = new TCanvas(Form("%sAsymmPanel_c", gorr), Form("%sAsymmPanel_c", gorr), 1);
  asymmPanel_p->Divide(3, 2, 0, 0);

  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 0, 100, 1);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 50, 100, 2);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 30, 50, 3);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 20, 30, 4);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 10, 20, 5);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 0, 10, 6);

  asymmPanel_p->Write();
  panelFile_p->Close();
  delete panelFile_p;
  delete asymmPanel_p;
}


void makeImbHist(TTree* getTree_p, const char* outName, const char* gorr, const char* perpProj, const char* FHL, Int_t nBins, Int_t histLow, Int_t histHi)
{
  inFile_p->cd();

  const char* title = Form("%sImb%s%s", gorr, perpProj, FHL);
  TH1F* imbHist_p;

  TString name = Form("%sImb%s%s >> %s_h(%d, %d, %d)", gorr, perpProj, FHL, title, nBins, histLow, histHi);

  getTree_p->Draw(name);
  imbHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  imbHist_p->SetYTitle("Events");
  imbHist_p->SetXTitle("<#slash{p}_{T}^{||}> (GeV/c)");

  outFile_p = new TFile(outName, "UPDATE");
  imbHist_p->Write();
  outFile_p->Close();

  delete outFile_p;
}


void makeAsymmImbProf(TTree* getTree_p, const char* outName, const char* gorr, const char* perpProj, const char* FHL, Int_t xBins, Int_t yBins, Int_t yLow, Int_t yHi, Int_t centLow, Int_t centHi, Int_t profLow, Int_t profHi, Int_t profXDiv, Int_t profYDiv)
{
  inFile_p->cd();

  const char* title = Form("%sAsymmImb%s%s_%d%d", gorr, perpProj, FHL, (Int_t)(centLow*2.5), (Int_t)((centHi + 1)*2.5));

  TH2F* asymmImbHist_p;
  TProfile* asymmImbHistProf_p;

  TString name = Form("%sImb%s%s:(gLeadJtPt - gSubLeadJtPt)/(gLeadJtPt + gSubLeadJtPt) >> %s_h(%d, 0., 0.5, %d, %d, %d)", gorr, perpProj, FHL, title, xBins, yBins, yLow, yHi);
  TCut centCut = Form("hiBin >= %d && hiBin <= %d", centLow, centHi);

  getTree_p->Draw(name, centCut);
  asymmImbHist_p = (TH2F*)inFile_p->Get(Form("%s_h", title));

  asymmImbHist_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  asymmImbHist_p->SetXTitle("A_{J}");

  asymmImbHistProf_p = (TProfile*)asymmImbHist_p->ProfileX(Form("%s_prof", title));

  niceTProf(asymmImbHistProf_p, profHi, profLow, profXDiv, profYDiv);

  outFile_p = new TFile(outName, "UPDATE");
  asymmImbHist_p->Write();
  asymmImbHistProf_p->Write();
  outFile_p->Close();

  delete outFile_p;
}


void addProfToPanel(TFile* file_p, TProfile* prof_p, TCanvas* canv_p, const char* gorr, const char* perpProj, const char* FHL, Int_t centLow, Int_t centHi, Int_t pos)
{
  prof_p = (TProfile*)file_p->Get(Form("%sAsymmImb%s%s_%d%d_prof", gorr, perpProj, FHL, centLow, centHi));
  canv_p->cd(pos);
  prof_p->Draw();

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->DrawLatex(.6, .3, Form("%d-%d%%", centLow, centHi));

  if(*gorr == 103 && pos == 1)
    label_p->DrawLatex(.2, .85, Form("<#slash{p}_{T}^{||}> as a Function of A_{J}, Truth, %s", perpProj));
  else if(*gorr == 114 && pos == 1)
    label_p->DrawLatex(.2, .85, Form("<#slash{p}_{T}^{||}> as a Function of A_{J}, Reco, %s", perpProj));

  delete label_p;
}


void makeAsymmImbPanel(const char* fileName, const char* gorr, const char* perpProj, const char* FHL)
{
  TFile* panelFile_p = new TFile(fileName, "UPDATE");
  TProfile* getProf_p;

  TCanvas* profPanel_p = new TCanvas(Form("%sAsymmImb%s%sPanel_c", gorr, perpProj, FHL), Form("%sAsymmImb%s%sPanel_c", gorr, perpProj, FHL), 1);
  profPanel_p->Divide(2, 1, 0, 0);

  addProfToPanel(panelFile_p, getProf_p, profPanel_p, gorr, perpProj, FHL, 30, 100, 1);
  TLine *zeroLine_p = new TLine(0., 0., .5, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw();
  addProfToPanel(panelFile_p, getProf_p, profPanel_p, gorr, perpProj, FHL, 0, 30, 2);
  zeroLine_p->Draw();

  profPanel_p->Write();
  claverCanvasSaving(profPanel_p, Form("%sAsymmImb%s%sPanel", gorr, perpProj, FHL), "png");
  panelFile_p->Close();
  delete panelFile_p;
  delete profPanel_p;
  delete zeroLine_p;
}


void cfmDiJetHist(const char* inName, bool montecarlo, const char* outName)
{
  inFile_p = new TFile(inName, "READ");
  inTree_p = (TTree*)inFile_p->Get("jetTree");
  inTree_p->AddFriend("trackTree");

  if(montecarlo)
    inTree_p->AddFriend("genTree");

  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 39);
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 11);
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 12, 39);

  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 0, 39);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 0, 3);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 4, 7);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 8, 11);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 12, 19);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 20, 39);
  
  makeAsymmPanel(outName, "r");

  makeImbHist(inTree_p, outName, "r", "Proj", "F", 21, -210, 210);
  makeImbHist(inTree_p, outName, "r", "Proj", "H", 21, -210, 210);
  makeImbHist(inTree_p, outName, "r", "Proj", "L", 21, -210, 210);

  makeImbHist(inTree_p, outName, "r", "Perp", "F", 21, -210, 210);
  makeImbHist(inTree_p, outName, "r", "Perp", "H", 21, -210, 210);
  makeImbHist(inTree_p, outName, "r", "Perp", "L", 21, -210, 210);

  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "H", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "H", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "H", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "L", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "L", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "L", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);

  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "H", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "H", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "H", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "L", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "L", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "L", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);

  makeAsymmImbPanel(outName, "r", "Proj", "F");
  makeAsymmImbPanel(outName, "r", "Proj", "H");
  makeAsymmImbPanel(outName, "r", "Proj", "L");

  makeAsymmImbPanel(outName, "r", "Perp", "F");
  makeAsymmImbPanel(outName, "r", "Perp", "H");
  makeAsymmImbPanel(outName, "r", "Perp", "L");


  //Monte

  if(montecarlo){
    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 0, 39);
    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 0, 11);
    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 12, 39);

    makeDivHist(outName, 0, 100);
    makeDivHist(outName, 0, 30);
    makeDivHist(outName, 30, 100);

    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 0, 39);
    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 0, 3);
    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 4, 7);
    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 8, 11);
    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 12, 19);
    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 20, 39);

    makeAsymmPanel(outName, "g");

    makeImbHist(inTree_p, outName, "g", "Proj", "F", 21, -210, 210);
    makeImbHist(inTree_p, outName, "g", "Proj", "H", 21, -210, 210);
    makeImbHist(inTree_p, outName, "g", "Proj", "L", 21, -210, 210);

    makeImbHist(inTree_p, outName, "g", "Perp", "F", 21, -210, 210);
    makeImbHist(inTree_p, outName, "g", "Perp", "H", 21, -210, 210);
    makeImbHist(inTree_p, outName, "g", "Perp", "L", 21, -210, 210);

    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "F", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "H", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "H", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "H", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "L", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "L", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "L", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);

    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "F", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "H", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "H", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "H", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "L", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "L", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404);
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "L", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404);

    makeAsymmImbPanel(outName, "g", "Proj", "F");
    makeAsymmImbPanel(outName, "g", "Proj", "H");
    makeAsymmImbPanel(outName, "g", "Proj", "L");

    makeAsymmImbPanel(outName, "g", "Perp", "F");
    makeAsymmImbPanel(outName, "g", "Perp", "H");
    makeAsymmImbPanel(outName, "g", "Perp", "L");
  }

  inFile_p->Close();
  delete inFile_p;
}


