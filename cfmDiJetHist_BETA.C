#include "../gammaJetAnalysis/commonUtility.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TDatime.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TProfile.h"

TFile* inFile_p = 0;
TFile* outFile_p = 0;

TTree* inTree_p = 0;

//append to every histo so know which sample via shorthand on workblog
const char* fileTag;

//shorthands

//auto ln -s title of address /net/hisrv0001/home/yenjie/scratch/tmp/Pythia80_HydjetDrum_mix01_HiForest2_v20.root
const char* Di80a = "Pythia80_HydjetDrum_mix01_HiForest2_v20.root";

//auto ln -s title of address /mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet80_HydjetDrum_v27_mergedV1.root
const char* Di80b = "Dijet80_HydjetDrum_v27_mergedV1.root";

//auto ln -s title of address /mnt/hadoop/cms/store/user/luck/PbPb_pythiaHYDJET_forest_EmEnrichedDijet/PbPb_pythiaHYDJET_forest_EmEnrichedDijet80.root
const char* EmDi80a = "PbPb_pythiaHYDJET_forest_EmEnrichedDijet80.root";


Float_t getDPHI( Float_t phi1, Float_t phi2) {
  Float_t dphi = phi1 - phi2;

  if ( dphi > TMath::Pi())
    dphi = dphi - 2.*(TMath::Pi());
  if ( dphi <= -(TMath::Pi()) )
    dphi = dphi + 2.*(TMath::Pi());

  if ( TMath::Abs(dphi) > TMath::Pi()) {
    std::cout << " commonUtility::getDPHI error!!! dphi is bigger than TMath::Pi() " << std::endl;
  }

  return dphi;
}

Float_t getAbsDphi( Float_t phi1, Float_t phi2) {
  return TMath::Abs( getDPHI(phi1, phi2) ) ;
}


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


void makePtSpectra(TTree* getTree_p, const char* outName, const char* gorr, const char* genTrk, Int_t nBins, Float_t histLow, Float_t histHi, Int_t centLow, Int_t centHi, Float_t delPhiLow, Float_t delPhiHi)
{
  inFile_p->cd();

  const char* delPhiTitle;

  if(delPhiHi < TMath::PiOver2() + .1)
    delPhiTitle = "Lead";
  else if(delPhiLow > TMath::PiOver2() - .1)
    delPhiTitle = "Away";
  else
    delPhiTitle = "All";


  const char* title = Form("%sPt_%d%d_%s_%s", gorr, (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5), delPhiTitle, fileTag);

  TCanvas* ptCanvas_p = new TCanvas(Form("%s_c", title), Form("%s_c", title), 1);
  TH1F* ptHist_p;

  TString name = Form("%s_h(%d, %f, %f)", title, nBins, histLow, histHi);
  TCut centCut = Form("hiBin >= %d && hiBin <= %d && %f < %sLeadDelPhi && %sLeadDelPhi < %f", centLow, centHi, delPhiLow, genTrk, genTrk, delPhiHi);
  //getTree_p->Draw(name, centCut);
  getTree_p->Project(name, Form("%sPt", genTrk), centCut);
  ptHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  TLatex* cent_p;

  if(*gorr == 103)
    cent_p = new TLatex(.6, .3, Form("Truth, %d-%d%%", (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5)));
  else if(*gorr == 114)
    cent_p = new TLatex(.6, .3, Form("Tracks, %d-%d%%", (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5)));
  else{
    std::cout << "Error: Neither 'g' or 'r' input." << std::endl;
    return;
  }

  ptHist_p->SetYTitle("Mult/1 GeV");
  ptHist_p->SetXTitle("p_{T} (GeV/c)");
  ptCanvas_p->SetLogy();

  outFile_p = new TFile(outName, "UPDATE");
  ptHist_p->Write();
  handsomeTH1N(ptHist_p);
  ptHist_p->Draw();
  cent_p->SetNDC();
  cent_p->Draw();
  ptCanvas_p->Write();
  outFile_p->Close();

  delete outFile_p;
  delete cent_p;
  delete ptCanvas_p;
}


void makeROnGHist(const char* fileName, Int_t centLow, Int_t centHi, const char* lAAll)
{
  outFile_p = new TFile(fileName, "UPDATE");

  TCanvas* rOnGCanvas_p = new TCanvas(Form("rOnGPt_%d%d_%s_%s_c", centLow, centHi, lAAll, fileTag), Form("rOnGPt_%d%d_%s_%s_c", centLow, centHi, lAAll, fileTag), 1);
  TH1F* rHist_p = (TH1F*)outFile_p->Get(Form("rPt_%d%d_%s_%s_h", centLow, centHi, lAAll, fileTag));
  TH1F* gHist_p = (TH1F*)outFile_p->Get(Form("gPt_%d%d_%s_%s_h", centLow, centHi, lAAll, fileTag));

  gHist_p->SetMarkerColor(kBlue);
  rHist_p->SetMarkerColor(kRed);
  rHist_p->SetMarkerStyle(3);

  gHist_p->Draw("P");
  rHist_p->Draw("SAME P");

  rOnGCanvas_p->SetLogy();

  
  TLegend *leg = new TLegend(0.55, 0.45, 0.9, 0.65, Form("%d-%d%%, %s", centLow, centHi, lAAll));
  leg->SetFillColor(0);
  leg->AddEntry(gHist_p, "Truth", "p");
  leg->AddEntry(rHist_p, "Tracks", "p");
  leg->Draw("SAME");
  

  rOnGCanvas_p->Write();
  claverCanvasSaving(rOnGCanvas_p, Form("pngDir/rOnGPt_%d%d_%s_%s", centLow, centHi, lAAll, fileTag), "png");
  outFile_p->Close();
  delete leg;
  delete rOnGCanvas_p;
  delete outFile_p;
}


void makeLDivAHist(const char* fileName, const char* gorr, Int_t centLow, Int_t centHi)
{
  outFile_p = new TFile(fileName, "UPDATE");

  TH1F* leadHist_p = (TH1F*)outFile_p->Get(Form("%sPt_%d%d_%s_%s_h", gorr, centLow, centHi, "Lead", fileTag));
  TH1F* awayHist_p = (TH1F*)outFile_p->Get(Form("%sPt_%d%d_%s_%s_h", gorr, centLow, centHi, "Away", fileTag));

  handsomeTH1N(leadHist_p);
  handsomeTH1N(awayHist_p);

  leadHist_p->Divide(awayHist_p);
  leadHist_p->Write(Form("%sLDivAPt_%d%d_%s_h", gorr, centLow, centHi, fileTag));
  outFile_p->Close();
  delete outFile_p;
}


//This is no good... FIX

void makeRDivGHist(const char* fileName, Int_t centLow, Int_t centHi, const char* lAAll)
{
  outFile_p = new TFile(fileName, "UPDATE");

  TCanvas* divCanvas_p = new TCanvas();
  TH1F* hEmpty_p = new TH1F(Form("hEmpty_%d%d_%s_%s_h", centLow, centHi, lAAll, fileTag), "hEmpty", 20, -0.5, 19.5);
  hEmpty_p->SetXTitle("p_{T} (GeV/c)");
  hEmpty_p->GetXaxis()->CenterTitle();
  hEmpty_p->SetYTitle("Mult Fraction");
  hEmpty_p->GetYaxis()->CenterTitle();
  hEmpty_p->Draw();

  TH1F* numHist_p = (TH1F*)outFile_p->Get(Form("rPt_%d%d_%s_%s_h", centLow, centHi, lAAll, fileTag));
  TH1F* denomHist_p = (TH1F*)outFile_p->Get(Form("gPt_%d%d_%s_%s_h", centLow, centHi, lAAll, fileTag));

  TGraphAsymmErrors* divHist_p = new TGraphAsymmErrors(numHist_p, denomHist_p, "cl=.683 b(1,1) mode");
  divHist_p->SetMarkerColor(1);
  divHist_p->DrawClone("same P");

  divHist_p->Write(Form("rDivGPt_%d%d_%s_%s_h", centLow, centHi, lAAll, fileTag));
  outFile_p->Close();

  delete divCanvas_p;
  delete outFile_p;
}


void makeAsymmHist(TTree* getTree_p, const char* outName, const char* gorr, Int_t nBins, Int_t histLow, Int_t histHi, Int_t centLow, Int_t centHi)
{
  inFile_p->cd();

  const char* title = Form("%sAsymm_%d%d_%s", gorr, (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5), fileTag);

  TCanvas* asymmCanvas_p = new TCanvas(Form("%s_c", title), Form("%s_c", title), 1);
  TH1F* asymmHist_p;

  TString name = Form("(%sLeadJtPt - %sSubLeadJtPt)/(%sLeadJtPt + %sSubLeadJtPt) >> %s_h(%d, %d, %d)", gorr, gorr, gorr, gorr, title, nBins, histLow, histHi);
  TCut centCut = Form("hiBin >= %d && hiBin <= %d && %sLeadJtPt > 120 && %sSubLeadJtPt > 50", centLow, centHi, gorr, gorr);
  getTree_p->Draw(name, centCut);
  asymmHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  TLatex* cent_p;

  if(*gorr == 103){
    niceTH1(asymmHist_p, .6, 0., 405, 506);
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
  asymmHist_p->Write(Form("%s_h", title));
  asymmCanvas_p->Write();
  outFile_p->Close();

  delete outFile_p;
  delete cent_p;
  delete asymmCanvas_p;
}


void addHistToPanel(TFile* file_p, TH1F* hist_p, TCanvas* canv_p, const char* gorr, Int_t centLow, Int_t centHi, Int_t pos)
{
  hist_p = (TH1F*)file_p->Get(Form("%sAsymm_%d%d_%s_h", gorr, centLow, centHi, fileTag));
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

  TCanvas* asymmPanel_p = new TCanvas(Form("%sAsymmPanel_%s_c", gorr, fileTag), Form("%sAsymmPanel_%s_c", gorr, fileTag), 1);
  asymmPanel_p->Divide(3, 2, 0, 0);

  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 0, 100, 1);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 50, 100, 2);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 30, 50, 3);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 20, 30, 4);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 10, 20, 5);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, 0, 10, 6);

  asymmPanel_p->Write();
  //  claverCanvasSaving(asymmPanel_p, Form("pngDir/%sAsymmPanel_%s", gorr, fileTag), "png");
  panelFile_p->Close();
  delete panelFile_p;
  delete asymmPanel_p;
}


void makeImbHist(TTree* getTree_p, const char* outName, const char* gorr, const char* perpProj, const char* FHL, Int_t nBins, Int_t histLow, Int_t histHi)
{
  inFile_p->cd();

  const char* title = Form("%sImb%s%s_%s", gorr, perpProj, FHL, fileTag);
  TH1F* imbHist_p;

  TString var = Form("%sImb%s%s", gorr, perpProj, FHL);
  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);

  getTree_p->Project(name, var);
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

  const char* title = Form("%sAsymmImb%s%s_%d%d_%s", gorr, perpProj, FHL, (Int_t)(centLow*2.5), (Int_t)((centHi + 1)*2.5), fileTag);

  TH2F* asymmImbHist_p;
  TProfile* asymmImbHistProf_p;

  TString var = Form("%sImb%s%s:(gLeadJtPt - gSubLeadJtPt)/(gLeadJtPt + gSubLeadJtPt)", gorr, perpProj, FHL);
  TString name = Form("%s_h(%d, 0., 0.5, %d, %d, %d)", title, xBins, yBins, yLow, yHi);
  TCut centCut = Form("hiBin >= %d && hiBin <= %d", centLow, centHi);

  getTree_p->Project(name, var, centCut);
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
  prof_p = (TProfile*)file_p->Get(Form("%sAsymmImb%s%s_%d%d_%s_prof", gorr, perpProj, FHL, centLow, centHi, fileTag));
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

  TCanvas* profPanel_p = new TCanvas(Form("%sAsymmImb%s%sPanel_%s_c", gorr, perpProj, FHL, fileTag), Form("%sAsymmImb%s%sPanel_%s_c", gorr, perpProj, FHL, fileTag), 1);
  profPanel_p->Divide(2, 1, 0, 0);

  addProfToPanel(panelFile_p, getProf_p, profPanel_p, gorr, perpProj, FHL, 30, 100, 1);
  TLine *zeroLine_p = new TLine(0., 0., .5, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw();
  addProfToPanel(panelFile_p, getProf_p, profPanel_p, gorr, perpProj, FHL, 0, 30, 2);
  zeroLine_p->Draw();

  profPanel_p->Write();
  if(*FHL == 70) 
    claverCanvasSaving(profPanel_p, Form("./pngDir/%sAsymmImb%s%sPanel_%s", gorr, perpProj, FHL, fileTag), "png");

  panelFile_p->Close();
  delete panelFile_p;
  delete profPanel_p;
  delete zeroLine_p;
}


void makeAsymmImbSamePanel(const char* fileName, const char* perpProj, const char* FHL)
{
  TFile* panelFile_p = new TFile(fileName, "UPDATE");
  TProfile* getProf_p;

  TCanvas* profPanel_p = new TCanvas(Form("gRAsymmImb%s%sPanel_%s_c", perpProj, FHL, fileTag), Form("gRAsymmImb%s%sPanel_%s_c", perpProj, FHL, fileTag), 1);
  profPanel_p->Divide(2, 1, 0, 0);

  getProf_p = (TProfile*)panelFile_p->Get(Form("gAsymmImb%s%s_30100_%s_prof", perpProj, FHL, fileTag));
  profPanel_p->cd(1);
  getProf_p->SetMarkerColor(kBlue);
  getProf_p->SetLineColor(kBlue);
  getProf_p->SetLineWidth(1);

  getProf_p->Draw();

  /*
  getProf_p = (TProfile*)panelFile_p->Get(Form("rAsymmImb%s%s_30100_prof", perpProj, FHL));
  profPanel_p->cd(1);
  getProf_p->SetMarkerColor(kRed);
  getProf_p->SetMarkerStyle(3);
  getProf_p->SetLineColor(kRed);
  getProf_p->SetLineWidth(1);

  getProf_p->Draw("SAME");
  */

  TLine *zeroLine_p = new TLine(0., 0., .5, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw();

  //  zeroLine_p->Draw();

  profPanel_p->Write();
  //  if(*FHL == 70)
  // claverCanvasSaving(profPanel_p, Form("./pngDir/gRAsymmImb%s%sPanel", perpProj, FHL), "png");

  panelFile_p->Close();
  delete panelFile_p;
  delete profPanel_p;
  delete zeroLine_p;
}


void cfmDiJetHist_BETA(const char* inName = "inFile_CFMHIST_BETA.root", bool montecarlo = 0, const char* outName = "outFile_CFMHIST_BETA.root")
{
  if(strcmp(inName, Di80a)){
    std::cout << Di80a << std::endl;
    fileTag = "Di80a";
  }
  else if(strcmp(inName, Di80b)){
    std::cout << Di80b << std::endl;
    fileTag = "Di80b";
  }
  else if(strcmp(inName, EmDi80a)){
    std::cout << EmDi80a << std::endl;
    fileTag = "EmDi80a";
  }

  std::cout << fileTag << std::endl;

  inFile_p = new TFile(inName, "READ");
  inTree_p = (TTree*)inFile_p->Get("jetTree");
  inTree_p->AddFriend("trackTree");

  if(montecarlo)
    inTree_p->AddFriend("genTree");

  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 39, 0., TMath::PiOver2());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 11, 0., TMath::PiOver2());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 12, 39, 0., TMath::PiOver2());

  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 39, TMath::PiOver2(), TMath::Pi());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 11, TMath::PiOver2(), TMath::Pi());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 12, 39, TMath::PiOver2(), TMath::Pi());

  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 39, 0., TMath::Pi());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 11, 0., TMath::Pi());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 12, 39, 0., TMath::Pi());

  makeLDivAHist(outName, "r", 0, 100);
  makeLDivAHist(outName, "r", 0, 30);
  makeLDivAHist(outName, "r", 30, 100);

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
    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 0, 39, 0., TMath::PiOver2());
    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 0, 11, 0., TMath::PiOver2());
    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 12, 39, 0., TMath::PiOver2());

    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 0, 39, TMath::PiOver2(), TMath::Pi());
    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 0, 11, TMath::PiOver2(), TMath::Pi());
    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 12, 39, TMath::PiOver2(), TMath::Pi());

    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 0, 39, 0., TMath::Pi());
    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 0, 11, 0., TMath::Pi());
    makePtSpectra(inTree_p, outName, "g", "gen", 20, -.5, 19.5, 12, 39, 0., TMath::Pi());

    makeLDivAHist(outName, "g", 0, 100);
    makeLDivAHist(outName, "g", 0, 30);
    makeLDivAHist(outName, "g", 30, 100);

    makeROnGHist(outName, 0, 100, "All");
    makeROnGHist(outName, 0, 30, "All");
    makeROnGHist(outName, 30, 100, "All");

    makeROnGHist(outName, 0, 100, "Lead");
    makeROnGHist(outName, 0, 30, "Lead");
    makeROnGHist(outName, 30, 100, "Lead");

    makeROnGHist(outName, 0, 100, "Away");
    makeROnGHist(outName, 0, 30, "Away");
    makeROnGHist(outName, 30, 100, "Away");

    makeRDivGHist(outName, 0, 100, "All");
    makeRDivGHist(outName, 0, 30, "All");
    makeRDivGHist(outName, 30, 100, "All");

    makeRDivGHist(outName, 0, 100, "Lead");
    makeRDivGHist(outName, 0, 30, "Lead");
    makeRDivGHist(outName, 30, 100, "Lead");

    makeRDivGHist(outName, 0, 100, "Away");
    makeRDivGHist(outName, 0, 30, "Away");
    makeRDivGHist(outName, 30, 100, "Away");

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

    makeAsymmImbSamePanel(outName, "Proj", "F");

  }

  inFile_p->Close();
  delete inFile_p;
}


