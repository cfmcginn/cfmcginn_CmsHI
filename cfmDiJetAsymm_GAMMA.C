#include "TStyle.h"
#include "commonUtility.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TDatime.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TProfile.h"

TFile* inFile_p = 0;
TFile* inFile2_p = 0;
TFile* outFile_p = 0;

TTree* inTree_p = 0;
TTree* inTree2_p = 0;


//append to every histo so know which sample via shorthand on workblog                                                          

const char* fileTag;

//shorthands, w/ _CFMSKIM.h
                                                                                                                                  
const char* Di80a = "Pythia80_HydjetDrum_mix01_HiForest2_v20_CFMSKIM.root";
                                                                                                                                  
const char* Di80b = "Dijet80_HydjetDrum_v27_mergedV1_CFMSKIM.root";
 
const char* Di80c = "HydjetDrum_Pyquen_Dijet80_Embedded_d20140122_Track7_v2_CFMSKIM.root";

const char* Di80d = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";
                                                                                                                                 
const char* Di100a = "Dijet100_HydjetDrum_v27_mergedV1_CFMSKIM.root";
                                                                                                                                 
const char* EmDi80a = "PbPb_pythiaHYDJET_forest_EmEnrichedDijet80_CFMSKIM.root";

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


Bool_t sameSign(Double_t num1, Double_t num2){
  Bool_t sameSign = false;

  if((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0))
    sameSign = true;

  return sameSign;
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

  uglyTProf->SetMarkerColor(1);
  uglyTProf->SetMarkerSize(1);
  uglyTProf->SetMarkerStyle(20);
  uglyTProf->SetLineColor(1);
}


void niceTGraphErrors(TGraphErrors* uglyTGraph, float max, float min)
{
  uglyTGraph->GetYaxis()->SetTitleOffset(1.25);
  uglyTGraph->GetXaxis()->CenterTitle();
  uglyTGraph->GetYaxis()->CenterTitle();
  uglyTGraph->SetMaximum(max);
  uglyTGraph->SetMinimum(min);

  uglyTGraph->SetMarkerColor(1);
  uglyTGraph->SetMarkerSize(1);
  uglyTGraph->SetMarkerStyle(20);
  uglyTGraph->SetLineColor(1);
}





TCut makeSetCut(const char* PFCaloT)
{
  TCut setCut = "";
  if(strcmp(PFCaloT, "T") == 0)
    setCut = Form("truthSet == 1");
  else if(strcmp(PFCaloT, "PF") == 0)
    setCut = Form("recoPFSet == 1");
  else if(strcmp(PFCaloT, "Calo") == 0)
    setCut = Form("recoCaloSet == 1");
  else
    std::cout << "makeSetCut: Cut unspecified, returning blank cut" << std::endl;

  return setCut;
}


TCut makeCentCut(Int_t centLow, Int_t centHi)
{
  TCut centCut = "";
  if(centLow >= 0 && centHi >= centLow && centHi <= 39)
    centCut = Form("hiBin >= %d && hiBin <= %d", centLow, centHi);
  else
    std::cout << "makeCentCut: centLow/centHi incorrectly specified, returning blank cut" << std::endl;

  return centCut;
}


TCut makeAsymmCut(const char* PFCaloT, Float_t asymmLow, Float_t asymmHi)
{
  TCut asymmCut = "";
  const char* asymmJetType = Form("%sJtAsymm", PFCaloT);

  if(asymmLow >= 0 && asymmHi >= asymmLow && asymmHi <= 1)
    asymmCut = Form("%s > %f && %s < %f ", asymmJetType, asymmLow, asymmJetType, asymmHi);
  else
    std::cout << "makeAsymmCut: asymmLow/asymmHi incorrectly specified, returning blank cut" << std::endl;

  return asymmCut;
}


TCut makeEtaCut(const char* PFCaloT, const char* GLN)
{
  TCut etaCut = "";
  if(strcmp(GLN, "N") == 0)
    return etaCut;

  const char* leadJt = Form("%sLeadJtEta", PFCaloT);
  const char* subLeadJt = Form("%sSubLeadJtEta", PFCaloT);

  if(strcmp(GLN, "G") == 0){
    etaCut = Form("TMath::Abs(%s) > 1.0 || TMath::Abs(%s) > 1.0", leadJt, subLeadJt);
  }
  else if(strcmp(GLN, "L") == 0){
    etaCut = Form("TMath::Abs(%s) < 1.0 && TMath::Abs(%s) < 1.0", leadJt, subLeadJt);
  }
  else
    std::cout << "makeEtaCut: GLN incorrectly specified, returning blank cut" << std::endl;

  return etaCut;
}


void makeAsymmHist(TTree* getTree_p, const char* outName, const char* PFCaloT, Int_t nBins, Int_t histLow, Int_t histHi, Int_t centLow, Int_t centHi)
{
  inFile_p->cd();

  const char* title = Form("%sAsymm_%d%d_%s", PFCaloT, (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5), fileTag);

  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(PFCaloT);
  TCut centCut = makeCentCut(centLow, centHi);

  getTree_p->Project(name, Form("%sJtAsymm", PFCaloT), centCut && setCut);
  asymmHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  niceTH1(asymmHist_p, .6, 0., 405, 506);

  asymmHist_p->SetYTitle("Event Fraction");
  asymmHist_p->SetXTitle("A_{J} = (p_{T,1} - p_{T,2})/(p_{T,1} + p_{T,2})");

  outFile_p = new TFile(outName, "UPDATE");
  asymmHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
}


void addHistToPanel(TFile* file_p, TCanvas* canv_p, const char* PFCaloT, Int_t centLow, Int_t centHi, Int_t pos)
{
  TH1F* hist_p = (TH1F*)file_p->Get(Form("%sAsymm_%d%d_%s_h", PFCaloT, centLow, centHi, fileTag));
  canv_p->cd(pos);
  hist_p->Draw();

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->DrawLatex(.6, .3, Form("%d-%d%%", centLow, centHi));

  if(strcmp(PFCaloT, "T") == 0 && pos == 1){
    label_p->DrawLatex(.3, .85, Form("Truth DiJet Asymmetry, Truth Set"));
    label_p->DrawLatex(.3, .8, Form("%s", fileTag));
  }
  else if(strcmp(PFCaloT, "PF") == 0 && pos == 1){
    label_p->DrawLatex(.3, .85, Form("RecoPF DiJet Asymmetry, RecoPF Set"));
    label_p->DrawLatex(.3, .8, Form("%s", fileTag));
  }
  else if(strcmp(PFCaloT, "Calo") == 0 && pos == 1){
    label_p->DrawLatex(.3, .85, Form("RecoCalo DiJet Asymmetry, RecoCalo Set"));
    label_p->DrawLatex(.3, .8, Form("%s", fileTag));
  }

  delete label_p;

}


void makeAsymmPanel(const char* fileName, const char* PFCaloT)
{
  TFile* panelFile_p = new TFile(fileName, "UPDATE");
  TCanvas* asymmPanel_p = new TCanvas(Form("%sAsymmPanel_%s_c", PFCaloT, fileTag), Form("%sAsymmPanel_%s_c", PFCaloT, fileTag), 1);
  asymmPanel_p->Divide(3, 2, 0, 0);

  addHistToPanel(panelFile_p, asymmPanel_p, PFCaloT, 0, 100, 1);
  addHistToPanel(panelFile_p, asymmPanel_p, PFCaloT, 50, 100, 2);
  addHistToPanel(panelFile_p, asymmPanel_p, PFCaloT, 30, 50, 3);
  addHistToPanel(panelFile_p, asymmPanel_p, PFCaloT, 20, 30, 4);
  addHistToPanel(panelFile_p, asymmPanel_p, PFCaloT, 10, 20, 5);
  addHistToPanel(panelFile_p, asymmPanel_p, PFCaloT, 0, 10, 6);

  asymmPanel_p->Write();
  panelFile_p->Close();
  delete panelFile_p;
  delete asymmPanel_p;
}


void makePtProjHist(TTree* getTree_p, const char* outName, const char* genTrk, const char* PFCaloT, Int_t ptLow, Int_t ptHi, Int_t centLow, Int_t centHi, const char* symmAsymm)
{
  inFile_p->cd();

  TH1F* ptProjHist_p;

  const char* title = Form("%sPt%s_%d_%d_%d%d_%s_%s_h", genTrk, PFCaloT, ptLow, ptHi, (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5), symmAsymm, fileTag);

  TCut setCut = makeSetCut(PFCaloT);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut asymmCut = "";

  if(strcmp(symmAsymm, "Symm") == 0)
    asymmCut = makeAsymmCut(PFCaloT, 0., 0.1);
  else if(strcmp(symmAsymm, "Asymm") == 0)
    asymmCut = makeAsymmCut(PFCaloT, 0.3, 1.0);
  else{
    std::cout << "makePtProjHist: symmAsymm incorrectly specified, returning w/o filling histogram" << std::endl;
    return;
  }

  TCut ptCut = Form("%sPt < %d && %sPt > %d", genTrk, ptHi, genTrk, ptLow);

  getTree_p->Project(Form("%s(100)", title), Form("%sPt%s", genTrk, PFCaloT), centCut && setCut && ptCut && asymmCut);
  ptProjHist_p = (TH1F*)inFile_p->Get(title);

  handsomeTH1(ptProjHist_p);
  ptProjHist_p->SetXTitle(Form("%sPt Proj. onto Lead %s Jet (GeV/c)", genTrk, PFCaloT));
  ptProjHist_p->SetYTitle(Form("N_{%s}", genTrk));

  outFile_p = new TFile(outName, "UPDATE");
  ptProjHist_p->Write();
  outFile_p->Close();
  delete outFile_p;
}


void makeImbAsymmGraph(TTree* getTree_p, const char* outName, const char* gorr, const char* PFCaloT, const char* perpProj, const char* FHLPT, Int_t centLow, Int_t centHi, Int_t graphLow, Int_t graphHi, const char* GLN = "N", const char* Corr = "")
{
  inFile_p->cd();

  const char* title = Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, FHLPT, Corr, (Int_t)(centLow*2.5), (Int_t)((centHi + 1)*2.5), GLN, fileTag);

  TGraphErrors* imbAsymmGraph_p = new TGraphErrors(4);
  imbAsymmGraph_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTGraphErrors(imbAsymmGraph_p, graphHi, graphLow);

  TH1F* getHist_p;

  TString var = Form("%s%sImb%s%s%s", gorr, PFCaloT, perpProj, FHLPT, Corr);

  TCut setCut = makeSetCut(PFCaloT);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(PFCaloT, GLN);
  TCut asymmCut = makeAsymmCut(PFCaloT, .00, .10);

  getTree_p->Project("0_1", var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get("0_1");

  imbAsymmGraph_p->SetPoint(0, 0.05, getHist_p->GetMean());
  imbAsymmGraph_p->SetPointError(0, 0.05, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));

  asymmCut = makeAsymmCut(PFCaloT, .10, .20);
  getTree_p->Project("1_2", var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get("1_2");
  imbAsymmGraph_p->SetPoint(1, 0.15, getHist_p->GetMean());
  imbAsymmGraph_p->SetPointError(1, 0.05, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));

  asymmCut = makeAsymmCut(PFCaloT, .20, .30);
  getTree_p->Project("2_3", var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get("2_3");
  imbAsymmGraph_p->SetPoint(2, 0.25, getHist_p->GetMean());
  imbAsymmGraph_p->SetPointError(2, 0.05, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));

  asymmCut = makeAsymmCut(PFCaloT, .30, .50);
  getTree_p->Project("3_5", var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get("3_5");
  imbAsymmGraph_p->SetPoint(3, 0.40, getHist_p->GetMean());
  imbAsymmGraph_p->SetPointError(3, 0.10, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));

  imbAsymmGraph_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTGraphErrors(imbAsymmGraph_p, graphHi, graphLow);

  outFile_p = new TFile(outName, "UPDATE");
  imbAsymmGraph_p->Write(title);
  outFile_p->Close();

  delete outFile_p;
  delete imbAsymmGraph_p;
}


Double_t sumYForPTStack(Double_t dIn = 0, Double_t comp1 = 0, Double_t comp2 = 0, Double_t comp3 = 0, Double_t comp4 = 0)
{
  Double_t dOut = dIn;

  if(sameSign(comp1, dOut))
    dOut += comp1;

  if(sameSign(comp2, dOut))
    dOut += comp2;

  if(sameSign(comp3, dOut))
    dOut += comp3;

  if(sameSign(comp4, dOut))
    dOut += comp4;

  return dOut;
}


void makeHistForPTStack(TGraph* g0_1_p, TGraph* g1_2_p, TGraph* g2_4_p, TGraph* g4_8_p, TGraph* g8_100_p, TGraph* gF_p, TH1F* h0_1_p, TH1F* h1_2_p, TH1F* h2_4_p, TH1F* h4_8_p, TH1F* h8_100_p, TH1F* hF_p, const char* CP = "C")
{
  Double_t x0_1[4] = {0, 0, 0, 0};
  Double_t y0_1[4] = {0, 0, 0, 0};
  Double_t x1_2[4] = {0, 0, 0, 0};
  Double_t y1_2[4] = {0, 0, 0, 0};
  Double_t x2_4[4] = {0, 0, 0, 0};
  Double_t y2_4[4] = {0, 0, 0, 0};
  Double_t x4_8[4] = {0, 0, 0, 0};
  Double_t y4_8[4] = {0, 0, 0, 0};
  Double_t x8_100[4] = {0, 0, 0, 0};
  Double_t y8_100[4] = {0, 0, 0, 0};
  Double_t xF[4] = {0, 0, 0, 0};
  Double_t yF[4] = {0, 0, 0, 0};

  for(Int_t iter = 0; iter < 4; iter++){
    g0_1_p->GetPoint(iter, x0_1[iter], y0_1[iter]);
    g1_2_p->GetPoint(iter, x1_2[iter], y1_2[iter]);
    g2_4_p->GetPoint(iter, x2_4[iter], y2_4[iter]);
    g4_8_p->GetPoint(iter, x4_8[iter], y4_8[iter]);
    g8_100_p->GetPoint(iter, x8_100[iter], y8_100[iter]);
    gF_p->GetPoint(iter, xF[iter], yF[iter]);

    h8_100_p->SetBinContent(iter + 1, y8_100[iter]);
    h8_100_p->SetBinError(iter + 1, g8_100_p->GetErrorY(iter));

    h4_8_p->SetBinContent(iter + 1, sumYForPTStack(y4_8[iter], y8_100[iter]));
    h4_8_p->SetBinError(iter + 1, g4_8_p->GetErrorY(iter));

    h2_4_p->SetBinContent(iter + 1, sumYForPTStack(y2_4[iter], y4_8[iter], y8_100[iter]));
    h2_4_p->SetBinError(iter + 1, g2_4_p->GetErrorY(iter));

    h1_2_p->SetBinContent(iter + 1, sumYForPTStack(y1_2[iter], y2_4[iter], y4_8[iter], y8_100[iter]));
    h1_2_p->SetBinError(iter + 1, g1_2_p->GetErrorY(iter));

    h0_1_p->SetBinContent(iter + 1, sumYForPTStack(y0_1[iter], y1_2[iter], y2_4[iter], y4_8[iter], y8_100[iter]));
    h0_1_p->SetBinError(iter + 1, g0_1_p->GetErrorY(iter));

    hF_p->SetBinContent(iter + 1, yF[iter]);
    hF_p->SetBinError(iter + 1, gF_p->GetErrorY(iter));
  }

  h8_100_p->SetXTitle("A_{J}");
  h4_8_p->SetXTitle("A_{J}");
  h2_4_p->SetXTitle("A_{J}");
  h1_2_p->SetXTitle("A_{J}");
  h0_1_p->SetXTitle("A_{J}");
  hF_p->SetXTitle("A_{J}");

  if(strcmp(CP, "P") == 0){
    h8_100_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h4_8_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h2_4_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h1_2_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h0_1_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    hF_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  }
}


void drawHistToPTStack(TH1F* drawHist_p, Int_t color, const char* drawOpt)
{
  drawHist_p->SetFillColor(color);
  drawHist_p->SetMarkerStyle(6);
  drawHist_p->SetMarkerSize(.5);
  drawHist_p->Draw(drawOpt);
  drawHist_p->Draw("E1 SAME");
}


void makeImbAsymmPTStack(const char* fileName, const char* gorr, const char* PFCaloT, const char* perpProj, const char* GLN, const char* Corr = "")
{
  TFile* panelFile_p = new TFile(fileName, "UPDATE");

  TGraphErrors* getGraphP0_1_p;
  TGraphErrors* getGraphP1_2_p;
  TGraphErrors* getGraphP2_4_p;
  TGraphErrors* getGraphP4_8_p;
  TGraphErrors* getGraphP8_100_p;
  TGraphErrors* getGraphPF_p;

  TGraphErrors* getGraphC0_1_p;
  TGraphErrors* getGraphC1_2_p;
  TGraphErrors* getGraphC2_4_p;
  TGraphErrors* getGraphC4_8_p;
  TGraphErrors* getGraphC8_100_p;
  TGraphErrors* getGraphCF_p;

  Float_t binArrayX[5] = {.00, .10, .20, .30, .50};

  TH1F* histP0_1_p = new TH1F("histP0_1_p", "histP0_1_p", 4, binArrayX);
  TH1F* histP1_2_p = new TH1F("histP1_2_p", "histP1_2_p", 4, binArrayX);
  TH1F* histP2_4_p = new TH1F("histP2_4_p", "histP2_4_p", 4, binArrayX);
  TH1F* histP4_8_p = new TH1F("histP4_8_p", "histP4_8_p", 4, binArrayX);
  TH1F* histP8_100_p = new TH1F("histP8_100_p", "histP8_100_p", 4, binArrayX);
  TH1F* histPF_p = new TH1F("histPF_p", "histPF_p", 4, binArrayX);

  TH1F* histC0_1_p = new TH1F("histC0_1_p", "histC0_1_p", 4, binArrayX);
  TH1F* histC1_2_p = new TH1F("histC1_2_p", "histC1_2_p", 4, binArrayX);
  TH1F* histC2_4_p = new TH1F("histC2_4_p", "histC2_4_p", 4, binArrayX);
  TH1F* histC4_8_p = new TH1F("histC4_8_p", "histC4_8_p", 4, binArrayX);
  TH1F* histC8_100_p = new TH1F("histC8_100_p", "histC8_100_p", 4, binArrayX);
  TH1F* histCF_p = new TH1F("histCF_p", "histCF_p", 4, binArrayX);

  niceTH1(histP0_1_p, 60, -60, 505, 406);
  niceTH1(histP1_2_p, 60, -60, 505, 406);
  niceTH1(histP2_4_p, 60, -60, 505, 406);
  niceTH1(histP4_8_p, 60, -60, 505, 406);
  niceTH1(histP8_100_p, 60, -60, 505, 406);
  niceTH1(histPF_p, 60, -60, 505, 406);

  niceTH1(histC0_1_p, 60, -60, 505, 406);
  niceTH1(histC1_2_p, 60, -60, 505, 406);
  niceTH1(histC2_4_p, 60, -60, 505, 406);
  niceTH1(histC4_8_p, 60, -60, 505, 406);
  niceTH1(histC8_100_p, 60, -60, 505, 406);
  niceTH1(histCF_p, 60, -60, 505, 406);

  getGraphP0_1_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, "0_1", Corr, 30, 100, GLN, fileTag));
  getGraphP1_2_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, "1_2", Corr, 30, 100, GLN, fileTag));
  getGraphP2_4_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, "2_4", Corr, 30, 100, GLN, fileTag));
  getGraphP4_8_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, "4_8", Corr, 30, 100, GLN, fileTag));
  getGraphP8_100_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, "8_100", Corr, 30, 100, GLN, fileTag));
  getGraphPF_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, "F", Corr, 30, 100, GLN, fileTag));

  makeHistForPTStack(getGraphP0_1_p, getGraphP1_2_p, getGraphP2_4_p, getGraphP4_8_p, getGraphP8_100_p, getGraphPF_p, histP0_1_p, histP1_2_p, histP2_4_p, histP4_8_p, histP8_100_p, histPF_p, "P");

  TCanvas* profPanel_p = new TCanvas(Form("%s%sImbAsymm%s%sPTStack_%s_%s_c", gorr, PFCaloT, perpProj, Corr, GLN, fileTag), Form("%s%sImbAsymm%s%sPTStack_%s_%s_c", gorr, PFCaloT, perpProj, Corr, GLN, fileTag), 1);
  profPanel_p->Divide(2, 1, 0, 0);

  profPanel_p->cd(1);

  drawHistToPTStack(histP0_1_p, kBlue - 9, "E1 HIST");
  drawHistToPTStack(histP1_2_p, kYellow - 9, "E1 HIST SAME");
  drawHistToPTStack(histP2_4_p, kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(histP4_8_p, kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(histP8_100_p, kRed + 1, "E1 HIST SAME");
  histPF_p->Draw("SAME E1");

  TLine* zeroLine_p = new TLine(0., 0., 0.5, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw();

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->DrawLatex(.2, .3, "30-100%");

  profPanel_p->cd(2);

  getGraphC0_1_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, "0_1", Corr, 0, 30, GLN, fileTag));
  getGraphC1_2_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, "1_2", Corr, 0, 30, GLN, fileTag));
  getGraphC2_4_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, "2_4", Corr, 0, 30, GLN, fileTag));
  getGraphC4_8_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, "4_8", Corr, 0, 30, GLN, fileTag));
  getGraphC8_100_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, "8_100", Corr, 0, 30, GLN, fileTag));
  getGraphCF_p = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%sF%s_%d%d_%s_%s_g", gorr, PFCaloT, perpProj, Corr, 0, 30, GLN, fileTag));

  makeHistForPTStack(getGraphC0_1_p, getGraphC1_2_p, getGraphC2_4_p, getGraphC4_8_p, getGraphC8_100_p, getGraphCF_p, histC0_1_p, histC1_2_p, histC2_4_p, histC4_8_p, histC8_100_p, histCF_p, "C");

  drawHistToPTStack(histC0_1_p, kBlue - 9, "E1 HIST");
  drawHistToPTStack(histC1_2_p, kYellow - 9, "E1 HIST SAME");
  drawHistToPTStack(histC2_4_p, kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(histC4_8_p, kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(histC8_100_p, kRed + 1, "E1 HIST SAME");

  histCF_p->Draw("SAME E1");

  zeroLine_p->Draw();
  label_p->DrawLatex(.2, .3, "0-30%");

  label_p->DrawLatex(.1, .92, "Lead Jet p_{T} > 120 GeV/c");
  label_p->DrawLatex(.1, .88, "Sublead Jet p_{T} > 50 GeV/c");
  label_p->DrawLatex(.1, .84, "Jet #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.1, .80, "Lead/Sublead Jet abs(#eta) < 1.6");

  profPanel_p->Write();

  delete label_p;
  delete zeroLine_p;
  delete profPanel_p;

  delete histCF_p;
  delete histC8_100_p;
  delete histC4_8_p;
  delete histC2_4_p;
  delete histC1_2_p;
  delete histC0_1_p;

  delete histPF_p;
  delete histP8_100_p;
  delete histP4_8_p;
  delete histP2_4_p;
  delete histP1_2_p;
  delete histP0_1_p;

  panelFile_p->Close();
  delete panelFile_p;
}


void cfmDiJetAsymm(const char* inName = "inFile_CFMHIST_GAMMA.root", bool montecarlo = 0, const char* outName = "outFile_CFMHIST_GAMMA.root")
{
  if(!strcmp(inName, Di80a)){
    std::cout << Di80a << std::endl;
    fileTag = "Di80a";
  }
  else if(!strcmp(inName, Di80b)){
    std::cout << Di80b << std::endl;
    fileTag = "Di80b";
  }
  else if(!strcmp(inName, Di100a)){
    std::cout << Di100a << std::endl;
    fileTag = "Di100a";
  }
  else if(!strcmp(inName, EmDi80a)){
    std::cout << EmDi80a << std::endl;
    fileTag = "EmDi80a";
  }
  else if(!strcmp(inName, Di80c)){
    std::cout << Di80c << std::endl;
    fileTag = "Di80c";
  }
  else if(!strcmp(inName, Di80d)){
    std::cout << Di80d << std::endl;
    fileTag = "Di80d";
  }

  std::cout << "Filetag is: " << fileTag << std::endl;

  inFile_p = new TFile(inName, "READ");
  inTree_p = (TTree*)inFile_p->Get("jetTree");
  inTree_p->AddFriend("trackTree");

  if(montecarlo)
    inTree_p->AddFriend("genTree");
  /*  
  //Asymm Hists, Full Reco
  makeAsymmHist(inTree_p, outName, "PF", 10, 0, 1, 0, 39);
  makeAsymmHist(inTree_p, outName, "PF", 10, 0, 1, 0, 3);
  makeAsymmHist(inTree_p, outName, "PF", 10, 0, 1, 4, 7);
  makeAsymmHist(inTree_p, outName, "PF", 10, 0, 1, 8, 11);
  makeAsymmHist(inTree_p, outName, "PF", 10, 0, 1, 12, 19);
  makeAsymmHist(inTree_p, outName, "PF", 10, 0, 1, 20, 39);

  makeAsymmPanel(outName, "PF");

  makeAsymmHist(inTree_p, outName, "Calo", 10, 0, 1, 0, 39);
  makeAsymmHist(inTree_p, outName, "Calo", 10, 0, 1, 0, 3);
  makeAsymmHist(inTree_p, outName, "Calo", 10, 0, 1, 4, 7);
  makeAsymmHist(inTree_p, outName, "Calo", 10, 0, 1, 8, 11);
  makeAsymmHist(inTree_p, outName, "Calo", 10, 0, 1, 12, 19);
  makeAsymmHist(inTree_p, outName, "Calo", 10, 0, 1, 20, 39);

  makeAsymmPanel(outName, "Calo");
  
  //Pt Proj Hists, Full Reco
  makePtProjHist(inTree_p, outName, "trk", "PF", 0., 1., 0, 11, "Symm");
  makePtProjHist(inTree_p, outName, "trk", "PF", 1., 8., 0, 11, "Symm");
  makePtProjHist(inTree_p, outName, "trk", "PF", 0., 1., 12, 39, "Symm");
  makePtProjHist(inTree_p, outName, "trk", "PF", 0., 1., 0, 11, "Asymm");

  makePtProjHist(inTree_p, outName, "trk", "Calo", 0., 1., 0, 11, "Symm");
  makePtProjHist(inTree_p, outName, "trk", "Calo", 1., 8., 0, 11, "Symm");
  makePtProjHist(inTree_p, outName, "trk", "Calo", 0., 1., 12, 39, "Symm");
  makePtProjHist(inTree_p, outName, "trk", "Calo", 0., 1., 0, 11, "Asymm");
  */
  //Imbalance v. Asymm, Graph, Full Reco

  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 0, 39, -40, 40, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 0, 11, -40, 40, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 12, 39, -40, 40, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 0, 39, -40, 40, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 0, 11, -40, 40, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 12, 39, -40, 40, "N", "Corr");
  /*  
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 0, 39, -40, 40, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 0, 11, -40, 40, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 12, 39, -40, 40, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 0, 39, -40, 40, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 0, 11, -40, 40, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 12, 39, -40, 40, "N", "Corr");
  
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 0, 11, -40, 40, "G");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 12, 39, -40, 40, "G");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 0, 11, -40, 40, "G", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 12, 39, -40, 40, "G", "Corr");

  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 0, 11, -40, 40, "G");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 12, 39, -40, 40, "G");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 0, 11, -40, 40, "G", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 12, 39, -40, 40, "G", "Corr");

  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 0, 11, -40, 40, "L");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 12, 39, -40, 40, "L");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 0, 11, -40, 40, "L", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "F", 12, 39, -40, 40, "L", "Corr");

  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 0, 11, -40, 40, "L");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 12, 39, -40, 40, "L");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 0, 11, -40, 40, "L", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Perp", "F", 12, 39, -40, 40, "L", "Corr");
  */  

  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "0_1", 0, 11, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "1_2", 0, 11, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "2_4", 0, 11, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "4_8", 0, 11, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "8_100", 0, 11, -60, 60, "N", "");

  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "0_1", 12, 39, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "1_2", 12, 39, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "2_4", 12, 39, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "4_8", 12, 39, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "8_100", 12, 39, -60, 60, "N", "");

  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "0_1", 0, 11, -60, 60, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "1_2", 0, 11, -60, 60, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "2_4", 0, 11, -60, 60, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "4_8", 0, 11, -60, 60, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "8_100", 0, 11, -60, 60, "N", "Corr");

  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "0_1", 12, 39, -60, 60, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "1_2", 12, 39, -60, 60, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "2_4", 12, 39, -60, 60, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "4_8", 12, 39, -60, 60, "N", "Corr");
  makeImbAsymmGraph(inTree_p, outName, "r", "PF", "Proj", "8_100", 12, 39, -60, 60, "N", "Corr");

  makeImbAsymmPTStack(outName, "r", "PF", "Proj", "N", "");
  makeImbAsymmPTStack(outName, "r", "PF", "Proj", "N", "Corr");

  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "F", 0, 39, -40, 40, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "F", 0, 11, -40, 40, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "F", 12, 39, -40, 40, "N", "");

  /*  
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Perp", "F", 0, 39, -40, 40, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Perp", "F", 0, 11, -40, 40, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Perp", "F", 12, 39, -40, 40, "N", "");
  
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "F", 0, 11, -40, 40, "G");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "F", 12, 39, -40, 40, "G");

  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Perp", "F", 0, 11, -40, 40, "G");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Perp", "F", 12, 39, -40, 40, "G");

  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "F", 0, 11, -40, 40, "L");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "F", 12, 39, -40, 40, "L");

  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Perp", "F", 0, 11, -40, 40, "L");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Perp", "F", 12, 39, -40, 40, "L");
  */  

  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "0_1", 0, 11, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "1_2", 0, 11, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "2_4", 0, 11, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "4_8", 0, 11, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "8_100", 0, 11, -60, 60, "N", "");

  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "0_1", 12, 39, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "1_2", 12, 39, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "2_4", 12, 39, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "4_8", 12, 39, -60, 60, "N", "");
  makeImbAsymmGraph(inTree_p, outName, "r", "Calo", "Proj", "8_100", 12, 39, -60, 60, "N", "");

  makeImbAsymmPTStack(outName, "r", "Calo", "Proj", "N", "");

  if(montecarlo){
    //Asymm Hists, Truth
    makeAsymmHist(inTree_p, outName, "T", 10, 0, 1, 0, 39);
    makeAsymmHist(inTree_p, outName, "T", 10, 0, 1, 0, 3);
    makeAsymmHist(inTree_p, outName, "T", 10, 0, 1, 4, 7);
    makeAsymmHist(inTree_p, outName, "T", 10, 0, 1, 8, 11);
    makeAsymmHist(inTree_p, outName, "T", 10, 0, 1, 12, 19);
    makeAsymmHist(inTree_p, outName, "T", 10, 0, 1, 20, 39);

    makeAsymmPanel(outName, "T");

    //Pt Proj Hists, Truth (either jets or particles) 
    makePtProjHist(inTree_p, outName, "trk", "T", 0., 1., 0, 11, "Symm");
    makePtProjHist(inTree_p, outName, "trk", "T", 1., 8., 0, 11, "Symm");
    makePtProjHist(inTree_p, outName, "trk", "T", 0., 1., 12, 39, "Symm");
    makePtProjHist(inTree_p, outName, "trk", "T", 0., 1., 0, 11, "Asymm");
    
    makePtProjHist(inTree_p, outName, "gen", "PF", 0., 1., 0, 11, "Symm");
    makePtProjHist(inTree_p, outName, "gen", "PF", 1., 8., 0, 11, "Symm");
    makePtProjHist(inTree_p, outName, "gen", "PF", 0., 1., 12, 39, "Symm");
    makePtProjHist(inTree_p, outName, "gen", "PF", 0., 1., 0, 11, "Asymm");
    
    makePtProjHist(inTree_p, outName, "gen", "Calo", 0., 1., 0, 11, "Symm");
    makePtProjHist(inTree_p, outName, "gen", "Calo", 1., 8., 0, 11, "Symm");
    makePtProjHist(inTree_p, outName, "gen", "Calo", 0., 1., 12, 39, "Symm");
    makePtProjHist(inTree_p, outName, "gen", "Calo", 0., 1., 0, 11, "Asymm");

    makePtProjHist(inTree_p, outName, "gen", "T", 0., 1., 0, 11, "Symm");
    makePtProjHist(inTree_p, outName, "gen", "T", 1., 8., 0, 11, "Symm");
    makePtProjHist(inTree_p, outName, "gen", "T", 0., 1., 12, 39, "Symm");
    makePtProjHist(inTree_p, outName, "gen", "T", 0., 1., 0, 11, "Asymm");
  }
  

  inFile_p->Close();
  delete inFile_p;
}
