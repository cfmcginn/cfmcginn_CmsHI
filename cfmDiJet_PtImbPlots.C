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

TFile* inFile1_p = 0;
TTree* inTree_p = 0;

TFile* outFile_p = 0;

const char* algType[5] = {"PuPF", "PuCalo", "VsPF", "VsCalo", "T"};

TCut fullVeto = "";

//append to every histo so know which sample via shorthand on workblog                                                          

const char* fileTag1;

//shorthands, w/ _CFMSKIM.h                                                                                                                      

const char* Di80a = "Pythia80_HydjetDrum_mix01_HiForest2_v20_CFMSKIM.root";

const char* Di80b = "Dijet80_HydjetDrum_v27_mergedV1_CFMSKIM.root";

const char* Di80c = "HydjetDrum_Pyquen_Dijet80_Embedded_d20140122_Track7_v2_CFMSKIM.root";

const char* Di80d = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";                            
const char* Di80e = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet14_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";

const char* Di80f = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_50k_CFMSKIM.root";

const char* Di80g = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_300k_CFMSKIM.root";

const char* Di80h = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0_300k_CFMSKIM.root";

const char* Di80i = "HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt_0_CFMSKIM.root";

const char* Di100a = "Dijet100_HydjetDrum_v27_mergedV1_CFMSKIM.root";

const char* Di100b = "HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHat_v0_0_CFMSKIM_20140323_4_0.root";

const char* Di120a = "HydjetDrum_Pyquen_Dijet120_Embedded_RECO_STARTHI53_LV1_Track8_Jet21_300k_v0_merged_0_CFMSKIM.root";

const char* Di120b = "HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHat_v0_0_CFMSKIM.root";

const char* EmDi80a = "PbPb_pythiaHYDJET_forest_EmEnrichedDijet80_CFMSKIM.root";

const char* DataA = "Track8_Jet17_GR_R_53_LV6_SUB_0_CFMSKIM.root";

const char* DataB = "hiForest_Jet80or95_GR_R_53_LV6_02Mar2014_1300CET_Track8_Jet15_0_1200k_CFMSKIM.root";

const char* DataC = "hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21_0_700k_CFMSKIM.root";

const char* DataD = "hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21_0_1200k_CFMSKIM.root";

const char* DataE = "hiForest_Jet80or95_GR_R_53_LV6_03Mar2014_1600CET_CMSSW_5_3_16_merged_0_CFMSKIM.root";

//const char* DataE = "hiForest_Jet80or95_GR_R_53_LV6_03Mar2014_1600CET_CMSSW_5_3_16_merged_0_CFMSKIM_20140328_JtCutDown_2_0.root";

Double_t quadSum(Double_t one, Double_t two)
{
  Double_t err = TMath::Sqrt(one*one + two*two);
  return err;
}


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
  if((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0)) return true;

  return false;
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

// 0 == PuPF, 1 == PuCalo, 2 == VsPF, 3 == VsCalo, 4 == Truth

TCut makeSetCut(Int_t setNum)
{
  TCut setCut = "";

  if(setNum > 4 || setNum < 0){
    std::cout << "makeSetCut: setNum must be between 0-4, empty cut returned" << std::endl;
    return setCut;
  }

  setCut = Form("eventSet[%d]", setNum);

  return setCut;
}


TCut makeCentCut(Int_t centLow, Int_t centHi)
{
  TCut centCut = "";
  if(centLow >= 0 && centHi >= centLow && centHi <= 199)
    centCut = Form("hiBin >= %d && hiBin <= %d", centLow, centHi);
  else
    std::cout << "makeCentCut: centLow/centHi incorrectly specified, empty cut returned" << std::endl;

  return centCut;
}


TCut makeAsymmCut(Int_t setNum, Float_t asymmLow, Float_t asymmHi, Bool_t ref = false)
{
  TCut asymmCut = "";

  if(setNum > 4 || setNum < 0){
    std::cout << "makeAsymmCut: setNum must be between 0-4, empty cut returned" << std::endl;
    return asymmCut;
  }

  if(asymmLow >= .00 && asymmHi >= asymmLow && asymmHi <= 1.){
    if(ref && setNum < 4){
      const char* refAsymm = Form("(AlgLeadRefPt[%d] - AlgSubLeadRefPt[%d])/(AlgLeadRefPt[%d] + AlgSubLeadRefPt[%d])", setNum, setNum, setNum, setNum);
      asymmCut = Form("%s > %f && %s < %f && AlgLeadRefPt[%d] > 0 && AlgSubLeadRefPt[%d] > 0", refAsymm, asymmLow, refAsymm, asymmHi, setNum, setNum);
    }
    else
      asymmCut = Form("AlgJtAsymm[%d] > %f && AlgJtAsymm[%d] < %f ", setNum, asymmLow, setNum, asymmHi);
  }
  else
    std::cout << "makeAsymmCut: asymmLow/asymmHi incorrectly specified, empty cut returned" << std::endl;

  return asymmCut;
}


TCut makeEtaCut(Int_t setNum, Float_t overallCut = 2.0, const char* GLN = "N")
{
  TCut etaCut = "";

  if(setNum > 4 || setNum < 0){
    std::cout << "makeEtaCut: setNum must be between 0-4, empty cut returned" << std::endl;
    return etaCut;
  }

  const char* leadJt = Form("AlgLeadJtEta[%d]", setNum);
  const char* subLeadJt = Form("AlgSubLeadJtEta[%d]", setNum);

  etaCut = Form("TMath::Abs(%s) < %f && TMath::Abs(%s) < %f", leadJt, overallCut, subLeadJt, overallCut);

  if(strcmp(GLN, "N") == 0)
    return etaCut;

  if(strcmp(GLN, "G") == 0){
    etaCut = Form("(TMath::Abs(%s) > 1.0 || TMath::Abs(%s) > 1.0) && TMath::Abs(%s) < %f && TMath::Abs(%s) < %f", leadJt, subLeadJt, leadJt, overallCut, subLeadJt, overallCut);
  }
  else if(strcmp(GLN, "L") == 0){
    etaCut = Form("TMath::Abs(%s) < 1.0 && TMath::Abs(%s) < 1.0", leadJt, subLeadJt);
  }
  else
    std::cout << "makeEtaCut: GLN incorrectly specified, returning blank cut" << std::endl;

  return etaCut;
}



TCut makeDelPhiCut(Int_t setNum, Float_t delPhiLow = 2*TMath::Pi()/3)
{
  TCut delPhiCut = "";

  if(setNum > 4 || setNum < 0){
    std::cout << "makeDelPhiCut: setNum must be between 0-4, empty cut returned" << std::endl;
    return delPhiCut;
  }

  const char* jtDelPhi = Form("AlgJtDelPhi[%d]", setNum);

  delPhiCut = Form("%s > %f", jtDelPhi, delPhiLow);

  return delPhiCut;
}


void makeImbAsymmGraph(TTree* getTree_p, const char* outName, const char* gorr, Int_t setNum, const char* perpProj, const char* CNC, const char* FPT, Int_t centLow, Int_t centHi, Int_t graphLow, Int_t graphHi, const char* GLN = "N", const char* Corr = "", Bool_t montecarlo = false)
{
  inFile1_p->cd();

  Int_t setCorrNum = setNum;
  if(strcmp("", Corr) != 0)
    setCorrNum = setNum + 5;

  const char* title = Form("%s%sImbAsymm%s%s%s%s_%d%d_%s_%s_g", gorr, algType[setNum], perpProj, CNC, FPT, Corr, (Int_t)(centLow*.5), (Int_t)((centHi + 1)*.5), GLN, fileTag1);

  TGraphErrors* imbAsymmGraph_p = new TGraphErrors(4);
  imbAsymmGraph_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTGraphErrors(imbAsymmGraph_p, graphHi, graphLow);

  TH1F* getHist_p;

  TString var = Form("%sAlgImb%s%s%s[%d]", gorr, perpProj, CNC, FPT, setCorrNum);

  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 1.6, GLN);

  TCut phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);
  if(strcmp("", CNC) != 0){
    etaCut = makeEtaCut(setNum, 1.6, GLN);
    phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);
  }

  TCut jetLCut = Form("AlgLeadJtPt[%d] > 120*.95", setNum);

  //  TCut isQuark = Form("isQuarkJet[%d]", setNum);
  //  TCut isGluon = Form("isGluonJet[%d]", setNum);       

  const char* name1[4] = {"0_1(10000, -10000, 10000)", "1_2(10000, -10000, 10000)", "2_3(10000, -10000, 10000)", "3_5(10000, -10000, 10000)"};
  const char* name2[4] = {"0_1", "1_2", "2_3", "3_5"};
  Float_t asymmBins[5] = {.00, .11, .22, .33, 1.};
  Float_t point[4] = {.055, .165, .275, .415};
  Float_t xErr[4] = {.055, .055, .055, .085};

  for(Int_t binIter = 0; binIter < 4; binIter++){
    TCut asymmCut = makeAsymmCut(setNum, asymmBins[binIter], asymmBins[binIter + 1]);

    if(montecarlo)
      getTree_p->Project(name1[binIter], var, Form("setWeight_2pi3[%d]", setNum)*(setCut && centCut && etaCut && phiCut && jetLCut && asymmCut && fullVeto));
    else
      getTree_p->Project(name1[binIter], var, setCut && centCut && etaCut && phiCut && jetLCut && asymmCut);

    getHist_p = (TH1F*)inFile1_p->Get(name2[binIter]);

    imbAsymmGraph_p->SetPoint(binIter, point[binIter], getHist_p->GetMean());
    imbAsymmGraph_p->SetPointError(binIter, xErr[binIter], getHist_p->GetMeanError());
  }

  imbAsymmGraph_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTGraphErrors(imbAsymmGraph_p, graphHi, graphLow);

  outFile_p = new TFile(outName, "UPDATE");
  imbAsymmGraph_p->Write(title);
  outFile_p->Close();

  delete outFile_p;
  delete imbAsymmGraph_p;
}


void makeImbAsymmGraph_Tight(TTree* getTree_p, const char* outName, const char* gorr, Int_t setNum, const char* perpProj, const char* CNC, const char* FPT, Int_t centLow, Int_t centHi, Int_t graphLow, Int_t graphHi, const char* GLN = "N", const char* Corr = "", Bool_t montecarlo = false)
{
  inFile1_p->cd();

  Int_t setCorrNum = setNum;
  if(strcmp("", Corr) != 0)
    setCorrNum = setNum + 5;

  const char* title = Form("%s%sImbAsymmTight%s%s%s%s_%d%d_%s_%s_g", gorr, algType[setNum], perpProj, CNC, FPT, Corr, (Int_t)(centLow*.5), (Int_t)((centHi + 1)*.5), GLN, fileTag1);


  TGraphErrors* imbAsymmGraph_p = new TGraphErrors(4);
  imbAsymmGraph_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTGraphErrors(imbAsymmGraph_p, graphHi, graphLow);

  TH1F* getHist_p;


  TString var = Form("%sAlgImb%s%s%s[%d]", gorr, perpProj, CNC, FPT, setCorrNum);

  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 1.6, GLN);

  TCut phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);
  if(strcmp(CNC, "") != 0){
    etaCut = makeEtaCut(setNum, 1.6, GLN);
    phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);
  }

  TCut jetLCut = Form("AlgLeadJtPt[%d] > 120*.95", setNum);

  //  TCut isQuark = Form("isQuarkJet[%d]", setNum);
  //  TCut isGluon = Form("isGluonJet[%d]", setNum);       

  const char* name1[8] = {"0_h(10000, -10000, 10000)", "1_h(10000, -10000, 10000)", "2_h(10000, -10000, 10000)", "3_h(10000, -10000, 10000)", "4_h(10000, -10000, 10000)", "5_h(10000, -10000, 10000)", "6_h(10000, -10000, 10000)", "7_h(10000, -10000, 10000)"};
  const char* name2[8] = {"0_h", "1_h", "2_h", "3_h", "4_h", "5_h", "6_h", "7_h"};
  Float_t asymmBins[9] = {.00, .055, .11, .165, .22, .275, .33, .415, 1.};
  Float_t point[8] = {.0275, .0825, .1375, .1925, .2475, .3025, .3725, .4575};
  Float_t xErr[8] = {.0275, .0275, .0275, .0275, .0275, .0275, .0425, .0425};

  for(Int_t binIter = 0; binIter < 8; binIter++){
    TCut asymmCut = makeAsymmCut(setNum, asymmBins[binIter], asymmBins[binIter + 1]);

    //    std::cout << setCut << ", " << centCut << ", " << etaCut << ", " << phiCut << ", " << jetLCut << ", " << asymmCut << ", " << fullVeto << std::endl;
    //    std::cout << std::endl;

    if(montecarlo)
      getTree_p->Project(name1[binIter], var, Form("setWeight_2pi3[%d]", setNum)*(setCut && centCut && etaCut && phiCut && jetLCut && asymmCut && fullVeto));
    else
      getTree_p->Project(name1[binIter], var, setCut && centCut && etaCut && phiCut && jetLCut && asymmCut);

    getHist_p = (TH1F*)inFile1_p->Get(name2[binIter]);

    imbAsymmGraph_p->SetPoint(binIter, point[binIter], getHist_p->GetMean());
    imbAsymmGraph_p->SetPointError(binIter, xErr[binIter], getHist_p->GetMeanError());
 }

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


void makeHistForPtStack(TGraph* g0_1_p, TGraph* g1_2_p, TGraph* g2_4_p, TGraph* g4_8_p, TGraph* g8_100_p, TGraph* gF_p, TH1F* h0_1_p, TH1F* h1_2_p, TH1F* h2_4_p, TH1F* h4_8_p, TH1F* h8_100_p, TH1F* hF_p, Int_t pos = 4, const char* Tight = "")
{
  Int_t points = 4;
  if(strcmp(Tight, "Tight") == 0)
    points = 8;

  std::cout << points << std::endl;

  Double_t x0_1[points];
  Double_t y0_1[points];
  Double_t x1_2[points];
  Double_t y1_2[points];
  Double_t x2_4[points];
  Double_t y2_4[points];
  Double_t x4_8[points];
  Double_t y4_8[points];
  Double_t x8_100[points];
  Double_t y8_100[points];
  Double_t xF[points];
  Double_t yF[points];

  for(Int_t iter = 0; iter < points; iter++){
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

    if(iter == 7){
      std::cout << std::endl;

      std::cout << h8_100_p->GetBinContent(iter+1) << std::endl;
      std::cout << h4_8_p->GetBinContent(iter+1) << std::endl;
      std::cout << h2_4_p->GetBinContent(iter+1) << std::endl;
      std::cout << h1_2_p->GetBinContent(iter+1) << std::endl;
      std::cout << h0_1_p->GetBinContent(iter+1) << std::endl;
      std::cout << hF_p->GetBinContent(iter+1) << std::endl;

      std::cout << std::endl;
    }
  }

  h8_100_p->SetXTitle("A_{J}");
  h4_8_p->SetXTitle("A_{J}");
  h2_4_p->SetXTitle("A_{J}");
  h1_2_p->SetXTitle("A_{J}");
  h0_1_p->SetXTitle("A_{J}");
  hF_p->SetXTitle("A_{J}");

  if(pos == 1 || pos == 7){
    h8_100_p->SetYTitle("<#slash{p}_{T}^{||}>, (GeV/c)");
    h4_8_p->SetYTitle("<#slash{p}_{T}^{||}>, (GeV/c)");
    h2_4_p->SetYTitle("<#slash{p}_{T}^{||}>, (GeV/c)");
    h1_2_p->SetYTitle("<#slash{p}_{T}^{||}>, (GeV/c)");
    h0_1_p->SetYTitle("<#slash{p}_{T}^{||}>, (GeV/c)");
    hF_p->SetYTitle("<#slash{p}_{T}^{||}>, (GeV/c)");
  }
  return;
}



void makeSubHistForPtStack(TGraph* g0_1_p, TGraph* g1_2_p, TGraph* g2_4_p, TGraph* g4_8_p, TGraph* g8_100_p, TGraph* gF_p, TH1F* h0_1_p, TH1F* h1_2_p, TH1F* h2_4_p, TH1F* h4_8_p, TH1F* h8_100_p, TH1F* hF_p, TGraph* pp0_1_p, TGraph* pp1_2_p, TGraph* pp2_4_p, TGraph* pp4_8_p, TGraph* pp8_100_p, TGraph* ppF_p, Int_t pos = 4, Bool_t disp = false)
{
  Int_t points = 4;

  Double_t x0_1[points];
  Double_t y0_1[points];
  Double_t x1_2[points];
  Double_t y1_2[points];
  Double_t x2_4[points];
  Double_t y2_4[points];
  Double_t x4_8[points];
  Double_t y4_8[points];
  Double_t x8_100[points];
  Double_t y8_100[points];
  Double_t xF[points];
  Double_t yF[points];

  Double_t xPP0_1[points];
  Double_t yPP0_1[points];
  Double_t xPP1_2[points];
  Double_t yPP1_2[points];
  Double_t xPP2_4[points];
  Double_t yPP2_4[points];
  Double_t xPP4_8[points];
  Double_t yPP4_8[points];
  Double_t xPP8_100[points];
  Double_t yPP8_100[points];
  Double_t xPPF[points];
  Double_t yPPF[points];

  for(Int_t iter = 0; iter < points; iter++){
    g0_1_p->GetPoint(iter, x0_1[iter], y0_1[iter]);
    g1_2_p->GetPoint(iter, x1_2[iter], y1_2[iter]);
    g2_4_p->GetPoint(iter, x2_4[iter], y2_4[iter]);
    g4_8_p->GetPoint(iter, x4_8[iter], y4_8[iter]);
    g8_100_p->GetPoint(iter, x8_100[iter], y8_100[iter]);
    gF_p->GetPoint(iter, xF[iter], yF[iter]);

    pp0_1_p->GetPoint(iter, xPP0_1[iter], yPP0_1[iter]);
    pp1_2_p->GetPoint(iter, xPP1_2[iter], yPP1_2[iter]);
    pp2_4_p->GetPoint(iter, xPP2_4[iter], yPP2_4[iter]);
    pp4_8_p->GetPoint(iter, xPP4_8[iter], yPP4_8[iter]);
    pp8_100_p->GetPoint(iter, xPP8_100[iter], yPP8_100[iter]);
    ppF_p->GetPoint(iter, xPPF[iter], yPPF[iter]);

    if(iter == 3 && disp){
      std::cout << yF[iter] << ", " << yPPF[iter] << ", " << yF[iter] - yPPF[iter] << std::endl;
    }


    y0_1[iter] += -yPP0_1[iter];
    y1_2[iter] += -yPP1_2[iter];
    y2_4[iter] += -yPP2_4[iter];
    y4_8[iter] += -yPP4_8[iter];
    y8_100[iter] += -yPP8_100[iter];
    yF[iter] += -yPPF[iter];

    if(iter == 3)
      std::cout << yF[iter] << std::endl;


    h8_100_p->SetBinContent(iter + 1, y8_100[iter]);
    h8_100_p->SetBinError(iter + 1, quadSum(g8_100_p->GetErrorY(iter), pp8_100_p->GetErrorY(iter)));
    
    h4_8_p->SetBinContent(iter + 1, sumYForPTStack(y4_8[iter], y8_100[iter]));
    h4_8_p->SetBinError(iter + 1, quadSum(g4_8_p->GetErrorY(iter), pp4_8_p->GetErrorY(iter)));
    
    h2_4_p->SetBinContent(iter + 1, sumYForPTStack(y2_4[iter], y4_8[iter], y8_100[iter]));
    h2_4_p->SetBinError(iter + 1, quadSum(g2_4_p->GetErrorY(iter), pp2_4_p->GetErrorY(iter)));

    h1_2_p->SetBinContent(iter + 1, sumYForPTStack(y1_2[iter], y2_4[iter], y4_8[iter], y8_100[iter]));
    h1_2_p->SetBinError(iter + 1, quadSum(g1_2_p->GetErrorY(iter), pp1_2_p->GetErrorY(iter)));

    h0_1_p->SetBinContent(iter + 1, sumYForPTStack(y0_1[iter], y1_2[iter], y2_4[iter], y4_8[iter], y8_100[iter]));
    h0_1_p->SetBinError(iter + 1, quadSum(g0_1_p->GetErrorY(iter), pp0_1_p->GetErrorY(iter)));

    hF_p->SetBinContent(iter + 1, yF[iter]);
    hF_p->SetBinError(iter + 1, quadSum(gF_p->GetErrorY(iter), ppF_p->GetErrorY(iter)));

  }
  
  h8_100_p->SetXTitle("A_{J}");
  h4_8_p->SetXTitle("A_{J}");
  h2_4_p->SetXTitle("A_{J}");
  h1_2_p->SetXTitle("A_{J}");
  h0_1_p->SetXTitle("A_{J}");
  hF_p->SetXTitle("A_{J}");
  
  if(pos == 1 || pos == 7){
    h8_100_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h4_8_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h2_4_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h1_2_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h0_1_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    hF_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  }
  return;
}



void drawHistToPTStack(TH1F* drawHist_p, Int_t color, const char* drawOpt)
{
  drawHist_p->SetFillColor(color);
  drawHist_p->SetMarkerStyle(6);
  drawHist_p->SetMarkerSize(.5);
  drawHist_p->DrawCopy(drawOpt);
  drawHist_p->DrawCopy("E1 SAME");
}



void makeImbAsymmPtStack(const char* fileName, const char* gorr, Int_t setNum, const char* perpProj, const char* GLN, const char* Corr = "", const char* Tight = "", const char* CNC = "", const char* ppName = "", const char* ppFileTag = "", Bool_t montecarlo = false)
{
  TFile* panelFile_p = new TFile(fileName, "UPDATE");
  TFile* ppFile_p;

  Int_t pos[4] = {1, 2, 3, 4};
  
  const char* mcLabelI[4] = {"PYTHIA", "PYTHIA + HYDJET", "PYTHIA + HYDJET", "(PYTHIA + HYDJET) - PYTHIA"};
  const char* dataLabelI[4] = {"pp, #sqrt{s_{NN}} = 2.76 TeV, 5.3/pb", "PbPb, #sqrt{s_{NN}} = 2.76 TeV, 150/#mub", "PbPb", "PbPb - pp"};

  const char* overLabel[4];
  Float_t overCoord[2];

  for(Int_t iter = 0; iter < 4; iter++){
    if(montecarlo)
      overLabel[iter] = mcLabelI[iter];
    else
      overLabel[iter] = dataLabelI[iter];
  }

  if(strcmp(ppFileTag, "") != 0){
    overCoord[0] = .90;
    overCoord[1] = .84;
  }
  else{
    overCoord[0] = .18;
    overCoord[1] = .13;
  }
  

  if(strcmp("", ppName) != 0){
    for(Int_t iter = 0; iter < 4; iter++){
      pos[iter] = pos[iter] + 1;
    }
  }

  TGraphErrors* getGraph1_p[6];
  TGraphErrors* getGraph2_p[6];
  TGraphErrors* getGraph3_p[6];
  TGraphErrors* getGraph4_p[6];

  TGraphErrors* getGraphpp_p[6];

  TH1F* hist1_p[6];
  TH1F* hist2_p[6];
  TH1F* hist3_p[6];
  TH1F* hist4_p[6];

  TH1F* histpp_p[6];

  TH1F* histAxis_p;

  Float_t binArrayX[5] = {.00, .11, .22, .33, .50};
  Float_t binArrayX_Tight[9] = {.00, .055, .11, .165, .22, .275, .33, .415, .5};

  const char* namePT[6] = {"0_1", "1_2", "2_4", "4_8", "8_100", "F"};

  histAxis_p = new TH1F("histAxis", "histAxis", 4, binArrayX);
  niceTH1(histAxis_p, 60, -60, 505, 406);

  histAxis_p->SetBinContent(1, 60);
  histAxis_p->SetBinContent(2, 60);
  histAxis_p->SetBinContent(3, 60);
  histAxis_p->SetBinContent(4, 60);

  histAxis_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");

  for(Int_t histIter = 0; histIter < 6; histIter++){
    if(strcmp(Tight, "") == 0){
      hist1_p[histIter] = new TH1F(Form("hist1%s_p", namePT[histIter]), Form("hist1%s_p", namePT[histIter]), 4, binArrayX);
      hist2_p[histIter] = new TH1F(Form("hist2%s_p", namePT[histIter]), Form("hist2%s_p", namePT[histIter]), 4, binArrayX);

      hist3_p[histIter] = new TH1F(Form("hist3%s_p", namePT[histIter]), Form("hist3%s_p", namePT[histIter]), 4, binArrayX);
      hist4_p[histIter] = new TH1F(Form("hist4%s_p", namePT[histIter]), Form("hist4%s_p", namePT[histIter]), 4, binArrayX);

      histpp_p[histIter] = new TH1F(Form("histpp%s_p", namePT[histIter]), Form("histpp%s_p", namePT[histIter]), 4, binArrayX);
    }
    else{
      hist1_p[histIter] = new TH1F(Form("hist1%s_p", namePT[histIter]), Form("hist1%s_p", namePT[histIter]), 8, binArrayX_Tight);
      hist2_p[histIter] = new TH1F(Form("hist2%s_p", namePT[histIter]), Form("hist2%s_p", namePT[histIter]), 8, binArrayX_Tight);

      hist3_p[histIter] = new TH1F(Form("hist3%s_p", namePT[histIter]), Form("hist3%s_p", namePT[histIter]), 8, binArrayX_Tight);
      hist4_p[histIter] = new TH1F(Form("hist4%s_p", namePT[histIter]), Form("hist4%s_p", namePT[histIter]), 8, binArrayX_Tight);

      histpp_p[histIter] = new TH1F(Form("histpp%s_p", namePT[histIter]), Form("histpp%s_p", namePT[histIter]), 8, binArrayX_Tight);
    }

    niceTH1(hist1_p[histIter], 60, -60, 505, 406);
    niceTH1(hist2_p[histIter], 60, -60, 505, 406);
    niceTH1(hist3_p[histIter], 60, -60, 505, 406);
    niceTH1(hist4_p[histIter], 60, -60, 505, 406);
    niceTH1(histpp_p[histIter], 60, -60, 505, 406);


    if(strcmp(CNC, "") == 0)
      getGraph1_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s%s%s_50100_%s_%s_g", gorr, algType[setNum], Tight, perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
    else
      getGraph1_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s%s%s_30100_%s_%s_g", gorr, algType[setNum], Tight, perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
    
  }

  makeHistForPtStack(getGraph1_p[0], getGraph1_p[1], getGraph1_p[2], getGraph1_p[3], getGraph1_p[4], getGraph1_p[5], hist1_p[0], hist1_p[1], hist1_p[2], hist1_p[3], hist1_p[4], hist1_p[5], pos[0], Tight);

  hist1_p[5]->SetMarkerStyle(5);
  hist1_p[5]->SetFillColor(0);
    
  TCanvas* profPanel_p;
  if(strcmp(ppFileTag, "") == 0){
    profPanel_p = new TCanvas(Form("%s%sImbAsymm%s%s%s%sPTStack_%s_%s_c", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), Form("%s%sImbAsymm%s%s%s%sPTStack_%s_%s_c", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), 700, 250);
    profPanel_p->Divide(4, 1, 0, 0);
    std::cout << "FourPanel Init" << std::endl;
  }
  else{
    if(strcmp(CNC, "") == 0){
      profPanel_p = new TCanvas(Form("%s%sImbAsymm%s%s%s%sPTStackPP_%s_%s_c", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), Form("%s%sImbAsymm%s%s%s%sPTStackPP_%s_%s_c", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), 1000, 500);
      profPanel_p->Divide(5, 2, 0, 0);
      std::cout << "FivePanel Init" << std::endl;
    }
    else{
      profPanel_p = new TCanvas(Form("%s%sImbAsymm%s%s%s%sPTStackPP_%s_%s_c", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), Form("%s%sImbAsymm%s%s%s%sPTStackPP_%s_%s_c", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), 600, 500);
      profPanel_p->Divide(3, 2, 0, 0);
      std::cout << "FivePanel Init" << std::endl;
    }
  }

  const char* qsquare;

  if(strcmp(fileTag1, "DataE") == 0)
    qsquare = Form(" ");
  else if(strcmp(fileTag1, "Di80i") == 0)
    qsquare = Form("(#hat{p}_{T} > 80)");
  else if(strcmp(fileTag1, "Di120b") == 0)
    qsquare = Form("(#hat{p}_{T} > 120)");


  TLegend* leg;
  if(strcmp(gorr, "g") == 0)
    leg = new TLegend(0.2, 0.60, 0.95, 0.95, Form("Truth #slash{p}_{T}^{||}  v. A_{J}   %s", qsquare));
  else 
    leg = new TLegend(0.2, 0.60, 0.95, 0.95, Form("%sTrk #slash{p}_{T}^{||} v. A_{J}   %s", Corr, qsquare));

  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSizePixels(12);
  leg->SetBorderSize(0);

  profPanel_p->cd(pos[0]);

  drawHistToPTStack(hist1_p[0], kBlue - 9, "E1 HIST");
  leg->AddEntry(hist1_p[0], ".5 < p_{T} < 1", "F");
  drawHistToPTStack(hist1_p[1], kYellow - 9, "E1 HIST SAME");
  leg->AddEntry(hist1_p[1], "1 < p_{T} < 2", "F");
  drawHistToPTStack(hist1_p[2], kOrange + 1, "E1 HIST SAME");
  leg->AddEntry(hist1_p[2], "2 < p_{T} < 4", "F");
  drawHistToPTStack(hist1_p[3], kGreen + 3, "E1 HIST SAME");
  leg->AddEntry(hist1_p[3], "4 < p_{T} < 8", "F");
  drawHistToPTStack(hist1_p[4], kRed + 1, "E1 HIST SAME");
  leg->AddEntry(hist1_p[4], "8 < p_{T}", "F");

  hist1_p[5]->DrawCopy("SAME E1");
  
  
  Int_t sysA[5] = {0, 3, 3, 3, 3};
  Int_t sysC[5] = {0, 4, 6, 8, 6};
  Int_t sysNC[5] = {0, 1, 4, 6, 6};

  for(Int_t iter = 1; iter <= hist1_p[5]->GetNbinsX(); iter++){
    Float_t yVal = hist1_p[5]->GetBinContent(iter);
    TLine* l = 0;

    if(strcmp(CNC, "") == 0){

      std::cout <<"Iter: " << iter << ", " << yVal << ", " << yVal*.04 << std::endl;

      l = new TLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter+1) - .01, yVal - TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal));
      l->SetLineColor(1);
      l->Draw();
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter+1) -.01, yVal + TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal));

      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal) - 2, hist1_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal));
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal) - 2, hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sysA[iter]*sysA[iter] + .04*yVal*.04*yVal));

    }
    else if(strcmp(CNC, "C") == 0){

      l = new TLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter+1) - .01, yVal - TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal));
      l->SetLineColor(1);
      l->Draw();
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter+1) -.01, yVal + TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal));

      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal) - 2, hist1_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal));
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal) - 2, hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sysC[iter]*sysC[iter] + .04*yVal*.04*yVal));

    }
    else{

      l = new TLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter+1) - .01, yVal - TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal));
      l->SetLineColor(1);
      l->Draw();
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal), hist1_p[5]->GetBinLowEdge(iter+1) -.01, yVal + TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal));

      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal) - 2, hist1_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal));
      l->DrawLine(hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal) - 2, hist1_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sysNC[iter]*sysNC[iter] + .04*yVal*.04*yVal));

    }
  }
  
  

  if(strcmp(ppName, "") == 0)
    leg->Draw("SAME");
  else{
    if(strcmp("", CNC) == 0)
      profPanel_p->cd(6);
    else
      profPanel_p->cd(4);

    histAxis_p->DrawCopy("HIST");
    leg->Draw("SAME");
  }

  profPanel_p->cd(pos[0]);

  TLine* zeroLine_p = new TLine(0., 0., 0.5, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw();

  TLatex* label1_p = new TLatex();
  label1_p->SetNDC();
  label1_p->SetTextSizePixels(18);

  TLatex* label2_p = new TLatex();
  label2_p->SetNDC();
  label2_p->SetTextSizePixels(8);

  label1_p->DrawLatex(.22, overCoord[0], Form("%s", overLabel[1]));
  if(strcmp("", CNC) == 0)
    label1_p->DrawLatex(.22, overCoord[1], "50-100%");
  else
    label1_p->DrawLatex(.22, overCoord[1], "30-100%");

  profPanel_p->cd(pos[1]);

  for(Int_t histIter = 0; histIter < 6; histIter++){
    if(strcmp(CNC, "") == 0)
      getGraph2_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s%s%s_3050_%s_%s_g", gorr, algType[setNum], Tight, perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
    else
      getGraph2_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s%s%s_030_%s_%s_g", gorr, algType[setNum], Tight, perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
  }

  makeHistForPtStack(getGraph2_p[0], getGraph2_p[1], getGraph2_p[2], getGraph2_p[3], getGraph2_p[4], getGraph2_p[5], hist2_p[0], hist2_p[1], hist2_p[2], hist2_p[3], hist2_p[4], hist2_p[5], pos[1], Tight);

  hist2_p[5]->SetMarkerStyle(5);
  hist2_p[5]->SetFillColor(0);

  drawHistToPTStack(hist2_p[0], kBlue - 9, "E1 HIST");
  drawHistToPTStack(hist2_p[1], kYellow - 9, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[2], kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[3], kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[4], kRed + 1, "E1 HIST SAME");

  hist2_p[5]->DrawCopy("SAME E1");
  
  Int_t sys2A[5] = {0, 3, 3, 3, 3};
  Int_t sys2C[5] = {0, 1, 5, 4, 5};
  Int_t sys2NC[5] = {0, 1, 7, 6, 4};

  for(Int_t iter = 1; iter <= hist2_p[5]->GetNbinsX(); iter++){
    Float_t yVal = hist2_p[5]->GetBinContent(iter);
    TLine* l = 0;

    if(strcmp(CNC, "") == 0){

      std::cout <<"Iter: " << iter << ", " << yVal << ", " << yVal*.04 << std::endl;

      l = new TLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter+1) - .01, yVal - TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal));
      l->SetLineColor(1);
      l->Draw();
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter+1) -.01, yVal + TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal));

      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal) - 2, hist2_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal));
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal) - 2, hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys2A[iter]*sys2A[iter] + .04*yVal*.04*yVal));

    }
    else if(strcmp(CNC, "C") == 0){

      l = new TLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter+1) - .01, yVal - TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal));
      l->SetLineColor(1);
      l->Draw();
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter+1) -.01, yVal + TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal));

      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal) - 2, hist2_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal));
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal) - 2, hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys2C[iter]*sys2C[iter] + .04*yVal*.04*yVal));

    }
    else{

      l = new TLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter+1) - .01, yVal - TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal));
      l->SetLineColor(1);
      l->Draw();
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal) + 2);
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal), hist2_p[5]->GetBinLowEdge(iter+1) -.01, yVal + TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal));

      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal) - 2, hist2_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal));
      l->DrawLine(hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal) - 2, hist2_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys2NC[iter]*sys2NC[iter] + .04*yVal*.04*yVal));

    }
  }
     
  zeroLine_p->Draw();


  label1_p->DrawLatex(.22, overCoord[0], Form("%s", overLabel[2]));
  if(strcmp("", CNC) == 0)
    label1_p->DrawLatex(.22, overCoord[1], "30-50%");
  else
    label1_p->DrawLatex(.22, overCoord[1], "0-30%");

  if(strcmp(ppName, "") == 0){
    label2_p->DrawLatex(.1, .92, Form("%s Leap_{T,1} > 120 GeV/c", algType[setNum]));
    label2_p->DrawLatex(.1, .86, Form("%s Sublead Jet p_{T} > 50 GeV/c", algType[setNum]));
    if(strcmp(CNC, "") != 0){
      label2_p->DrawLatex(.1, .80, Form("%s Jet |#eta| < 1.6", algType[setNum]));
      label2_p->DrawLatex(.1, .74, Form("%s Jet #Delta #phi > 5#pi/6", algType[setNum]));
    }
    else{
      label2_p->DrawLatex(.1, .80, Form("%s Jet |#eta| < 1.6", algType[setNum]));
      label2_p->DrawLatex(.1, .74, Form("%s Jet #Delta #phi > 5#pi/6", algType[setNum]));
    }
  }
  else{
    if(strcmp(CNC, "") == 0)
      profPanel_p->cd(6);
    else
      profPanel_p->cd(4);

    label2_p->DrawLatex(.24, .52, "Jet R = .3");
    label2_p->DrawLatex(.24, .47, "p_{T,1} > 120, p_{T,2} > 50 GeV/c");
    if(strcmp(CNC, "") == 0){
      label2_p->DrawLatex(.24, .42, "Jet |#eta|_{1}, |#eta|_{2} < 1.6");
      label2_p->DrawLatex(.24, .37, "Jet #Delta #phi > 5#pi/6");
    }
    else{
      label2_p->DrawLatex(.24, .42, "Jet |#eta|_{1}, |#eta|_{2} < 1.6");
      label2_p->DrawLatex(.24, .37, "Jet #Delta #phi > 5#pi/6");
    }

    if(strcmp(CNC, "C") == 0)
      label2_p->DrawLatex(.26, .32, Form("In-Cone, #Delta R < .8"));
    else if(strcmp(CNC, "NC") == 0)
      label2_p->DrawLatex(.26, .32, Form("Out-of-Cone, .8 < #Delta R"));
    else if(strcmp(CNC, "NCCut") == 0)
      label2_p->DrawLatex(.26, .32, Form("Out-of-Cone, .8 < #Delta R < #pi/2"));
  }

  if(strcmp(CNC, "") == 0){
    profPanel_p->cd(pos[2]);

    for(Int_t histIter = 0; histIter < 6; histIter++){
      getGraph3_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s%s%s_1030_%s_%s_g", gorr, algType[setNum], Tight, perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
    }

    makeHistForPtStack(getGraph3_p[0], getGraph3_p[1], getGraph3_p[2], getGraph3_p[3], getGraph3_p[4], getGraph3_p[5], hist3_p[0], hist3_p[1], hist3_p[2], hist3_p[3], hist3_p[4], hist3_p[5], pos[2], Tight);

    hist3_p[5]->SetMarkerStyle(5);
    hist3_p[5]->SetFillColor(0);

    drawHistToPTStack(hist3_p[0], kBlue - 9, "E1 HIST");
    drawHistToPTStack(hist3_p[1], kYellow - 9, "E1 HIST SAME");
    drawHistToPTStack(hist3_p[2], kOrange + 1, "E1 HIST SAME");
    drawHistToPTStack(hist3_p[3], kGreen + 3, "E1 HIST SAME");
    drawHistToPTStack(hist3_p[4], kRed + 1, "E1 HIST SAME");
  
    hist3_p[5]->DrawCopy("SAME E1");
    
    
    Int_t sys3A[5] = {0, 5, 5, 5, 5};
    
    for(Int_t iter = 1; iter <= hist3_p[5]->GetNbinsX(); iter++){
      Float_t yVal = hist3_p[5]->GetBinContent(iter);
      TLine* l = 0;
      
      if(strcmp(CNC, "") == 0){
	std::cout <<"Iter: " << iter << ", " << yVal << ", " << yVal*.04 << std::endl;
	
	l = new TLine(hist3_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal), hist3_p[5]->GetBinLowEdge(iter+1) - .01, yVal - TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal));
	l->SetLineColor(1);
	l->Draw();
	l->DrawLine(hist3_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal), hist3_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal) + 2);
	l->DrawLine(hist3_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal), hist3_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal) + 2);
	l->DrawLine(hist3_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal), hist3_p[5]->GetBinLowEdge(iter+1) -.01, yVal + TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal));
      
	l->DrawLine(hist3_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal) - 2, hist3_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal));
	l->DrawLine(hist3_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal) - 2, hist3_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys3A[iter]*sys3A[iter] + .04*yVal*.04*yVal));
	
      }
    }
        
    zeroLine_p->Draw();
    
    label1_p->DrawLatex(.22, overCoord[0], Form("%s", overLabel[2]));
    label1_p->DrawLatex(.22, overCoord[1], "10-30%");  
    
    profPanel_p->cd(pos[3]);

    for(Int_t histIter = 0; histIter < 6; histIter++){
      getGraph4_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbAsymm%s%s%s%s%s_010_%s_%s_g", gorr, algType[setNum], Tight, perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
    }

    makeHistForPtStack(getGraph4_p[0], getGraph4_p[1], getGraph4_p[2], getGraph4_p[3], getGraph4_p[4], getGraph4_p[5], hist4_p[0], hist4_p[1], hist4_p[2], hist4_p[3], hist4_p[4], hist4_p[5], pos[3], Tight);

    hist4_p[5]->SetMarkerStyle(5);
    hist4_p[5]->SetFillColor(0);

    drawHistToPTStack(hist4_p[0], kBlue - 9, "E1 HIST");
    drawHistToPTStack(hist4_p[1], kYellow - 9, "E1 HIST SAME");
    drawHistToPTStack(hist4_p[2], kOrange + 1, "E1 HIST SAME");
    drawHistToPTStack(hist4_p[3], kGreen + 3, "E1 HIST SAME");
    drawHistToPTStack(hist4_p[4], kRed + 1, "E1 HIST SAME");
  
    hist4_p[5]->DrawCopy("SAME E1");
    
    Int_t sys4A[5] = {0, 5, 5, 5, 5};

    
    for(Int_t iter = 1; iter <= hist4_p[5]->GetNbinsX(); iter++){
      Float_t yVal = hist4_p[5]->GetBinContent(iter);
      TLine* l = 0;
      
      if(strcmp(CNC, "") == 0){

	std::cout <<"Iter: " << iter << ", " << yVal << ", " << yVal*.04 << std::endl;


	l = new TLine(hist4_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal), hist4_p[5]->GetBinLowEdge(iter+1) - .01, yVal - TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal));
	l->SetLineColor(1);
	l->Draw();
	l->DrawLine(hist4_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal), hist4_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal) + 2);
	l->DrawLine(hist4_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal), hist4_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal) + 2);
	l->DrawLine(hist4_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal), hist4_p[5]->GetBinLowEdge(iter+1) -.01, yVal + TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal));
	
	l->DrawLine(hist4_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal) - 2, hist4_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal));
	l->DrawLine(hist4_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal) - 2, hist4_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys4A[iter]*sys4A[iter] + .04*yVal*.04*yVal));
	
      }
    }
           

    zeroLine_p->Draw();
    
    label1_p->DrawLatex(.22, overCoord[0], Form("%s", overLabel[2]));
    label1_p->DrawLatex(.22, overCoord[1], "0-10%");
  }

  if(strcmp(ppFileTag, "") != 0){
    ppFile_p = new TFile(ppName, "READ");

    profPanel_p->cd(1);

    for(Int_t histIter = 0; histIter < 6; histIter++){
      getGraphpp_p[histIter] = (TGraphErrors*)ppFile_p->Get(Form("%s%sImbAsymm%s%s%s%s%s_%s_%s_g", gorr, algType[setNum], Tight, perpProj, CNC, namePT[histIter], Corr, GLN, ppFileTag));
    }
    
    makeHistForPtStack(getGraphpp_p[0], getGraphpp_p[1], getGraphpp_p[2], getGraphpp_p[3], getGraphpp_p[4], getGraphpp_p[5], histpp_p[0], histpp_p[1], histpp_p[2], histpp_p[3], histpp_p[4], histpp_p[5], 1, Tight);
    
    histpp_p[5]->SetMarkerStyle(25);
    histpp_p[5]->SetFillColor(0);

    drawHistToPTStack(histpp_p[0], kBlue - 9, "E1 HIST");
    drawHistToPTStack(histpp_p[1], kYellow - 9, "E1 HIST SAME");
    drawHistToPTStack(histpp_p[2], kOrange + 1, "E1 HIST SAME");
    drawHistToPTStack(histpp_p[3], kGreen + 3, "E1 HIST SAME");
    drawHistToPTStack(histpp_p[4], kRed + 1, "E1 HIST SAME");
    
    histpp_p[5]->DrawCopy("SAME E1");
    
    
    Int_t sys5A[5] = {0, 3, 3, 3, 3};
    Int_t sys5C[5] = {0, 4, 8, 11, 14};
    Int_t sys5NC[5] = {0, 1, 7, 8, 10};
   

    for(Int_t iter = 1; iter <= histpp_p[5]->GetNbinsX(); iter++){
      Float_t yVal = histpp_p[5]->GetBinContent(iter);
      TLine* l = 0;
      
      if(strcmp(CNC, "") == 0){
	
	std::cout <<"Iter: " << iter << ", " << yVal << ", " << yVal*.04 << std::endl;

	l = new TLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter+1) - .01, yVal - TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal));
	l->SetLineColor(1);
	l->Draw();
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal) + 2);
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal) + 2);
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter+1) -.01, yVal + TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal));
	
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal) - 2, histpp_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal));
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal) - 2, histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys5A[iter]*sys5A[iter] + .04*yVal*.04*yVal));
	
      }
      else if(strcmp(CNC, "C") == 0){
	
	l = new TLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter+1) - .01, yVal - TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal));
	l->SetLineColor(1);
	l->Draw();
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal) + 2);
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal) + 2);
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter+1) -.01, yVal + TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal));
	
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal) - 2, histpp_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal));
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal) - 2, histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys5C[iter]*sys5C[iter] + .04*yVal*.04*yVal));
	
      }
      else{
	
	l = new TLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter+1) - .01, yVal - TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal));
	l->SetLineColor(1);
	l->Draw();
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter) + .01, yVal - TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal) + 2);
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal - TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal) + 2);
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal), histpp_p[5]->GetBinLowEdge(iter+1) -.01, yVal + TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal));
	
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal) - 2, histpp_p[5]->GetBinLowEdge(iter) + .01, yVal + TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal));
	l->DrawLine(histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal) - 2, histpp_p[5]->GetBinLowEdge(iter + 1) - .01, yVal + TMath::Sqrt(sys5NC[iter]*sys5NC[iter] + .04*yVal*.04*yVal));
	
      }
    }
          
    zeroLine_p->Draw();
    
    label1_p->DrawLatex(.24, overCoord[0], Form("%s", overLabel[0]));
    label1_p->DrawLatex(.24, overCoord[1], "CMS Preliminary");
    //    label1_p->DrawLatex(.24, .78, "Gluon Jets");

  
    TLegend* leg2 = new TLegend(.2, .1, .55, .4);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(43);
    leg2->SetTextSize(14);

    hist1_p[5]->SetMarkerStyle(5);

    if(montecarlo){
      leg2->SetTextSize(10);
      leg2->AddEntry(histpp_p[5], "PYTHIA", "p");
      leg2->AddEntry(hist1_p[5], "PYTHIA+HYDJET", "p");
    }
    else{
      leg2->AddEntry(histpp_p[5], "pp", "p");
      leg2->AddEntry(hist1_p[5], "PbPb", "p");
    }

    leg2->Draw("SAME");

    TH1::SetDefaultSumw2();

    float histBound;

    if(strcmp("", CNC) == 0)
      histBound = 60;
    else
      histBound = 60;
   
    for(Int_t histIter = 0; histIter < 6; histIter++){
      niceTH1(hist1_p[histIter], histBound, -histBound, 505, 406);
      niceTH1(hist2_p[histIter], histBound, -histBound, 505, 406);
      niceTH1(hist3_p[histIter], histBound, -histBound, 505, 406);
      niceTH1(hist4_p[histIter], histBound, -histBound, 505, 406);
    }
 
    makeSubHistForPtStack(getGraph1_p[0], getGraph1_p[1], getGraph1_p[2], getGraph1_p[3], getGraph1_p[4], getGraph1_p[5], hist1_p[0], hist1_p[1], hist1_p[2], hist1_p[3], hist1_p[4], hist1_p[5], getGraphpp_p[0], getGraphpp_p[1], getGraphpp_p[2], getGraphpp_p[3], getGraphpp_p[4], getGraphpp_p[5], 7);
       
    hist1_p[5]->SetMarkerStyle(5);

    if(strcmp(CNC, "") == 0)
      profPanel_p->cd(7);
    else
      profPanel_p->cd(5);
    
    drawHistToPTStack(hist1_p[0], kBlue - 9, "E1 HIST");
    drawHistToPTStack(hist1_p[1], kYellow - 9, "E1 HIST SAME");
    drawHistToPTStack(hist1_p[2], kOrange + 1, "E1 HIST SAME");
    drawHistToPTStack(hist1_p[3], kGreen + 3, "E1 HIST SAME");
    drawHistToPTStack(hist1_p[4], kRed + 1, "E1 HIST SAME");

    hist1_p[5]->SetMarkerStyle(5);

    //    hist1_p[5]->DrawCopy("SAME E1");    

    zeroLine_p->Draw();
    
    label1_p->DrawLatex(.22, overCoord[0], Form("%s", overLabel[3]));

    if(strcmp(CNC, "") == 0)
      label1_p->DrawLatex(.22, overCoord[1], "50-100%");    
    else
      label1_p->DrawLatex(.22, overCoord[1], "30-100%");    

    makeSubHistForPtStack(getGraph2_p[0], getGraph2_p[1], getGraph2_p[2], getGraph2_p[3], getGraph2_p[4], getGraph2_p[5], hist2_p[0], hist2_p[1], hist2_p[2], hist2_p[3], hist2_p[4], hist2_p[5], getGraphpp_p[0], getGraphpp_p[1], getGraphpp_p[2], getGraphpp_p[3], getGraphpp_p[4], getGraphpp_p[5], 8);
    
    if(strcmp(CNC, "") == 0)
      profPanel_p->cd(8);
    else
      profPanel_p->cd(6);    

    drawHistToPTStack(hist2_p[0], kBlue - 9, "E1 HIST");
    drawHistToPTStack(hist2_p[1], kYellow - 9, "E1 HIST SAME");
    drawHistToPTStack(hist2_p[2], kOrange + 1, "E1 HIST SAME");
    drawHistToPTStack(hist2_p[3], kGreen + 3, "E1 HIST SAME");
    drawHistToPTStack(hist2_p[4], kRed + 1, "E1 HIST SAME");
    
    //    hist2_p[5]->DrawCopy("SAME E1");
    
    zeroLine_p->Draw();
    
    label1_p->DrawLatex(.22, overCoord[0], Form("%s", overLabel[3]));
    if(strcmp(CNC, "") == 0)
      label1_p->DrawLatex(.22, overCoord[1], "30-50%");    
    else
      label1_p->DrawLatex(.22, overCoord[1], "0-30%");    


    if(strcmp("", CNC) == 0){
      makeSubHistForPtStack(getGraph3_p[0], getGraph3_p[1], getGraph3_p[2], getGraph3_p[3], getGraph3_p[4], getGraph3_p[5], hist3_p[0], hist3_p[1], hist3_p[2], hist3_p[3], hist3_p[4], hist3_p[5], getGraphpp_p[0], getGraphpp_p[1], getGraphpp_p[2], getGraphpp_p[3], getGraphpp_p[4], getGraphpp_p[5], 9);

      profPanel_p->cd(9);
      
      drawHistToPTStack(hist3_p[0], kBlue - 9, "E1 HIST");
      drawHistToPTStack(hist3_p[1], kYellow - 9, "E1 HIST SAME");
      drawHistToPTStack(hist3_p[2], kOrange + 1, "E1 HIST SAME");
      drawHistToPTStack(hist3_p[3], kGreen + 3, "E1 HIST SAME");
      drawHistToPTStack(hist3_p[4], kRed + 1, "E1 HIST SAME");
    
      //      hist3_p[5]->DrawCopy("SAME E1");
     
      zeroLine_p->Draw();
      
      label1_p->DrawLatex(.22, overCoord[0], Form("%s", overLabel[3]));
      label1_p->DrawLatex(.22, overCoord[1], "10-30%");    
      
      makeSubHistForPtStack(getGraph4_p[0], getGraph4_p[1], getGraph4_p[2], getGraph4_p[3], getGraph4_p[4], getGraph4_p[5], hist4_p[0], hist4_p[1], hist4_p[2], hist4_p[3], hist4_p[4], hist4_p[5], getGraphpp_p[0], getGraphpp_p[1], getGraphpp_p[2], getGraphpp_p[3], getGraphpp_p[4], getGraphpp_p[5], 10);
      
      profPanel_p->cd(10);
      
      drawHistToPTStack(hist4_p[0], kBlue - 9, "E1 HIST");
      drawHistToPTStack(hist4_p[1], kYellow - 9, "E1 HIST SAME");
      drawHistToPTStack(hist4_p[2], kOrange + 1, "E1 HIST SAME");
      drawHistToPTStack(hist4_p[3], kGreen + 3, "E1 HIST SAME");
      drawHistToPTStack(hist4_p[4], kRed + 1, "E1 HIST SAME");
    
      //      hist4_p[5]->DrawCopy("SAME E1");

      zeroLine_p->Draw();
      
      label1_p->DrawLatex(.22, overCoord[0], Form("%s", overLabel[3]));
      label1_p->DrawLatex(.22, overCoord[1], "0-10%");
    }
  }
  panelFile_p->cd();

  profPanel_p->Write();

  if(strcmp(ppFileTag, "") == 0)
    claverCanvasSaving(profPanel_p, Form("../pngDir/%s%sImbAsymm%s%s%s%sPTStack_%s_%s", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), "pdf");
  else
    claverCanvasSaving(profPanel_p, Form("../pngDir/%s%sImbAsymm%s%s%s%sPTStackPP_%s_%s", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), "pdf");

  delete label1_p;
  delete label2_p;
  delete zeroLine_p;
  delete profPanel_p;
  delete histAxis_p;
  for(Int_t histIter = 0; histIter < 6; histIter++){
    delete hist1_p[histIter];
    delete hist2_p[histIter];
    delete hist3_p[histIter];
    delete hist4_p[histIter];
    delete histpp_p[histIter];
  }


  panelFile_p->Close();
  if(strcmp(ppName, "") !=0){
    ppFile_p->Close();
    delete ppFile_p;
  }

  delete panelFile_p;
  delete leg;
}


void makeImbAsymmPtStack_RECOGEN(const char* HIName, const char* ppName, const char* ppFileTag, const char* perpProj, const char* GLN, const char* Corr = "",  const char* CNC = "", Bool_t montecarlo = false)
{
  TFile* panelFile_p = new TFile(HIName, "UPDATE");
  TFile* ppFile_p = new TFile(ppName, "UPDATE");

  TH1::SetDefaultSumw2();

  TGraphErrors* getGraph1Reco_p[6];
  TGraphErrors* getGraph2Reco_p[6];
  TGraphErrors* getGraph3Reco_p[6];
  TGraphErrors* getGraph4Reco_p[6];

  TGraphErrors* getGraph1Gen_p[6];
  TGraphErrors* getGraph2Gen_p[6];
  TGraphErrors* getGraph3Gen_p[6];
  TGraphErrors* getGraph4Gen_p[6];

  TGraphErrors* getGraphppReco_p[6];
  TGraphErrors* getGraphppGen_p[6];

  TH1F* hist1_p[6];
  TH1F* hist2_p[6];
  TH1F* hist3_p[6];
  TH1F* hist4_p[6];
  TH1F* histpp_p[6];

  TH1F* histAxis_p;

  Float_t binArrayX[5] = {.00, .11, .22, .33, .50};

  const char* namePT[6] = {"0_1", "1_2", "2_4", "4_8", "8_100", "F"};

  histAxis_p = new TH1F("histAxis", "histAxis", 4, binArrayX);
  niceTH1(histAxis_p, 40, -40, 505, 406);

  histAxis_p->SetBinContent(1, 60);
  histAxis_p->SetBinContent(2, 60);
  histAxis_p->SetBinContent(3, 60);
  histAxis_p->SetBinContent(4, 60);

  for(Int_t histIter = 0; histIter < 6; histIter++){

    hist1_p[histIter] = new TH1F(Form("hist1%s_p", namePT[histIter]), Form("hist1%s_p", namePT[histIter]), 4, binArrayX);
    hist2_p[histIter] = new TH1F(Form("hist2%s_p", namePT[histIter]), Form("hist2%s_p", namePT[histIter]), 4, binArrayX); 
    hist3_p[histIter] = new TH1F(Form("hist3%s_p", namePT[histIter]), Form("hist3%s_p", namePT[histIter]), 4, binArrayX);
    hist4_p[histIter] = new TH1F(Form("hist4%s_p", namePT[histIter]), Form("hist4%s_p", namePT[histIter]), 4, binArrayX);
    histpp_p[histIter] = new TH1F(Form("histpp%s_p", namePT[histIter]), Form("histpp%s_p", namePT[histIter]), 4, binArrayX);

    niceTH1(hist1_p[histIter], 60, -60, 505, 406);
    niceTH1(hist2_p[histIter], 60, -60, 505, 406);
    niceTH1(hist3_p[histIter], 60, -60, 505, 406);
    niceTH1(hist4_p[histIter], 60, -60, 505, 406);
    niceTH1(histpp_p[histIter], 60, -60, 505, 406);

    if(strcmp(CNC, "") == 0){
      getGraph1Reco_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("rVsCaloImbAsymm%s%s%s%s_50100_%s_%s_g",  perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
      getGraph1Gen_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("gTImbAsymm%s%s%s_50100_%s_%s_g", perpProj, CNC, namePT[histIter], GLN, fileTag1));
    }
    else{
      getGraph1Reco_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("rVsCaloImbAsymm%s%s%s%s_30100_%s_%s_g", perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
      getGraph1Gen_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("gTImbAsymm%s%s%s_30100_%s_%s_g", perpProj, CNC, namePT[histIter], GLN, fileTag1));
    }    
  }

  makeSubHistForPtStack(getGraph1Reco_p[0], getGraph1Reco_p[1], getGraph1Reco_p[2], getGraph1Reco_p[3], getGraph1Reco_p[4], getGraph1Reco_p[5], hist1_p[0], hist1_p[1], hist1_p[2], hist1_p[3], hist1_p[4], hist1_p[5], getGraph1Gen_p[0], getGraph1Gen_p[1], getGraph1Gen_p[2], getGraph1Gen_p[3], getGraph1Gen_p[4], getGraph1Gen_p[5], 2);

  hist1_p[5]->SetMarkerStyle(5);

  TCanvas* profPanel_p;
  
  if(strcmp(CNC, "") == 0){
    profPanel_p = new TCanvas(Form("rMinGVsCaloImbAsymm%s%s%sPTStackPP_%s_%s_c", perpProj, CNC, Corr, GLN, fileTag1), Form("rMinGVsCaloImbAsymm%s%s%sPTStackPP_%s_%s_c", perpProj, CNC, Corr, GLN, fileTag1), 1000, 500);
    profPanel_p->Divide(5, 2, 0, 0);
    std::cout << "FivePanel Init" << std::endl;
  }
  else{
    profPanel_p = new TCanvas(Form("rMinGVsCaloImbAsymm%s%s%sPTStackPP_%s_%s_c", perpProj, CNC, Corr, GLN, fileTag1), Form("rMinGVsCaloImbAsymm%s%s%sPTStackPP_%s_%s_c",  perpProj, CNC, Corr, GLN, fileTag1), 600, 500);
    profPanel_p->Divide(3, 2, 0, 0);
    std::cout << "ThreePanel Init" << std::endl;
  }

  const char* qsquare;

  if(strcmp(fileTag1, "DataE") == 0)
    qsquare = Form(" ");
  else if(strcmp(fileTag1, "Di80i") == 0)
    qsquare = Form("(#hat{p}_{T} > 80)");
  else if(strcmp(fileTag1, "Di120b") == 0)
    qsquare = Form("(#hat{p}_{T} > 120)");


  TLegend* leg;
  leg = new TLegend(0.2, 0.60, 0.95, 0.95, Form("Corr. Trk - Truth #slash{p}_{T}^{||}  v. A_{J}   %s", qsquare));


  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSizePixels(12);
  leg->SetBorderSize(0);

  profPanel_p->cd(2);

  drawHistToPTStack(hist1_p[0], kBlue - 9, "E1 HIST");
  leg->AddEntry(hist1_p[0], ".5 < p_{T} < 1", "F");
  drawHistToPTStack(hist1_p[1], kYellow - 9, "E1 HIST SAME");
  leg->AddEntry(hist1_p[1], "1 < p_{T} < 2", "F");
  drawHistToPTStack(hist1_p[2], kOrange + 1, "E1 HIST SAME");
  leg->AddEntry(hist1_p[2], "2 < p_{T} < 4", "F");
  drawHistToPTStack(hist1_p[3], kGreen + 3, "E1 HIST SAME");
  leg->AddEntry(hist1_p[3], "4 < p_{T} < 8", "F");
  drawHistToPTStack(hist1_p[4], kRed + 1, "E1 HIST SAME");
  leg->AddEntry(hist1_p[4], "8 < p_{T}", "F");

  hist1_p[5]->DrawCopy("SAME E1");

  if(strcmp("", CNC) == 0)
    profPanel_p->cd(6);
  else
    profPanel_p->cd(4);
  
  histAxis_p->DrawCopy("HIST");
  leg->Draw("SAME");


  profPanel_p->cd(2);

  TLine* zeroLine_p = new TLine(0., 0., 0.5, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw();

  TLatex* label1_p = new TLatex();
  label1_p->SetNDC();
  label1_p->SetTextFont(42);
  label1_p->SetTextSizePixels(18);

  TLatex* label2_p = new TLatex();
  label2_p->SetNDC();
  label2_p->SetTextFont(42);
  label2_p->SetTextSizePixels(8);

  label1_p->DrawLatex(.22, .90, "#Delta_{Reco,Gen}(PYTHIA+HYDJET)");
  if(strcmp("", CNC) == 0)
    label1_p->DrawLatex(.22, .84, "50-100%");
  else
    label1_p->DrawLatex(.22, .84, "30-100%");

  profPanel_p->cd(3);

  for(Int_t histIter = 0; histIter < 6; histIter++){
    if(strcmp(CNC, "") == 0){
      getGraph2Reco_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("rVsCaloImbAsymm%s%s%s%s_3050_%s_%s_g", perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
      getGraph2Gen_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("gTImbAsymm%s%s%s_3050_%s_%s_g", perpProj, CNC, namePT[histIter], GLN, fileTag1));
    }
    else{
      getGraph2Reco_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("rVsCaloImbAsymm%s%s%s%s_030_%s_%s_g", perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
      getGraph2Gen_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("gTImbAsymm%s%s%s_030_%s_%s_g", perpProj, CNC, namePT[histIter], GLN, fileTag1));
    }
  }

  makeSubHistForPtStack(getGraph2Reco_p[0], getGraph2Reco_p[1], getGraph2Reco_p[2], getGraph2Reco_p[3], getGraph2Reco_p[4], getGraph2Reco_p[5], hist2_p[0], hist2_p[1], hist2_p[2], hist2_p[3], hist2_p[4], hist2_p[5], getGraph2Gen_p[0], getGraph2Gen_p[1], getGraph2Gen_p[2], getGraph2Gen_p[3], getGraph2Gen_p[4], getGraph2Gen_p[5], 3);

  hist2_p[5]->SetMarkerStyle(5);

  drawHistToPTStack(hist2_p[0], kBlue - 9, "E1 HIST");
  drawHistToPTStack(hist2_p[1], kYellow - 9, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[2], kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[3], kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[4], kRed + 1, "E1 HIST SAME");

  hist2_p[5]->DrawCopy("SAME E1");

  zeroLine_p->Draw();

  label1_p->DrawLatex(.22, .90, "#Delta_{Reco,Gen}(PYTHIA+HYDJET)");
  if(strcmp("", CNC) == 0)
    label1_p->DrawLatex(.22, .84, "30-50%");
  else
    label1_p->DrawLatex(.22, .84, "0-30%");

    

  if(strcmp(CNC, "") == 0)
    profPanel_p->cd(6);
  else
    profPanel_p->cd(4);

  label2_p->DrawLatex(.24, .52, "Jet R = .3");
  label2_p->DrawLatex(.24, .47, "p_{T,1} > 120, p_{T,2} > 50 GeV/c");
  if(strcmp(CNC, "") == 0){
    label2_p->DrawLatex(.24, .42, "Jet |#eta|_{1}, |#eta|_{2} < 1.6");
    label2_p->DrawLatex(.24, .37, "Jet #Delta #phi > 5#pi/6");
  }
  else{
    label2_p->DrawLatex(.24, .42, "Jet |#eta|_{1}, |#eta|_{2} < 1.6");
    label2_p->DrawLatex(.24, .37, "Jet #Delta #phi > 5#pi/6");
  }

  if(strcmp(CNC, "C") == 0)
    label2_p->DrawLatex(.26, .32, Form("In-Cone, #Delta R < .8"));
  else if(strcmp(CNC, "NC") == 0)
    label2_p->DrawLatex(.26, .32, Form("Out-of-Cone, .8 < #Delta R"));
  else if(strcmp(CNC, "NCCut") == 0)
    label2_p->DrawLatex(.26, .32, Form("Out-of-Cone, .8 < #Delta R < #pi/2"));


  if(strcmp(CNC, "") == 0){
    profPanel_p->cd(4);

    for(Int_t histIter = 0; histIter < 6; histIter++){
      getGraph3Reco_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("rVsCaloImbAsymm%s%s%s%s_1030_%s_%s_g", perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
      getGraph3Gen_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("gTImbAsymm%s%s%s_1030_%s_%s_g", perpProj, CNC, namePT[histIter], GLN, fileTag1));
    }
    
    makeSubHistForPtStack(getGraph3Reco_p[0], getGraph3Reco_p[1], getGraph3Reco_p[2], getGraph3Reco_p[3], getGraph3Reco_p[4], getGraph3Reco_p[5], hist3_p[0], hist3_p[1], hist3_p[2], hist3_p[3], hist3_p[4], hist3_p[5], getGraph3Gen_p[0], getGraph3Gen_p[1], getGraph3Gen_p[2], getGraph3Gen_p[3], getGraph3Gen_p[4], getGraph3Gen_p[5], 4);

    hist3_p[5]->SetMarkerStyle(5);
    
    drawHistToPTStack(hist3_p[0], kBlue - 9, "E1 HIST");
    drawHistToPTStack(hist3_p[1], kYellow - 9, "E1 HIST SAME");
    drawHistToPTStack(hist3_p[2], kOrange + 1, "E1 HIST SAME");
    drawHistToPTStack(hist3_p[3], kGreen + 3, "E1 HIST SAME");
    drawHistToPTStack(hist3_p[4], kRed + 1, "E1 HIST SAME");
    
    hist3_p[5]->DrawCopy("SAME E1");
    
    zeroLine_p->Draw();
    
    label1_p->DrawLatex(.22, .90, "#Delta_{Reco,Gen}(PYTHIA+HYDJET)");
    label1_p->DrawLatex(.22, .84, "10-30%");  
    
    profPanel_p->cd(5);
    
    for(Int_t histIter = 0; histIter < 6; histIter++){
      getGraph4Reco_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("rVsCaloImbAsymm%s%s%s%s_010_%s_%s_g", perpProj, CNC, namePT[histIter], Corr, GLN, fileTag1));
      getGraph4Gen_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("gTImbAsymm%s%s%s_010_%s_%s_g", perpProj, CNC, namePT[histIter], GLN, fileTag1));
    }
    
    makeSubHistForPtStack(getGraph4Reco_p[0], getGraph4Reco_p[1], getGraph4Reco_p[2], getGraph4Reco_p[3], getGraph4Reco_p[4], getGraph4Reco_p[5], hist4_p[0], hist4_p[1], hist4_p[2], hist4_p[3], hist4_p[4], hist4_p[5], getGraph4Gen_p[0], getGraph4Gen_p[1], getGraph4Gen_p[2], getGraph4Gen_p[3], getGraph4Gen_p[4], getGraph4Gen_p[5], 4);
    
    hist4_p[5]->SetMarkerStyle(5);

    drawHistToPTStack(hist4_p[0], kBlue - 9, "E1 HIST");
    drawHistToPTStack(hist4_p[1], kYellow - 9, "E1 HIST SAME");
    drawHistToPTStack(hist4_p[2], kOrange + 1, "E1 HIST SAME");
    drawHistToPTStack(hist4_p[3], kGreen + 3, "E1 HIST SAME");
    drawHistToPTStack(hist4_p[4], kRed + 1, "E1 HIST SAME");
    
    hist4_p[5]->DrawCopy("SAME E1");
    
    zeroLine_p->Draw();
    
    label1_p->DrawLatex(.22, .90, "#Delta_{Reco,Gen}(PYTHIA+HYDJET)");
    label1_p->DrawLatex(.22, .84, "0-10%");
  }

  
  profPanel_p->cd(1);

  for(Int_t histIter = 0; histIter < 6; histIter++){
    getGraphppReco_p[histIter] = (TGraphErrors*)ppFile_p->Get(Form("rVsCaloImbAsymm%s%s%s%s_%s_%s_g", perpProj, CNC, namePT[histIter], Corr, GLN, ppFileTag));
    getGraphppGen_p[histIter] = (TGraphErrors*)ppFile_p->Get(Form("gTImbAsymm%s%s%s_%s_%s_g", perpProj, CNC, namePT[histIter], GLN, ppFileTag));
  }

  makeSubHistForPtStack(getGraphppReco_p[0], getGraphppReco_p[1], getGraphppReco_p[2], getGraphppReco_p[3], getGraphppReco_p[4], getGraphppReco_p[5], histpp_p[0], histpp_p[1], histpp_p[2], histpp_p[3], histpp_p[4], histpp_p[5], getGraphppGen_p[0], getGraphppGen_p[1], getGraphppGen_p[2], getGraphppGen_p[3], getGraphppGen_p[4], getGraphppGen_p[5], 1, true);
  
  histpp_p[5]->SetMarkerStyle(25);

  drawHistToPTStack(histpp_p[0], kBlue - 9, "E1 HIST");
  drawHistToPTStack(histpp_p[1], kYellow - 9, "E1 HIST SAME");
  drawHistToPTStack(histpp_p[2], kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(histpp_p[3], kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(histpp_p[4], kRed + 1, "E1 HIST SAME");
  
  histpp_p[5]->DrawCopy("SAME E1");
  
  zeroLine_p->Draw();
  label1_p->DrawLatex(.24, .90, "#Delta_{Reco,Gen}(PYTHIA)");
  label1_p->DrawLatex(.24, .84, "CMS Preliminary");
  //  label1_p->DrawLatex(.24, .78, "Gluon Jets");
  
  TLegend* leg2 = new TLegend(.2, .1, .55, .4);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(14);

  TH1F* dummyHist_p = new TH1F("dummy", "dummy", 4, binArrayX);
  dummyHist_p->SetMarkerStyle(5);

  if(montecarlo){
    leg2->SetTextSize(10);
    leg2->AddEntry(histpp_p[5], "PYTHIA", "p");
    leg2->AddEntry(dummyHist_p, "PYTHIA+HYDJET", "p");
    leg2->AddEntry(hist1_p[5], "(PYTHIA+HYDJET) - PYTHIA", "p");
  }
  else{
    leg2->AddEntry(histpp_p[5], "pp", "p");
    leg2->AddEntry(dummyHist_p, "PbPb", "p");
    leg2->AddEntry(hist1_p[5], "PbPb - pp", "p");
  }

  leg2->Draw("SAME");  

   
  if(strcmp(CNC, "") == 0)
    profPanel_p->cd(7);
  else
    profPanel_p->cd(5);
  

  hist1_p[5]->Add(histpp_p[5], -1);
  hist1_p[5]->SetMaximum(40);
  hist1_p[5]->SetMinimum(-40);

  hist1_p[5]->SetMarkerStyle(24);
  hist1_p[5]->DrawCopy("E1");
     
  zeroLine_p->Draw();
     
  label1_p->DrawLatex(.12, .90, "#Delta_{Reco,Gen}((PYTHIA+HYDJET) - PYTHIA)");
     
  if(strcmp(CNC, "") == 0)
    label1_p->DrawLatex(.22, .84, "50-100%, expanded axis range");    
  else
    label1_p->DrawLatex(.22, .84, "30-100%, expanded axis range");    
     
  
  if(strcmp(CNC, "") == 0)
    profPanel_p->cd(8);
  else
    profPanel_p->cd(6);    
     
  hist2_p[5]->Add(histpp_p[5], -1);  
  hist2_p[5]->SetMaximum(40);
  hist2_p[5]->SetMinimum(-40);
  hist2_p[5]->SetMarkerStyle(24);
  hist2_p[5]->DrawCopy("E1");
    
  zeroLine_p->Draw();
    
  label1_p->DrawLatex(.12, .90, "#Delta_{Reco,Gen}((PYTHIA+HYDJET) - PYTHIA)");
  if(strcmp(CNC, "") == 0)
    label1_p->DrawLatex(.22, .84, "30-50%");    
  else
    label1_p->DrawLatex(.22, .84, "0-30%");    


  if(strcmp("", CNC) == 0){
    
    profPanel_p->cd(9);

    hist3_p[5]->Add(histpp_p[5], -1);    
    hist3_p[5]->SetMaximum(40);
    hist3_p[5]->SetMinimum(-40);
    hist3_p[5]->SetMarkerStyle(24);
    hist3_p[5]->DrawCopy("E1");
      
    zeroLine_p->Draw();
      
    label1_p->DrawLatex(.12, .90, "#Delta_{Reco,Gen}((PYTHIA+HYDJET) - PYTHIA)");
    label1_p->DrawLatex(.22, .84, "10-30%");    
      
    profPanel_p->cd(10);
 
    hist4_p[5]->Add(histpp_p[5], -1);  
    hist4_p[5]->SetMaximum(40);
    hist4_p[5]->SetMinimum(-40);
    hist4_p[5]->SetMarkerStyle(24);
    hist4_p[5]->DrawCopy("E1");
    
    zeroLine_p->Draw();

    label1_p->DrawLatex(.12, .90, "#Delta_{Reco,Gen}((PYTHIA+HYDJET) - PYTHIA)");    
    label1_p->DrawLatex(.22, .84, "0-10%");
  }
  
  panelFile_p->cd();
  profPanel_p->Write();

  claverCanvasSaving(profPanel_p, Form("../pngDir/rMinGVsCaloImbAsymm%s%s%sPTStackPP_%s_%s", perpProj, CNC, Corr, GLN, fileTag1), "pdf");

  delete label1_p;
  delete label2_p;
  delete zeroLine_p;
  delete profPanel_p;
  delete histAxis_p;

  for(Int_t histIter = 0; histIter < 6; histIter++){
    delete hist1_p[histIter];
    delete hist2_p[histIter];
    delete hist3_p[histIter];
    delete hist4_p[histIter];
    delete histpp_p[histIter];
  }


  panelFile_p->Close();
  if(strcmp(ppName, "") !=0){
    ppFile_p->Close();
    delete ppFile_p;
  }

  delete panelFile_p;
  delete leg;
}






void cfmDiJet_PtImbPlots(const char* inName, const char* outName, Bool_t montecarlo = false, const char* ppName = "", const char* ppFileTag = "", Bool_t Veto = false)
{
  //Some personal shorthand for myself, can be ignored as long as you feed it something for filetag

  TH1::SetDefaultSumw2();

  if(!strcmp(inName, Di80a)){
    std::cout << Di80a << std::endl;
    fileTag1 = "Di80a";
  }
  else if(!strcmp(inName, Di80b)){
    std::cout << Di80b << std::endl;
    fileTag1 = "Di80b";
  }
  else if(!strcmp(inName, Di100a)){
    std::cout << Di100a << std::endl;
    fileTag1 = "Di100a";
  }
  else if(!strcmp(inName, Di100b)){
    std::cout << Di100b << std::endl;
    fileTag1 = "Di100b";
  }
  else if(!strcmp(inName, Di120a)){
    std::cout << Di120a << std::endl;
    fileTag1 = "Di120a";
  }
  else if(!strcmp(inName, Di120b)){
    std::cout << Di120b << std::endl;
    fileTag1 = "Di120b";
  }
  else if(!strcmp(inName, EmDi80a)){
    std::cout << EmDi80a << std::endl;
    fileTag1 = "EmDi80a";
  }
  else if(!strcmp(inName, Di80c)){
    std::cout << Di80c << std::endl;
    fileTag1 = "Di80c";
  }
  else if(!strcmp(inName, Di80d)){
    std::cout << Di80d << std::endl;
    fileTag1 = "Di80d";
  }
  else if(!strcmp(inName, Di80e)){
    std::cout << Di80e << std::endl;
    fileTag1 = "Di80e";
  }
  else if(!strcmp(inName, Di80f)){
    std::cout << Di80f << std::endl;
    fileTag1 = "Di80f";
  }
  else if(!strcmp(inName, Di80g)){
    std::cout << Di80g << std::endl;
    fileTag1 = "Di80g";
  }
  else if(!strcmp(inName, Di80h)){
    std::cout << Di80h << std::endl;
    fileTag1 = "Di80h";
  }
  else if(!strcmp(inName, Di80i)){
    std::cout << Di80i << std::endl;
    fileTag1 = "Di80i";
  }
  else if(!strcmp(inName, DataA)){
    std::cout << DataA << std::endl;
    fileTag1 = "DataA";
  }
  else if(!strcmp(inName, DataB)){
    std::cout << DataB << std::endl;
    fileTag1 = "DataB";
  }
  else if(!strcmp(inName, DataC)){
    std::cout << DataC << std::endl;
    fileTag1 = "DataC";
  }
  else if(!strcmp(inName, DataD)){
    std::cout << DataD << std::endl;
    fileTag1 = "DataD";
  }
  else if(!strcmp(inName, DataE)){
    std::cout << DataE << std::endl;
    fileTag1 = "DataE";
  }

  std::cout << "Filetag1 is: " << fileTag1 << std::endl;

  inFile1_p = new TFile(inName, "READ");
  inTree_p = (TTree*)inFile1_p->Get("jetTree");
  inTree_p->AddFriend("trackTree");

  Int_t jetAlgMax = 4;

  if(montecarlo){
    inTree_p->AddFriend("genTree");
    jetAlgMax = 5;
  }

  const char* corr[2] = {"", "Corr"};
  const char* CNC[4] = {"", "C", "NC", "NCCut"};
  const char* ptBins[5] = {"0_1", "1_2", "2_4", "4_8", "8_100"};

  Int_t centLow[6] = {0, 20, 60, 100, 0, 60};
  Int_t centHi[6] = {19, 59, 99, 199, 59, 199};

  for(Int_t algIter = 0; algIter < jetAlgMax; algIter++){
    if(algIter != 3 && algIter != 4) continue;

    std::cout << "Algorithm: " << algType[algIter] << std::endl;

    if(Veto && algIter !=4)
      fullVeto = Form("AlgLeadRefPt[%d] > 0 && AlgSubLeadRefPt[%d] > 0 && AlgLeadRefPt[%d] > AlgSubLeadRefPt[%d]", algIter, algIter, algIter, algIter);
    if(Veto && algIter == 4)
      fullVeto = "";
    //      fullVeto = Form("AlgThirdJtPt[4] < 50");

    for(Int_t corrIter = 0; corrIter < 2; corrIter++){
      for(Int_t CNCIter = 0; CNCIter < 4; CNCIter++){
	  
	for(Int_t centIter = 0; centIter < 6; centIter++){
	  makeImbAsymmGraph(inTree_p, outName, "r", algIter, "ProjA", CNC[CNCIter], "F", centLow[centIter], centHi[centIter], -60, 60, "N", corr[corrIter], montecarlo);
	    
	  if(montecarlo && corrIter == 0){
	    makeImbAsymmGraph(inTree_p, outName, "g", algIter, "ProjA", CNC[CNCIter], "F", centLow[centIter], centHi[centIter], -60, 60, "N", "", montecarlo);
	  }
   
	  for(Int_t ptBinIter = 0; ptBinIter < 5; ptBinIter++){
	    makeImbAsymmGraph(inTree_p, outName, "r", algIter, "ProjA", CNC[CNCIter], ptBins[ptBinIter], centLow[centIter], centHi[centIter], -60, 60, "N", corr[corrIter], montecarlo);

	    if(montecarlo && corrIter == 0){
	      makeImbAsymmGraph(inTree_p, outName, "g", algIter, "ProjA", CNC[CNCIter], ptBins[ptBinIter], centLow[centIter], centHi[centIter], -60, 60, "N", "", montecarlo);

	    }
	  }
	}	
      }
	

      makeImbAsymmPtStack(outName, "r", algIter, "ProjA", "N", corr[corrIter], "", "", ppName, ppFileTag, montecarlo);
      
      makeImbAsymmPtStack(outName, "r", algIter, "ProjA", "N", corr[corrIter], "", "C", ppName, ppFileTag, montecarlo);
      
      makeImbAsymmPtStack(outName, "r", algIter, "ProjA", "N", corr[corrIter], "", "NC", ppName, ppFileTag, montecarlo);
      
      makeImbAsymmPtStack(outName, "r", algIter, "ProjA", "N", corr[corrIter], "", "NCCut", ppName, ppFileTag, montecarlo);


    }


    if(montecarlo){

      makeImbAsymmPtStack(outName, "g", algIter, "ProjA", "N", "", "", "", ppName, ppFileTag, montecarlo);      
      makeImbAsymmPtStack(outName, "g", algIter, "ProjA", "N", "", "", "C", ppName, ppFileTag, montecarlo);
      makeImbAsymmPtStack(outName, "g", algIter, "ProjA", "N", "", "", "NC", ppName, ppFileTag, montecarlo);
      makeImbAsymmPtStack(outName, "g", algIter, "ProjA", "N", "", "", "NCCut", ppName, ppFileTag, montecarlo);
     
    }
  
  }
  if(montecarlo){
    makeImbAsymmPtStack_RECOGEN(outName, ppName, ppFileTag, "ProjA", "N", "Corr", "", montecarlo);
    makeImbAsymmPtStack_RECOGEN(outName, ppName, ppFileTag, "ProjA", "N", "Corr", "C", montecarlo);
    makeImbAsymmPtStack_RECOGEN(outName, ppName, ppFileTag, "ProjA", "N", "Corr", "NC", montecarlo);
  }

  inFile1_p->Close();
  delete inFile1_p;
  return;
}
