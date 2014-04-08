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

TTree* inFile2_p = 0;

TFile* outFile_p = 0;

const char* algType[5] = {"PuPF", "PuCalo", "VsPF", "VsCalo", "T"};
const char* algType_pp[5] = {"PuPF", "PuCalo", "PF", "Calo", "T"};
TCut thirdJtVeto = "";

//Float_t setDelPhiCut = 5*TMath::Pi()/6;
Float_t setDelPhiCut = 0;

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

const char* Di120a = "HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHat_v0_0_CFMSKIM_20140323.root";

const char* Di120b = "HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHat_v0_0_CFMSKIM.root";

const char* EmDi80a = "PbPb_pythiaHYDJET_forest_EmEnrichedDijet80_CFMSKIM.root";

const char* DataA = "Track8_Jet17_GR_R_53_LV6_SUB_0_CFMSKIM.root";

const char* DataB = "hiForest_Jet80or95_GR_R_53_LV6_02Mar2014_1300CET_Track8_Jet15_0_1200k_CFMSKIM.root";

const char* DataC = "hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21_0_700k_CFMSKIM.root";

const char* DataD = "hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21_0_1200k_CFMSKIM.root";

const char* DataE = "hiForest_Jet80or95_GR_R_53_LV6_03Mar2014_1600CET_CMSSW_5_3_16_merged_0_CFMSKIM.root";

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


void niceTH1(TH1F* uglyTH1, float max , float min, float ndivX, float ndivY, Bool_t norm = true)
{
  if(norm)
    handsomeTH1N(uglyTH1);
  else 
    handsomeTH1(uglyTH1);

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


TCut makeAsymmCut(Int_t setNum, Float_t asymmLow, Float_t asymmHi)
{
  TCut asymmCut = "";

  if(setNum > 4 || setNum < 0){
    std::cout << "makeAsymmCut: setNum must be between 0-4, empty cut returned" << std::endl;
    return asymmCut;
  }

  if(asymmLow >= 0 && asymmHi >= asymmLow && asymmHi <= 1)
    asymmCut = Form("AlgJtAsymm[%d] > %f && AlgJtAsymm[%d] < %f ", setNum, asymmLow, setNum, asymmHi);
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


void makeCentHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, const char* fileTag2 = "")
{
  inFile_p->cd();

  const char* title = Form("hiBin_%s_%s", algType[setNum], fileTag1);

  TH1F* centHist_p;

  TString name = Form("%s_h(200, -.5, 199.5)", title);

  TCut setCut = makeSetCut(setNum);
  TCut jtDelPhiCut = makeDelPhiCut(setNum, setDelPhiCut);

  std::cout << setCut << ", " << jtDelPhiCut << std::endl;

  getTree_p->Project(name, "hiBin", setCut && jtDelPhiCut);
  
  centHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  centHist_p->Sumw2();
  centHist_p->Scale(1./centHist_p->GetEntries());
  centHist_p->SetXTitle(title);

  outFile_p = new TFile(outName, "UPDATE");
  centHist_p->Write(Form("%s_h", title));
  if(strcmp(fileTag2, "") != 0 && setNum < 4){
    TH1F* divHist_p = (TH1F*)outFile_p->Get(Form("hiBin_%s_%s_h", algType[setNum], fileTag2));
    divHist_p->Sumw2();
    divHist_p->Divide(centHist_p);
    divHist_p->SetXTitle(Form("hiBin_%s_DataOverMC_h", algType[setNum]));
    divHist_p->Write(Form("hiBin_%s_DataOverMC_h", algType[setNum]));
  }
  outFile_p->Close();

  delete outFile_p;
  return;
}


void makeAsymmHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Float_t histLow, Float_t histHi, Int_t centLow, Int_t centHi, Bool_t montecarlo = false)
{
  inFile_p->cd();

  const char* title = Form("%sAsymm_%d%d_%s", algType[setNum], (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag1);

  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %f, %f)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, 5*TMath::Pi()/6);
  TCut etaCut = makeEtaCut(setNum, 0.5);


  std::cout << delPhiCut << std::endl;

  if(montecarlo)
    getTree_p->Project(name, Form("AlgJtAsymm[%d]", setNum), Form("setWeight[%d]", setNum)*(centCut && setCut && thirdJtVeto && delPhiCut && etaCut));
  else 
    getTree_p->Project(name, Form("AlgJtAsymm[%d]", setNum), centCut && setCut && thirdJtVeto && delPhiCut && etaCut);

  asymmHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  outFile_p = new TFile(outName, "UPDATE");

  asymmHist_p->Sumw2();

  niceTH1(asymmHist_p, .55, 0., 405, 506);

  asymmHist_p->SetYTitle("Event Fraction");
  asymmHist_p->SetXTitle("A_{J} = (p_{T,1} - p_{T,2})/(p_{T,1} + p_{T,2})");

  asymmHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
  return;
}


void makePtRatioHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Int_t histLow, Int_t histHi, Int_t centLow, Int_t centHi, Bool_t montecarlo = false)
{
  inFile_p->cd();

  const char* title = Form("%sRatio_%d%d_%s", algType[setNum], (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag1);

  TH1F* ratioHist_p;

  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, setDelPhiCut);
  TCut etaCut = makeEtaCut(setNum, 0.5);

  std::cout << delPhiCut << std::endl;

  if(montecarlo)
    getTree_p->Project(name, Form("AlgSubLeadJtPt[%d]/AlgLeadJtPt[%d]", setNum, setNum), Form("setWeight[%d]", setNum)*(centCut && setCut && thirdJtVeto && delPhiCut && etaCut));
  else 
    getTree_p->Project(name, Form("AlgSubLeadJtPt[%d]/AlgLeadJtPt[%d]", setNum, setNum), centCut && setCut && thirdJtVeto && delPhiCut && etaCut);

  ratioHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  outFile_p = new TFile(outName, "UPDATE");

  ratioHist_p->Sumw2();

  niceTH1(ratioHist_p, .25, 0., 405, 506);

  ratioHist_p->SetYTitle("Event Fraction");
  ratioHist_p->SetXTitle("p_{T,2}/p_{T,1}");

  ratioHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
  return;
}



void makeAsymmDelPhiHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t centLow, Int_t centHi, Bool_t montecarlo = false)
{
  inFile1_p->cd();

  const char* title = Form("%sAsymmDelPhi_%d%d_%s", algType[setNum], (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag1);
  TH1F* getValHist_p[6];

  TH1F* asymmDelPhiHist_p = new TH1F(Form("%s_h", title), Form("%s_h", title), 6, 2*TMath::Pi()/3, TMath::Pi());

  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 0.5);

  for(Int_t binIter = 0; binIter < 6; binIter++){
    TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f && AlgJtDelPhi[%d] < %f", setNum, 2*TMath::Pi()/3 + binIter*TMath::Pi()/18, setNum, 2*TMath::Pi()/3 + (binIter + 1)*TMath::Pi()/18);
    TString name = Form("%d_h(10, 0, 1)", binIter);

    if(montecarlo)
      getTree_p->Project(name, Form("AlgJtAsymm[%d]", setNum), Form("setWeight[%d]", setNum)*(centCut && setCut && delPhiCut && etaCut));
    else
      getTree_p->Project(name, Form("AlgJtAsymm[%d]", setNum), centCut && setCut && delPhiCut && etaCut);

    getValHist_p[binIter] = (TH1F*)inFile_p->Get(Form("%d_h", binIter));

    asymmDelPhiHist_p->SetBinContent(binIter + 1, getValHist_p[binIter]->GetMean());
    asymmDelPhiHist_p->SetBinError(binIter + 1, getValHist_p[binIter]->GetMeanError());
  }

  outFile_p = new TFile(outName, "UPDATE");

  asymmDelPhiHist_p->Sumw2();

  niceTH1(asymmDelPhiHist_p, .45, 0., 405, 506, false);

  asymmDelPhiHist_p->SetYTitle("<A_{J}>");
  asymmDelPhiHist_p->SetXTitle("#Delta #phi_{1,2}");

  asymmDelPhiHist_p->Write();
  outFile_p->Close();

  delete outFile_p;
  delete asymmDelPhiHist_p;
  return;
}



void makeRatioDelPhiHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t centLow, Int_t centHi, Bool_t montecarlo = false)
{
  inFile1_p->cd();

  const char* title = Form("%sRatioDelPhi_%d%d_%s", algType[setNum], (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag1);
  TH1F* getValHist_p[6];

  TH1F* ratioDelPhiHist_p = new TH1F(Form("%s_h", title), Form("%s_h", title), 6, 2*TMath::Pi()/3, TMath::Pi());

  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 0.5);

  for(Int_t binIter = 0; binIter < 6; binIter++){
    TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f && AlgJtDelPhi[%d] < %f", setNum, 2*TMath::Pi()/3 + binIter*TMath::Pi()/18, setNum, 2*TMath::Pi()/3 + (binIter + 1)*TMath::Pi()/18);
    TString name = Form("%d_h(10, 0, 1)", binIter);

    if(montecarlo)
      getTree_p->Project(name, Form("AlgSubLeadJtPt[%d]/AlgLeadJtPt[%d]", setNum, setNum), Form("setWeight[%d]", setNum)*(centCut && setCut && delPhiCut && etaCut));
    else
      getTree_p->Project(name, Form("AlgSubLeadJtPt[%d]/AlgLeadJtPt[%d]", setNum, setNum), centCut && setCut && delPhiCut && etaCut);

    getValHist_p[binIter] = (TH1F*)inFile_p->Get(Form("%d_h", binIter));

    ratioDelPhiHist_p->SetBinContent(binIter + 1, getValHist_p[binIter]->GetMean());
    ratioDelPhiHist_p->SetBinError(binIter + 1, getValHist_p[binIter]->GetMeanError());
  }

  outFile_p = new TFile(outName, "UPDATE");

  ratioDelPhiHist_p->Sumw2();

  niceTH1(ratioDelPhiHist_p, .85, .40, 405, 506, false);

  ratioDelPhiHist_p->SetYTitle("<p_{T,2}/p_{T,1}>");
  ratioDelPhiHist_p->SetXTitle("#Delta #phi_{1,2}");

  ratioDelPhiHist_p->Write();
  outFile_p->Close();

  delete outFile_p;
  delete ratioDelPhiHist_p;
  return;
}


void makeAsymmHist_ThirdJet(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Float_t histLow, Float_t histHi, Int_t centLow, Int_t centHi, Bool_t montecarlo = false)
{
  inFile_p->cd();

  const char* title = Form("%sAsymm_3_%d%d_%s", algType[setNum], (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag1);

  TH1F* asymmHist3_p;

  TString name = Form("%s_h(%d, %f, %f)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, 5*TMath::Pi()/6);
  TCut isThirdJet = Form("AlgThirdJtPt[%d] > 50", setNum);

  TCut etaCut = makeEtaCut(setNum, 0.5);
  TCut thirdEtaCut = Form("TMath::Abs(AlgThirdJtEta[%d]) < 0.5", setNum);

  const char* var =  Form("(AlgLeadJtPt[%d] + (AlgSubLeadJtPt[%d]*cos(AlgJtDelPhi[%d]) + AlgThirdJtPt[%d]*cos(getDPHI(AlgThirdJtPhi[%d], AlgLeadJtPhi[%d]))))/(AlgLeadJtPt[%d] - (AlgSubLeadJtPt[%d]*cos(AlgJtDelPhi[%d]) + AlgThirdJtPt[%d]*cos(getDPHI(AlgThirdJtPhi[%d], AlgLeadJtPhi[%d]))))", setNum, setNum, setNum, setNum, setNum, setNum, setNum, setNum, setNum, setNum, setNum, setNum);

  if(montecarlo)
    getTree_p->Project(name, var, Form("setWeight[%d]", setNum)*(centCut && setCut && etaCut && thirdEtaCut  && delPhiCut && isThirdJet));
  else 
    getTree_p->Project(name, var, centCut && setCut && etaCut && thirdEtaCut && delPhiCut && isThirdJet);

  asymmHist3_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  outFile_p = new TFile(outName, "UPDATE");

  asymmHist3_p->Sumw2();

  niceTH1(asymmHist3_p, .55, 0., 405, 506);

  asymmHist3_p->SetYTitle("Event Fraction");
  asymmHist3_p->SetXTitle("A_{3} = (p_{T,1} + p_{T,2}^{||} + p_{T,3}^{||})/(p_{T,1} - p_{T,2}^{||} - p_{T,3}^{||})");

  asymmHist3_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
  return;
}



void makeRatioDelPhiHist_ThirdJet(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t centLow, Int_t centHi, Bool_t montecarlo = false)
{
  inFile_p->cd();

  const char* title12 = Form("%sRatioDelPhi_12_%d%d_%s", algType[setNum], (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag1);
  const char* title13 = Form("%sRatioDelPhi_13_%d%d_%s", algType[setNum], (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag1);

  std::cout << "Titles: " << title12 << ", " << title13 << std::endl;

  TH1F* getValHist12_p[6];
  TH1F* getValHist13_p[6];

  TH1F* ratioDelPhiHist12_p = new TH1F(Form("%s_h", title12), Form("%s_h", title12), 6, 2*TMath::Pi()/3, TMath::Pi());
  TH1F* ratioDelPhiHist13_p = new TH1F(Form("%s_h", title13), Form("%s_h", title13), 6, 2*TMath::Pi()/3, TMath::Pi());

  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut isThirdJet = Form("AlgThirdJtPt[%d] > 50", setNum);
  TCut etaCut = makeEtaCut(setNum, 0.5);
  TCut thirdEtaCut = Form("TMath::Abs(AlgThirdJtEta[%d]) < 0.5", setNum);

  const char* var12 = Form("AlgSubLeadJtPt[%d]/AlgLeadJtPt[%d]", setNum, setNum);
  const char* var13 = Form("AlgThirdJtPt[%d]/AlgLeadJtPt[%d]", setNum, setNum);

  for(Int_t binIter = 0; binIter < 6; binIter++){
    TCut delPhiCut12 = Form("AlgJtDelPhi[%d] > %f && AlgJtDelPhi[%d] < %f", setNum, 2*TMath::Pi()/3 + binIter*TMath::Pi()/18, setNum, 2*TMath::Pi()/3 + (binIter + 1)*TMath::Pi()/18);
    TCut delPhiCut13 = Form("TMath::Abs(getDPHI(AlgLeadJtPhi[%d], AlgThirdJtPhi[%d])) > %f && TMath::Abs(getDPHI(AlgLeadJtPhi[%d], AlgThirdJtPhi[%d])) < %f", setNum, setNum, 2*TMath::Pi()/3 + binIter*TMath::Pi()/18, setNum, setNum, 2*TMath::Pi()/3 + (binIter + 1)*TMath::Pi()/18);
    TString name12 = Form("%d12_h(10, 0, 1)", binIter);
    TString name13 = Form("%d13_h(10, 0, 1)", binIter);

    if(montecarlo){
      getTree_p->Project(name12, var12, Form("setWeight[%d]", setNum)*(centCut && setCut && etaCut && thirdEtaCut && isThirdJet && delPhiCut12));
      getTree_p->Project(name13, var13, Form("setWeight[%d]", setNum)*(centCut && setCut && etaCut && thirdEtaCut && isThirdJet && delPhiCut13));
    }
    else{
      getTree_p->Project(name12, var12, centCut && setCut && etaCut && thirdEtaCut && isThirdJet && delPhiCut12);
      getTree_p->Project(name13, var13, centCut && setCut && etaCut && thirdEtaCut && isThirdJet && delPhiCut13);
    }

    getValHist12_p[binIter] = (TH1F*)inFile_p->Get(Form("%d12_h", binIter));
    getValHist13_p[binIter] = (TH1F*)inFile_p->Get(Form("%d13_h", binIter));

    ratioDelPhiHist12_p->SetBinContent(binIter + 1, getValHist12_p[binIter]->GetMean());
    ratioDelPhiHist12_p->SetBinError(binIter + 1, getValHist12_p[binIter]->GetMeanError());

    ratioDelPhiHist13_p->SetBinContent(binIter + 1, getValHist13_p[binIter]->GetMean());
    ratioDelPhiHist13_p->SetBinError(binIter + 1, getValHist13_p[binIter]->GetMeanError());
  }

  outFile_p = new TFile(outName, "UPDATE");

  ratioDelPhiHist12_p->Sumw2();
  ratioDelPhiHist13_p->Sumw2();

  niceTH1(ratioDelPhiHist12_p, .85, .40, 405, 506, false);
  niceTH1(ratioDelPhiHist13_p, .65, .20, 405, 506, false);

  ratioDelPhiHist12_p->SetYTitle("<p_{T,2}/p_{T,1}>");
  ratioDelPhiHist12_p->SetXTitle("#Delta #phi_{1,2}");

  ratioDelPhiHist13_p->SetYTitle("<p_{T,3}/p_{T,1}>");
  ratioDelPhiHist13_p->SetXTitle("#Delta #phi_{1,3}");

  ratioDelPhiHist12_p->Write();
  ratioDelPhiHist13_p->Write();
  outFile_p->Close();

  delete outFile_p;
  return;
}


void makeDelPhiHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Float_t histLow, Float_t histHi, Int_t centLow, Int_t centHi, Bool_t montecarlo = false)
{
  inFile_p->cd();

  const char* title = Form("%sDelPhi_%d%d_%s", algType[setNum], (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag1);

  TH1F* delPhiHist_p;

  TString name = Form("%s_h(%d, %f, %f)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut lJtCut = Form("AlgLeadJtPt[%d] > 120", setNum);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, setDelPhiCut);

  if(montecarlo)
    getTree_p->Project(name, Form("AlgJtDelPhi[%d]", setNum), Form("setWeight[%d]", setNum)*(centCut && setCut && lJtCut && thirdJtVeto && delPhiCut));
  else
    getTree_p->Project(name, Form("AlgJtDelPhi[%d]", setNum), centCut && setCut && lJtCut && thirdJtVeto && delPhiCut);

  delPhiHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  delPhiHist_p->Sumw2();
  if(setDelPhiCut < .5)
    niceTH1(delPhiHist_p, 1., .00001, 405, 506);
  else
    niceTH1(delPhiHist_p, 1., .001, 405, 506);

  delPhiHist_p->SetYTitle("Event Fraction");
  delPhiHist_p->SetXTitle("#Delta #phi_{1,2}");

  outFile_p = new TFile(outName, "UPDATE");
  delPhiHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
  return;
}


void makeJtPtHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Int_t histLow, Int_t histHi, Int_t centLow, Int_t centHi, const char* Sub_Lead = "", Bool_t montecarlo = false)
{
  inFile_p->cd();

  const char* title = Form("%s%sLeadJtPt_%d%d_%s", algType[setNum], Sub_Lead, (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag1);

  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, setDelPhiCut);

  if(montecarlo)
    getTree_p->Project(name, Form("Alg%sLeadJtPt[%d]", Sub_Lead, setNum), Form("setWeight[%d]", setNum)*(centCut && setCut && thirdJtVeto && delPhiCut));
  else
    getTree_p->Project(name, Form("Alg%sLeadJtPt[%d]", Sub_Lead, setNum), centCut && setCut && thirdJtVeto && delPhiCut);

  asymmHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  asymmHist_p->Sumw2();
  niceTH1(asymmHist_p, 1., .0001, 405, 506);

  asymmHist_p->SetYTitle("Event Fraction");

  if(strcmp("", Sub_Lead) == 0)
    asymmHist_p->SetXTitle("p_{T,1} (GeV/c)");
  else if(strcmp("Sub", Sub_Lead) == 0)
    asymmHist_p->SetXTitle("p_{T,2} (GeV/c)");

  outFile_p = new TFile(outName, "UPDATE");
  asymmHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
  return;
}


void makeJtEtaHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Float_t histLow, Float_t histHi, Int_t centLow, Int_t centHi, const char* Sub_Lead = "", Bool_t montecarlo = false)
{
  inFile_p->cd();

  const char* title = Form("%s%sLeadJtEta_%d%d_%s", algType[setNum], Sub_Lead, (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag1);

  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %f, %f)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, setDelPhiCut);

  if(montecarlo)
    getTree_p->Project(name, Form("Alg%sLeadJtEta[%d]", Sub_Lead, setNum), Form("setWeight[%d]", setNum)*(centCut && setCut && thirdJtVeto && delPhiCut));
  else
    getTree_p->Project(name, Form("Alg%sLeadJtEta[%d]", Sub_Lead, setNum), centCut && setCut && thirdJtVeto && delPhiCut);
  
  asymmHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  asymmHist_p->Sumw2();
  niceTH1(asymmHist_p, .25, 0., 405, 506);

  asymmHist_p->SetYTitle("Event Fraction");

  if(strcmp("", Sub_Lead) == 0)
    asymmHist_p->SetXTitle("#eta_{1}");
  else if(strcmp("Sub", Sub_Lead) == 0)
    asymmHist_p->SetXTitle("#eta_{2}");

  outFile_p = new TFile(outName, "UPDATE");
  asymmHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
  return;
}


void addJtHistToPanel(TFile* file_p, TCanvas* canv_p, Int_t setNum, const char* jtVarIn, Int_t centLow, Int_t centHi, Int_t pos, Bool_t logY = false, const char* subOrLead = "")
{
  file_p->cd();

  std::cout << Form("%s%s%s_%d%d_%s_h", algType[setNum], subOrLead, jtVarIn, centLow, centHi, fileTag1) << std::endl;

  TH1F* hist_p = (TH1F*)file_p->Get(Form("%s%s%s_%d%d_%s_h", algType[setNum], subOrLead, jtVarIn, centLow, centHi, fileTag1));
  canv_p->cd(pos);

  if(logY)
    gPad->SetLogy();

  if(pos == 2 || pos == 3 || pos == 4)
    hist_p->SetYTitle("");

  hist_p->SetFillColor(17);
  hist_p->SetMarkerStyle(6);
  hist_p->SetMarkerSize(.5);

  hist_p->DrawCopy("HIST");
  hist_p->DrawCopy("E1 SAME");

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  if(strcmp("Asymm", jtVarIn) == 0 || strcmp("LeadJtEta", jtVarIn) == 0  || strcmp("DelPhi", jtVarIn) == 0 || strcmp("Asymm_3", jtVarIn) == 0)
    label_p->DrawLatex(.7, .3, Form("%d-%d%%", centLow, centHi));
  else
    label_p->DrawLatex(.3, .3, Form("%d-%d%%", centLow, centHi));

  
  if(pos == 1){
    label_p->DrawLatex(.20, .94, "CMS Preliminary");
  }
  else if(pos == 2){
    if(strcmp("LeadJtPt", jtVarIn) != 0){
      if(strcmp(jtVarIn, "Asymm") == 0){
	label_p->DrawLatex(.12, .875, Form("Anti-k_{T} %s, R = 0.3", algType[setNum]));
	label_p->DrawLatex(.12, .80, "p_{T,1} > 120, p_{T,2} > 50 GeV/c");
      }
      else if(strcmp(jtVarIn, "AsymmDelPhi") == 0 || strcmp(jtVarIn, "RatioDelPhi") == 0){
	label_p->DrawLatex(.12, .875, Form("Anti-k_{T} %s, R = 0.3", algType[setNum]));
	label_p->DrawLatex(.12, .80, "p_{T,1} > 120, p_{T,2} > 50 GeV/c");
      }
      else if(strcmp(jtVarIn, "RatioDelPhi_12") == 0 || strcmp(jtVarIn, "RatioDelPhi_13") == 0){
	label_p->DrawLatex(.12, .875, Form("Anti-k_{T} %s, R = 0.3", algType[setNum]));
	label_p->DrawLatex(.12, .80, "p_{T,1} > 120, p_{T,2(3)} > 50 GeV/c");
      }
      else if(strcmp(jtVarIn, "Asymm_3") == 0){
	label_p->DrawLatex(.12, .875, Form("Anti-k_{T} %s, R = 0.3", algType[setNum]));
	label_p->DrawLatex(.12, .80, "p_{T,1} > 120, p_{T,2} > 50 GeV/c");
      }
      else{
	label_p->DrawLatex(.15, .875, Form("Anti-k_{T} %s, R = 0.3", algType[setNum]));
	label_p->DrawLatex(.15, .80, "p_{T,1} > 120, p_{T,2} > 50 GeV/c");
	label_p->DrawLatex(.15, .725, "|#eta|_{1,2} < 2.0");
	//      label_p->DrawLatex(.15, .65, "#Delta #phi_{1,2} > 2#pi/3");
      }
    }
    else{
      label_p->DrawLatex(.42, .875, Form("Anti-k_{T} %s, R = 0.3", algType[setNum]));
      label_p->DrawLatex(.42, .80, "p_{T,1} > 120, p_{T,2} > 50 GeV/c");
      label_p->DrawLatex(.42, .725, "|#eta|_{1,2} < 2.0");
      //      label_p->DrawLatex(.45, .65, "#Delta #phi_{1,2} > 2#pi/3");
    }
  }
  else if(pos == 3){
    if(strcmp(jtVarIn, "Asymm") == 0){
      label_p->DrawLatex(.12, .875, "|#eta|_{1,2} < 0.5");
      label_p->DrawLatex(.12, .80, "#Delta #phi_{1,2} > 5#pi/6");
    }
    else if(strcmp(jtVarIn, "Asymm_3") == 0){
      label_p->DrawLatex(.12, .875, "|#eta|_{1,2,3} < 0.5");
      label_p->DrawLatex(.12, .80, "#Delta #phi_{1,2} > 5#pi/6");
    }
    else if(strcmp(jtVarIn, "AsymmDelPhi") == 0 || strcmp(jtVarIn, "RatioDelPhi") == 0){
      label_p->DrawLatex(.12, .875, "|#eta|_{1,2} < 0.5");
    }
    else if(strcmp(jtVarIn, "RatioDelPhi_12") == 0 || strcmp(jtVarIn, "RatioDelPhi_13") == 0){
      label_p->DrawLatex(.12, .875, "|#eta|_{1,2,3} < 0.5");
    }
  }
  else if(pos == 4){
    if(strcmp(jtVarIn, "Asymm_3") == 0)
      label_p->DrawLatex(.12, .875, "Require Third Jet w/ p_{T} > 50 GeV/c");

    if(strcmp(jtVarIn, "LeadJtPt") == 0)
      label_p->DrawLatex(.12, .875, "#sqrt{s_{NN}} = 2.76 TeV, 150/#mub");
    else
      label_p->DrawLatex(.12, .80, "#sqrt{s_{NN}} = 2.76 TeV, 150/#mub");
  }

  delete label_p;
  return;
}


void addJtSubHistToPanel(TFile* file_p, TCanvas* canv_p, Int_t setNum, const char* jtVarIn, Int_t centLow, Int_t centHi, Int_t pos, Bool_t logY = false)
{
  file_p->cd();

  TH1F* hist1_p = (TH1F*)file_p->Get(Form("%s%s_%d%d_%s_h", algType[setNum], jtVarIn, centLow, centHi, fileTag1));
  TH1F* hist2_p = (TH1F*)file_p->Get(Form("%sSub%s_%d%d_%s_h", algType[setNum], jtVarIn, centLow, centHi, fileTag1));

  hist1_p->Sumw2();
  hist2_p->Sumw2();

  hist1_p->Add(hist2_p, -1);

  canv_p->cd(pos);

  if(logY)
    gPad->SetLogy();

  if(pos == 1)
    hist1_p->SetYTitle("Fraction Lead - Fraction SubLead");
  else if(pos == 2 || pos == 3 || pos == 4)
    hist1_p->SetYTitle("");

  hist1_p->SetFillColor(17);
  hist1_p->SetMarkerStyle(6);
  hist1_p->SetMarkerSize(.5);

  hist1_p->DrawCopy("HIST");
  hist1_p->DrawCopy("Same E1");

  TLine* zeroLine_p = new TLine(0., 0., 0.5, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw("SAME");

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  if(strcmp("Asymm", jtVarIn) == 0 || strcmp("LeadJtEta", jtVarIn) == 0  || strcmp("DelPhi", jtVarIn) == 0 || strcmp("Asymm_3", jtVarIn) == 0)
    label_p->DrawLatex(.7, .3, Form("%d-%d%%", centLow, centHi));
  else
    label_p->DrawLatex(.3, .3, Form("%d-%d%%", centLow, centHi));

  if(pos == 2){
    if(strcmp("LeadJtPt", jtVarIn) != 0){
      label_p->DrawLatex(.15, .875, Form("Anti-k_{T} %s, R = 0.3", algType[setNum]));
      label_p->DrawLatex(.15, .80, "p_{T,1} > 120, p_{T,2} > 50 GeV/c");
      label_p->DrawLatex(.15, .725, "|#eta|_{1,2} < 0.5");
      //      label_p->DrawLatex(.15, .65, "#Delta #phi_{1,2} > 2#pi/3");
    }
    else{
      label_p->DrawLatex(.45, .875, Form("Anti-k_{T} %s, R = 0.3", algType[setNum]));
      label_p->DrawLatex(.45, .80, "p_{T,1} > 120, p_{T,2} > 50 GeV/c");
      label_p->DrawLatex(.45, .725, "|#eta|_{1,2} < 0.5");
      //      label_p->DrawLatex(.45, .65, "#Delta #phi_{1,2} > 2#pi/3");
    }
  }    

  delete label_p;
  return;
}


void addJtHistToPanel_OVER(TFile* file_p, const char* fileTag2, TCanvas* canv_p, Int_t setNum, const char* jtVarIn, Int_t centLow, Int_t centHi, Int_t pos, Bool_t logY = false, const char* subOrLead = "")
{
  file_p->cd();

  std::cout << Form("%s%s%s_%d%d_%s_h", algType[setNum], subOrLead, jtVarIn, centLow, centHi, fileTag2) << std::endl;;

  TH1F* hist_p = (TH1F*)file_p->Get(Form("%s%s%s_%d%d_%s_h", algType[setNum], subOrLead, jtVarIn, centLow, centHi, fileTag2));
  canv_p->cd(pos);

  if(logY)
    gPad->SetLogy();

  if(pos == 2 || pos == 3 || pos == 4)
    hist_p->SetYTitle("");

  hist_p->SetMarkerColor(kRed);
  hist_p->SetLineColor(kRed);

  hist_p->DrawCopy("E1 SAME");
  return;
}


void addJtHistToPanel_pp(TFile* file_p, const char* fileTag2, TCanvas* canv_p, Int_t setNum, const char* jtVarIn, Bool_t logY = false, const char* subOrLead = "")
{
  file_p->cd();

  std::cout << Form("%s%s%s_pp_%s_h", algType_pp[setNum], subOrLead, jtVarIn, fileTag2) << std::endl;

  TH1F* hist_p = (TH1F*)file_p->Get(Form("%s%s%s_pp_%s_h", algType_pp[setNum], subOrLead, jtVarIn, fileTag2));
  canv_p->cd(1);

  if(logY)
    gPad->SetLogy();

  hist_p->SetMarkerColor(kBlue);
  hist_p->SetLineColor(kBlue);
  hist_p->DrawCopy("E1 SAME");

  canv_p->cd(2);
  hist_p->DrawCopy("E1 SAME");

  canv_p->cd(3);
  hist_p->DrawCopy("E1 SAME");

  canv_p->cd(4);
  hist_p->DrawCopy("E1 SAME");

  canv_p->cd(1);

  return;
}


void makeJtVarPanel(const char* fileName1, Int_t setNum, const char* jtVarIn, Bool_t logY = false, const char* subOrLead = "", const char* fileName2 = "", const char* fileTag2 = "", const char* ppName = "", const char* ppFileTag = "")
{
  TCanvas* jtVarPanel_p = new TCanvas(Form("%s%s%sPanel_%s", algType[setNum], subOrLead, jtVarIn, fileTag1), Form("%s%s%sPanel_%s", algType[setNum], subOrLead, jtVarIn, fileTag1), 900, 300);

  jtVarPanel_p->Divide(4,1,0,0);

  TFile* panelFile_p = new TFile(fileName1, "UPDATE");

  std::cout << "1" << std::endl;

  addJtHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 50, 100, 1, logY, subOrLead);
  addJtHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 30, 50, 2, logY, subOrLead);
  addJtHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 10, 30, 3, logY, subOrLead);
  addJtHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 0, 10, 4, logY, subOrLead);

  std::cout << "2" << std::endl;

  if(strcmp(fileName2, "") != 0){
    TFile* overFile_p = new TFile(fileName2, "UPDATE");
    TFile* ppFile_p;

    std::cout << "3" << std::endl;

    addJtHistToPanel_OVER(overFile_p, fileTag2, jtVarPanel_p, setNum, jtVarIn, 50, 100, 1, logY, subOrLead);
    addJtHistToPanel_OVER(overFile_p, fileTag2, jtVarPanel_p, setNum, jtVarIn, 30, 50, 2, logY, subOrLead);
    addJtHistToPanel_OVER(overFile_p, fileTag2, jtVarPanel_p, setNum, jtVarIn, 10, 30, 3, logY, subOrLead);
    addJtHistToPanel_OVER(overFile_p, fileTag2, jtVarPanel_p, setNum, jtVarIn, 0, 10, 4, logY, subOrLead);

    std::cout << "4" << std::endl;

    TLegend* leg;
    if(strcmp(jtVarIn, "LeadJtPt") != 0)
      leg = new TLegend(0.20, 0.75, 0.40, 0.90);
    else
      leg = new TLegend(0.45, 0.75, 0.75, 0.90);

    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(.06);
    leg->SetBorderSize(0);

    TH1F* legHist_p = new TH1F("legHist", "legHist", 10, 0., 1.);
    legHist_p->SetMarkerColor(kRed);
    leg->AddEntry(legHist_p, "PbPb", "P");

    TH1F* legHist2_p = new TH1F("legHist2", "legHist2", 10, 0., 1.);
    legHist2_p->SetFillColor(17);
    leg->AddEntry(legHist2_p, "PYTHIA+HYDJET", "F");

    jtVarPanel_p->cd(1);

    if(strcmp(ppName, "") != 0){
      ppFile_p = new TFile(ppName, "READ");

      std::cout << "5" << std::endl;

      addJtHistToPanel_pp(ppFile_p, ppFileTag, jtVarPanel_p, setNum, jtVarIn, logY, subOrLead);

      std::cout << "6" << std::endl;

      TH1F* legHist3_p = new TH1F("legHist3", "legHist3", 10, 0., 1.);
      legHist3_p->SetMarkerColor(kBlue);
      leg->AddEntry(legHist3_p, "PP", "P");
    }

    leg->Draw("SAME");    

    claverCanvasSaving(jtVarPanel_p, Form("../pdfDir/%s%s%sPanelOver_%s_c", algType[setNum], subOrLead, jtVarIn, fileTag1), "pdf");
    overFile_p->cd();
    jtVarPanel_p->Write(Form("%s%s%sPanelOver_%s", algType[setNum], subOrLead, jtVarIn, fileTag1));

    TFile* sepFile_p = new TFile("jtVarDataOnMC.root", "RECREATE");
    jtVarPanel_p->Write(Form("%s%s%sPanelOver_%s_%s", algType[setNum], subOrLead, jtVarIn, fileTag1, fileTag2));
    sepFile_p->Close();

    delete sepFile_p;

    delete legHist2_p;
    delete legHist_p;
    delete leg;

    overFile_p->Close();
    delete overFile_p;
  }
  else{
    claverCanvasSaving(jtVarPanel_p, Form("../pdfDir/%s%s%sPanel_%s_c", algType[setNum], subOrLead, jtVarIn, fileTag1), "pdf");
    jtVarPanel_p->Write();
  }


  panelFile_p->Close();
  delete panelFile_p;
  delete jtVarPanel_p;
  return;
}



void makeJtVarSubPanel(const char* fileName1, Int_t setNum, const char* jtVarIn, Bool_t logY = false, const char* fileName2 = "", const char* fileTag2 = "", const char* ppName = "", const char* ppFileTag = "")
{
  TCanvas* jtVarPanel_p = new TCanvas(Form("%s%sMinPanel_%s", algType[setNum], jtVarIn, fileTag1), Form("%s%sMinPanel_%s", algType[setNum], jtVarIn, fileTag1), 900, 300);

  jtVarPanel_p->Divide(4,1,0,0);

  TFile* panelFile_p = new TFile(fileName1, "UPDATE");

  addJtSubHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 50, 100, 1, logY);
  addJtSubHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 30, 50, 2, logY);
  addJtSubHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 10, 30, 3, logY);
  addJtSubHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 0, 10, 4, logY);

  claverCanvasSaving(jtVarPanel_p, Form("../pdfDir/%s%sMinPanel_%s_c", algType[setNum], jtVarIn, fileTag1), "pdf");
  jtVarPanel_p->Write();
  
  panelFile_p->Close();
  delete panelFile_p;
  delete jtVarPanel_p;
  return;
}


void cfmDiJet_JtVarPlots(const char* inName, const char* outName, Bool_t montecarlo = false, const char* inName2 = "", const char* fileTag2 = "", const char* ppName = "", const char* ppFileTag = "", Bool_t Veto = false)
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
  std::cout << "FileTag2 is: " << fileTag2 << std::endl;

  inFile1_p = new TFile(inName, "READ");
  inTree_p = (TTree*)inFile1_p->Get("jetTree");
  inTree_p->AddFriend("trackTree");

  Int_t jetAlgMax = 4;

  if(montecarlo){
    inTree_p->AddFriend("genTree");
    jetAlgMax = 5;
  }


  Int_t centLow[5] = {0, 0, 20, 60, 100};
  Int_t centHi[5] = {199, 19, 59, 99, 199};

  const char* subOrLead[2] = {"", "Sub"};
  Int_t leadOrSubBins[6] = {18, 120, 300, 25, 50, 300};

  const char* jtVar[10] = {"Asymm", "Asymm_3", "AsymmDelPhi", "RatioDelPhi", "RatioDelPhi_12", "RatioDelPhi_13", "Ratio", "DelPhi", "LeadJtPt", "LeadJtEta"};
  Bool_t isLogY[10] = {false, false, false, false, false, false, false, true, true, false};

  for(Int_t algIter = 0; algIter < jetAlgMax; algIter++){

    makeCentHist(inFile1_p, inTree_p, "centHist_eventSelect_120_5pi6.root", algIter, fileTag2);

    if(Veto)
      thirdJtVeto = Form("(AlgThirdJtPt[%d] > 50 || AlgThirdJtPt[%d] < 0)", algIter, algIter);

    std::cout << thirdJtVeto << std::endl;

    for(Int_t centIter = 0; centIter < 5; centIter++){
      makeAsymmHist(inFile1_p, inTree_p, outName, algIter, 10, 0, 1, centLow[centIter], centHi[centIter], montecarlo);
      makePtRatioHist(inFile1_p, inTree_p, outName, algIter, 10, 0, 1, centLow[centIter], centHi[centIter], montecarlo);
      makeAsymmHist_ThirdJet(inFile1_p, inTree_p, outName, algIter, 15, -0.5, 1., centLow[centIter], centHi[centIter], montecarlo);
      makeDelPhiHist(inFile1_p, inTree_p, outName, algIter, 20, 0, TMath::Pi(), centLow[centIter], centHi[centIter]);
      //      makeDelPhiHist(inFile1_p, inTree_p, outName, algIter, 20, 2*TMath::Pi()/3, TMath::Pi(), centLow[centIter], centHi[centIter], montecarlo);
      makeAsymmDelPhiHist(inFile1_p, inTree_p, outName, algIter, centLow[centIter], centHi[centIter], montecarlo);
      makeRatioDelPhiHist(inFile1_p, inTree_p, outName, algIter, centLow[centIter], centHi[centIter], montecarlo);
      makeRatioDelPhiHist_ThirdJet(inFile1_p, inTree_p, outName, algIter, centLow[centIter], centHi[centIter], montecarlo);

      for(Int_t subIter = 0; subIter < 2; subIter++){
	makeJtPtHist(inFile1_p, inTree_p, outName, algIter, leadOrSubBins[subIter*3], leadOrSubBins[subIter*3 + 1], leadOrSubBins[subIter*3 + 2], centLow[centIter], centHi[centIter], subOrLead[subIter], montecarlo);
	makeJtEtaHist(inFile1_p, inTree_p, outName, algIter, 10, -2.0, 2.0, centLow[centIter], centHi[centIter], subOrLead[subIter], montecarlo);
      }      
    }

    for(Int_t varIter = 0; varIter < 10; varIter++){

      makeJtVarPanel(outName, algIter, jtVar[varIter], isLogY[varIter], "");

      if(varIter > 7)
	makeJtVarPanel(outName, algIter, jtVar[varIter], isLogY[varIter], "Sub");

      if(strcmp("", fileTag2) != 0 && algIter != 4){
	makeJtVarPanel(outName, algIter, jtVar[varIter], isLogY[varIter], "", inName2, fileTag2, ppName, ppFileTag);

	if(varIter > 7)
	  makeJtVarPanel(outName, algIter, jtVar[varIter], isLogY[varIter], "Sub", inName2, fileTag2, ppName, ppFileTag);
      }
    }

    makeJtVarSubPanel(outName, algIter, "LeadJtEta", false);

  }

  inFile1_p->Close();
  delete inFile1_p;
  return;
}
