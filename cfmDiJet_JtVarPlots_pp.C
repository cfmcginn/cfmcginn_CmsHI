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

const char* algType_pp[5] = {"PuPF", "PuCalo", "PF", "Calo", "T"};
TCut thirdJtVeto = "";
Float_t setDelPhiCut = 5*TMath::Pi()/6;

//append to every histo so know which sample via shorthand on workblog                                                          

const char* fileTag1;

//shorthands, w/ _CFMSKIM.h                                                                                                                      

const char* PPDataA = "PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82_CFMSKIM.root";

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

// 0 == PuPF, 1 == PuCalo, 2 == PF, 3 == Calo, 4 == Truth

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


void makeAsymmHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Int_t histLow, Int_t histHi)
{
  inFile_p->cd();

  const char* title = Form("%sAsymm_pp_%s", algType_pp[setNum], fileTag1);

  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, 5*TMath::Pi()/6);
  TCut etaCut = makeEtaCut(setNum, 0.5);

  getTree_p->Project(name, Form("AlgJtAsymm[%d]", setNum), setCut && thirdJtVeto && delPhiCut && etaCut);
  asymmHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  asymmHist_p->Sumw2();
  niceTH1(asymmHist_p, .55, 0., 405, 506);

  asymmHist_p->SetYTitle("Event Fraction");
  asymmHist_p->SetXTitle("A_{J} = (p_{T,1} - p_{T,2})/(p_{T,1} + p_{T,2})");

  outFile_p = new TFile(outName, "UPDATE");
  asymmHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
  return;
}


void makePtRatioHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Int_t histLow, Int_t histHi)
{
  inFile_p->cd();

  const char* title = Form("%sRatio_pp_%s", algType_pp[setNum], fileTag1);

  TH1F* ratioHist_p;

  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, setDelPhiCut);
  TCut etaCut = makeEtaCut(setNum, 0.5);

  getTree_p->Project(name, Form("AlgSubLeadJtPt[%d]/AlgLeadJtPt[%d]", setNum, setNum), setCut && thirdJtVeto && delPhiCut && etaCut);

  ratioHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  outFile_p = new TFile(outName, "UPDATE");

  ratioHist_p->Sumw2();

  niceTH1(ratioHist_p, .25, 0., 405, 506);

  ratioHist_p->SetYTitle("Event Fraction");
  ratioHist_p->SetXTitle("p_{T,2})/p_{T,1}");

  ratioHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
  return;
}


void makeAsymmDelPhiHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum)
{
  inFile1_p->cd();

  const char* title = Form("%sAsymmDelPhi_pp_%s", algType_pp[setNum], fileTag1);
  TH1F* getValHist_p[4];

  TH1F* asymmDelPhiHist_p = new TH1F(Form("%s_h", title), Form("%s_h", title), 6, 2*TMath::Pi()/3, TMath::Pi());

  TCut setCut = makeSetCut(setNum);
  TCut etaCut = makeEtaCut(setNum, 0.5);

  for(Int_t binIter = 0; binIter < 6; binIter++){
    TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f && AlgJtDelPhi[%d] < %f", setNum, 2*TMath::Pi()/3 + binIter*TMath::Pi()/18, setNum, 2*TMath::Pi()/3 + (binIter + 1)*TMath::Pi()/18);
    TString name = Form("%d_h(10, 0, 1)", binIter);

    getTree_p->Project(name, Form("AlgJtAsymm[%d]", setNum), setCut && delPhiCut && etaCut);

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



void makeRatioDelPhiHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum)
{
  inFile1_p->cd();

  const char* title = Form("%sRatioDelPhi_pp_%s", algType_pp[setNum], fileTag1);
  TH1F* getValHist_p[4];

  TH1F* ratioDelPhiHist_p = new TH1F(Form("%s_h", title), Form("%s_h", title), 6, 2*TMath::Pi()/3, TMath::Pi());

  TCut setCut = makeSetCut(setNum);
  TCut etaCut = makeEtaCut(setNum, 0.5);

  for(Int_t binIter = 0; binIter < 6; binIter++){
    TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f && AlgJtDelPhi[%d] < %f", setNum, 2*TMath::Pi()/3 + binIter*TMath::Pi()/18, setNum, 2*TMath::Pi()/3 + (binIter + 1)*TMath::Pi()/18);
    TString name = Form("%d_h(10, 0, 1)", binIter);

    getTree_p->Project(name, Form("AlgSubLeadJtPt[%d]/AlgLeadJtPt[%d]", setNum, setNum), setCut && delPhiCut && etaCut);

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



void makeAsymmHist_ThirdJet(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Float_t histLow, Float_t histHi)
{
  inFile_p->cd();

  const char* title = Form("%sAsymm_3_pp_%s", algType_pp[setNum], fileTag1);

  TH1F* asymmHist3_p;

  TString name = Form("%s_h(%d, %f, %f)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, 5*TMath::Pi()/6);
  //  TCut thirdPhiCut = Form("getDPHI(AlgLeadJtPhi[%d], AlgThirdJtPhi[%d]) > %f", setNum, setNum, TMath::Pi()/2);                                                 
  TCut isThirdJet = Form("AlgThirdJtPt[%d] > 50", setNum);
  //  TCut noFourthJet = Form("AlgFourthJtPt[%d] < 30", setNum);                                                                                                   

  TCut etaCut = makeEtaCut(setNum, 0.5);
  TCut thirdEtaCut = Form("TMath::Abs(AlgThirdJtEta[%d]) < 0.5", setNum);

  //  std::cout << delPhiCut << std::endl;                                                                                                                         

  const char* var =  Form("(AlgLeadJtPt[%d] + (AlgSubLeadJtPt[%d]*cos(AlgJtDelPhi[%d]) + AlgThirdJtPt[%d]*cos(getDPHI(AlgThirdJtPhi[%d], AlgLeadJtPhi[%d]))))/(AlgLeadJtPt[%d] - (AlgSubLeadJtPt[%d]*cos(AlgJtDelPhi[%d]) + AlgThirdJtPt[%d]*cos(getDPHI(AlgThirdJtPhi[%d], AlgLeadJtPhi[%d]))))", setNum, setNum, setNum, setNum, setNum, setNum, setNum, setNum, setNum, setNum, setNum, setNum);

  getTree_p->Project(name, var, setCut && etaCut && thirdEtaCut && delPhiCut && isThirdJet/* && thirdPhiCut && noFourthJet*/);

  std::cout << setCut << ", " << etaCut << ", " <</* delPhiCut <<*/ ", " << isThirdJet << ", "/* << thirdPhiCut*/ << std::endl;

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


void makeRatioDelPhiHist_ThirdJet(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum)
{
  inFile_p->cd();

  const char* title12 = Form("%sRatioDelPhi_12_pp_%s", algType_pp[setNum], fileTag1);
  const char* title13 = Form("%sRatioDelPhi_13_pp_%s", algType_pp[setNum], fileTag1);
  TH1F* getValHist12_p[4];
  TH1F* getValHist13_p[4];

  TH1F* ratioDelPhiHist12_p = new TH1F(Form("%s_h", title12), Form("%s_h", title12), 6, 2*TMath::Pi()/3, TMath::Pi());
  TH1F* ratioDelPhiHist13_p = new TH1F(Form("%s_h", title13), Form("%s_h", title13), 6, 2*TMath::Pi()/3, TMath::Pi());

  TCut setCut = makeSetCut(setNum);
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

    getTree_p->Project(name12, var12, setCut && etaCut && thirdEtaCut && isThirdJet && delPhiCut12);
    getTree_p->Project(name13, var13, setCut && etaCut && thirdEtaCut && isThirdJet && delPhiCut13);
    
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

  ratioDelPhiHist12_p->Write(Form("%s_h", title12));
  ratioDelPhiHist13_p->Write(Form("%s_h", title13));

  outFile_p->Close();

  delete outFile_p;
  return;
}



void makeDelPhiHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Float_t histLow, Float_t histHi)
{
  inFile_p->cd();

  const char* title = Form("%sDelPhi_pp_%s", algType_pp[setNum], fileTag1);

  TH1F* delPhiHist_p;

  TString name = Form("%s_h(%d, %f, %f)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut lJtCut = Form("AlgLeadJtPt[%d] > 120", setNum);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, setDelPhiCut);

  getTree_p->Project(name, Form("AlgJtDelPhi[%d]", setNum), setCut && lJtCut && thirdJtVeto && delPhiCut);
  delPhiHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  delPhiHist_p->Sumw2();
  if(setDelPhiCut < .5)
    niceTH1(delPhiHist_p, 1., .00001, 405, 506);
  else
    niceTH1(delPhiHist_p, .2, .0, 405, 506);

  delPhiHist_p->SetYTitle("Event Fraction");
  delPhiHist_p->SetXTitle("#Delta #phi_{1,2}");

  outFile_p = new TFile(outName, "UPDATE");
  delPhiHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
  return;
}


void makeJtPtHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Int_t histLow, Int_t histHi, const char* Sub_Lead = "")
{
  inFile_p->cd();

  const char* title = Form("%s%sLeadJtPt_pp_%s", algType_pp[setNum], Sub_Lead, fileTag1);

  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, setDelPhiCut);


  getTree_p->Project(name, Form("Alg%sLeadJtPt[%d]", Sub_Lead, setNum), setCut && thirdJtVeto && delPhiCut);
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


void makeJtEtaHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Float_t histLow, Float_t histHi, const char* Sub_Lead = "")
{
  inFile_p->cd();

  const char* title = Form("%s%sLeadJtEta_pp_%s", algType_pp[setNum], Sub_Lead, fileTag1);

  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %f, %f)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut delPhiCut = Form("AlgJtDelPhi[%d] > %f", setNum, setDelPhiCut);

  getTree_p->Project(name, Form("Alg%sLeadJtEta[%d]", Sub_Lead, setNum), setCut && thirdJtVeto && delPhiCut);
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


void cfmDiJet_JtVarPlots_pp(const char* inName, const char* outName, Bool_t montecarlo = false, Bool_t Veto = false)
{
  //Some personal shorthand for myself, can be ignored as long as you feed it something for filetag

  if(!strcmp(inName, PPDataA)){
    std::cout << PPDataA << std::endl;
    fileTag1 = "PPDataA";
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

  if(Veto)
    thirdJtVeto = "!thirdJtVeto";

  std::cout << thirdJtVeto << std::endl;

  const char* subOrLead[2] = {"", "Sub"};
  Int_t leadOrSubBins[6] = {18, 120, 300, 25, 50, 300};

  for(Int_t algIter = 0; algIter < jetAlgMax; algIter++){

    makeAsymmHist(inFile1_p, inTree_p, outName, algIter, 10, 0, 1);
    makePtRatioHist(inFile1_p, inTree_p, outName, algIter, 10, 0, 1);
    makeAsymmHist_ThirdJet(inFile1_p, inTree_p, outName, algIter, 15, -.5, 1.);
    makeAsymmDelPhiHist(inFile1_p, inTree_p, outName, algIter);
    makeRatioDelPhiHist(inFile1_p, inTree_p, outName, algIter);

    makeRatioDelPhiHist_ThirdJet(inFile1_p, inTree_p, outName, algIter);

    makeDelPhiHist(inFile1_p, inTree_p, outName, algIter, 20, 5*TMath::Pi()/6, TMath::Pi());

    for(Int_t subIter = 0; subIter < 2; subIter++){
      makeJtPtHist(inFile1_p, inTree_p, outName, algIter, leadOrSubBins[subIter*3], leadOrSubBins[subIter*3 + 1], leadOrSubBins[subIter*3 + 2], subOrLead[subIter]);
      makeJtEtaHist(inFile1_p, inTree_p, outName, algIter, 10, -2.0, 2.0, subOrLead[subIter]);
    }

  }

  inFile1_p->Close();
  delete inFile1_p;
  return;
}
