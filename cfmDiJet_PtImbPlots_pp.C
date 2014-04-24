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
//const char* PPDataA = "PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82_CFMSKIM.root";

const char* PPMC80A = "pt80_pp2013_P01_prod22_v81_merged_forest_0_CFMSKIM.root";

const char* PPMC80B = "QCDpT80_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet22_1GeVcut_merged_0_CFMSKIM.root";

const char* PPMC120A = "pt120_pp2013_P01_prod22_v81_merged_forest_0_CFMSKIM_20140323_4_0.root";

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


void makeImbAsymmGraph(TTree* getTree_p, const char* outName, const char* gorr, Int_t setNum, const char* perpProj, const char* CNC, const char* FPT, Int_t graphLow, Int_t graphHi, const char* GLN = "N", const char* Corr = "")
{
  inFile1_p->cd();

  Int_t setCorrNum = setNum;
  if(strcmp("", Corr) != 0)
    setCorrNum = setNum + 5;

  const char* title = Form("%s%sImbAsymm%s%s%s%s_%s_%s_g", gorr, algType[setNum], perpProj, CNC, FPT, Corr, GLN, fileTag1);

  TGraphErrors* imbAsymmGraph_p = new TGraphErrors(4);
  imbAsymmGraph_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTGraphErrors(imbAsymmGraph_p, graphHi, graphLow);

  TH1F* getHist_p;

  TString var = Form("%sAlgImb%s%s%s[%d]", gorr, perpProj, CNC, FPT, setCorrNum);

  TCut setCut = makeSetCut(setNum);
  TCut etaCut = makeEtaCut(setNum, 1.6, GLN);

  TCut phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);
  if(strcmp(CNC, "") != 0){
    etaCut = makeEtaCut(setNum, 1.6, GLN);
    phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);
  }

  TCut jetLCut = Form("AlgLeadJtPt[%d] > 120*.95", setNum);

  //  TCut isQuark = Form("isQuarkJet[%d]", setNum);
  //  TCut isGluon = Form("isGluonJet[%d]", setNum);


  const char* name1[4] = {"0_1(100000, -100000, 100000)", "1_2(100000, -100000, 100000)", "2_3(100000, -100000, 100000)", "3_5(100000, -100000, 100000)"};
  const char* name2[4] = {"0_1", "1_2", "2_3", "3_5"};
  Float_t asymmBins[5] = {.00, .11, .22, .33, 1.};
  Float_t point[4] = {.055, .165, .275, .415};
  Float_t xErr[4] = {.055, .055, .055, .085};

  for(Int_t binIter = 0; binIter < 4; binIter++){
    TCut asymmCut = makeAsymmCut(setNum, asymmBins[binIter], asymmBins[binIter + 1]);

    if(strcmp("", CNC) != 0)
      std::cout << phiCut << ", " << fullVeto << std::endl;

    getTree_p->Project(name1[binIter], var, setCut && etaCut && phiCut && jetLCut && asymmCut && fullVeto);
    
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


void makeImbAsymmGraph_Tight(TTree* getTree_p, const char* outName, const char* gorr, Int_t setNum, const char* perpProj, const char* CNC, const char* FPT, Int_t graphLow, Int_t graphHi, const char* GLN = "N", const char* Corr = "")
{
  inFile1_p->cd();

  Int_t setCorrNum = setNum;
  if(strcmp("", Corr) != 0)
    setCorrNum = setNum + 5;

  const char* title = Form("%s%sImbAsymmTight%s%s%s%s_%s_%s_g", gorr, algType[setNum], perpProj, CNC, FPT, Corr, GLN, fileTag1);


  TGraphErrors* imbAsymmGraph_p = new TGraphErrors(4);
  imbAsymmGraph_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTGraphErrors(imbAsymmGraph_p, graphHi, graphLow);

  TH1F* getHist_p;


  TString var = Form("%sAlgImb%s%s%s[%d]", gorr, perpProj, CNC, FPT, setCorrNum);

  TCut setCut = makeSetCut(setNum);
  TCut etaCut = makeEtaCut(setNum, .5, GLN);

  TCut phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);
  if(strcmp(CNC, "") != 0){
    phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);
    etaCut = makeEtaCut(setNum, .5, GLN);
  }

  TCut jetLCut = Form("AlgLeadJtPt[%d] > 120*.95", setNum);

  //  TCut isQuark = Form("isQuarkJet[%d]", setNum);
  //  TCut isGluon = Form("isGluonJet[%d]", setNum);

  const char* name1[8] = {"0_h(100000, -100000, 100000)", "1_h(100000, -100000, 100000)", "2_h(100000, -100000, 100000)", "3_h(100000, -100000, 100000)", "4_h(100000, -100000, 100000)", "5_h(100000, -100000, 100000)", "6_h(100000, -100000, 100000)", "7_h(100000, -100000, 100000)"};
  const char* name2[8] = {"0_h", "1_h", "2_h", "3_h", "4_h", "5_h", "6_h", "7_h"};
  Float_t asymmBins[9] = {.00, .055, .11, .165, .22, .275, .33, .415, 1.};
  Float_t point[8] = {.0275, .0825, .1375, .1925, .2475, .3025, .3725, .4575};
  Float_t xErr[8] = {.0275, .0275, .0275, .0275, .0275, .0275, .0425, .0425};

  for(Int_t binIter = 0; binIter < 8; binIter++){
    TCut asymmCut = makeAsymmCut(setNum, asymmBins[binIter], asymmBins[binIter + 1]);

    getTree_p->Project(name1[binIter], var, setCut && etaCut && phiCut && jetLCut && asymmCut && fullVeto);
    
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


void cfmDiJet_PtImbPlots_pp(const char* inName, const char* outName, Bool_t montecarlo = false, Bool_t Veto = false)
{
  //Some personal shorthand for myself, can be ignored as long as you feed it something for filetag

  TH1::SetDefaultSumw2();

  if(!strcmp(inName, PPDataA)){
    std::cout << PPDataA << std::endl;
    fileTag1 = "PPDataA";
  }
  else if(!strcmp(inName, PPMC80A)){
    std::cout << PPMC80A << std::endl;
    fileTag1 = "PPMC80A";
  }
  else if(!strcmp(inName, PPMC80B)){
    std::cout << PPMC80B << std::endl;
    fileTag1 = "PPMC80B";
  }
  else if(!strcmp(inName, PPMC120A)){
    std::cout << PPMC120A << std::endl;
    fileTag1 = "PPMC120A";
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

  for(Int_t algIter = 0; algIter < jetAlgMax; algIter++){
    if(algIter != 3 && algIter != 4) continue;

    if(Veto && algIter !=4)
      fullVeto = Form("AlgLeadRefPt[%d] > 0 && AlgSubLeadRefPt[%d] > 0 && AlgLeadRefPt[%d] > AlgSubLeadRefPt[%d] && AlgThirdJtPt[%d] < 50", algIter, algIter, algIter, algIter, algIter);
    if(Veto && algIter == 4)
      fullVeto = Form("AlgThirdJtPt[4] < 50");

    for(Int_t corrIter = 0; corrIter < 2; corrIter++){

      for(Int_t CNCIter = 0; CNCIter < 4; CNCIter++){
	makeImbAsymmGraph(inTree_p, outName, "r", algIter, "ProjA", CNC[CNCIter], "F", -60, 60, "N", corr[corrIter]);
	
	if(montecarlo && corrIter == 0){
	  makeImbAsymmGraph(inTree_p, outName, "g", algIter, "ProjA", CNC[CNCIter], "F", -60, 60, "N", "");
	}
   
	for(Int_t ptBinIter = 0; ptBinIter < 5; ptBinIter++){
	  makeImbAsymmGraph(inTree_p, outName, "r", algIter, "ProjA", CNC[CNCIter], ptBins[ptBinIter], -60, 60, "N", corr[corrIter]);

	  if(montecarlo && corrIter == 0){
	    makeImbAsymmGraph(inTree_p, outName, "g", algIter, "ProjA", CNC[CNCIter], ptBins[ptBinIter], -60, 60, "N", "");
	  }

	}

      }

    }

  }


  inFile1_p->Close();
  delete inFile1_p;
  return;
}
