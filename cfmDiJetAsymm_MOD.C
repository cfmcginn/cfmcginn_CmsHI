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



//append to every histo so know which sample via shorthand on workblog                                                          

const char* fileTag;

//shorthands, w/ _CFMSKIM.h                                                                                                                      

const char* Di80a = "Pythia80_HydjetDrum_mix01_HiForest2_v20_CFMSKIM.root";

const char* Di80b = "Dijet80_HydjetDrum_v27_mergedV1_CFMSKIM.root";

const char* Di80c = "HydjetDrum_Pyquen_Dijet80_Embedded_d20140122_Track7_v2_CFMSKIM.root";

const char* Di80d = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";                                      \


const char* Di80e = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet14_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";

const char* Di80f = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_50k_CFMSKIM.root";

const char* Di80g = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_300k_CFMSKIM.root";

const char* Di100a = "Dijet100_HydjetDrum_v27_mergedV1_CFMSKIM.root";

const char* EmDi80a = "PbPb_pythiaHYDJET_forest_EmEnrichedDijet80_CFMSKIM.root";

const char* DataA = "Track8_Jet17_GR_R_53_LV6_SUB_0_CFMSKIM.root";

const char* DataB = "hiForest_Jet80or95_GR_R_53_LV6_02Mar2014_1300CET_Track8_Jet15_0_1200k_CFMSKIM.root";

const char* DataC = "hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21_0_700k_CFMSKIM.root";



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

  const char* leadJt = Form("algLeadJtEta[%d]", setNum);
  const char* subLeadJt = Form("algSubLeadJtEta[%d]", setNum);

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


void makeAsymmHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Int_t histLow, Int_t histHi, Int_t centLow, Int_t centHi)
{
  inFile_p->cd();

  const char* title = Form("%sAsymm_%d%d_%s", algType[setNum], (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag);

  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);

  getTree_p->Project(name, Form("AlgJtAsymm[%d]", setNum), centCut && setCut);
  asymmHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  asymmHist_p->Sumw2();
  niceTH1(asymmHist_p, .65, 0., 405, 506);

  asymmHist_p->SetYTitle("Event Fraction");
  asymmHist_p->SetXTitle("A_{J} = (p_{T,1} - p_{T,2})/(p_{T,1} + p_{T,2})");

  outFile_p = new TFile(outName, "UPDATE");
  asymmHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
}


void makeDelPhiHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Float_t histLow, Float_t histHi, Int_t centLow, Int_t centHi)
{
  inFile_p->cd();

  const char* title = Form("%sDelPhi_%d%d_%s", algType[setNum], (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag);

  TH1F* delPhiHist_p;

  TString name = Form("%s_h(%d, %f, %f)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);

  getTree_p->Project(name, Form("AlgJtDelPhi[%d]", setNum), centCut && setCut);
  delPhiHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  delPhiHist_p->Sumw2();
  niceTH1(delPhiHist_p, 1., .001, 405, 506);

  delPhiHist_p->SetYTitle("Event Fraction");
  delPhiHist_p->SetXTitle("#Delta #phi_{1,2}");

  outFile_p = new TFile(outName, "UPDATE");
  delPhiHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
}


void makeJtPtHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Int_t histLow, Int_t histHi, Int_t centLow, Int_t centHi, const char* Sub_Lead = "")
{
  inFile_p->cd();

  const char* title = Form("%s%sLeadJtPt_%d%d_%s", algType[setNum], Sub_Lead, (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag);

  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);

  getTree_p->Project(name, Form("Alg%sLeadJtPt[%d]", Sub_Lead, setNum), centCut && setCut);
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
}


void makeJtEtaHist(TFile* inFile_p, TTree* getTree_p, const char* outName, Int_t setNum, Int_t nBins, Float_t histLow, Float_t histHi, Int_t centLow, Int_t centHi, const char* Sub_Lead = "")
{
  inFile_p->cd();

  const char* title = Form("%s%sLeadJtEta_%d%d_%s", algType[setNum], Sub_Lead, (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5), fileTag);

  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %f, %f)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);

  getTree_p->Project(name, Form("Alg%sLeadJtEta[%d]", Sub_Lead, setNum), centCut && setCut);
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
}


void addJtHistToPanel(TFile* file_p, TCanvas* canv_p, Int_t setNum, const char* jtVarIn, Int_t centLow, Int_t centHi, Int_t pos, Bool_t logY = false, const char* subOrLead = "")
{
  TH1F* hist_p = (TH1F*)file_p->Get(Form("%s%s%s_%d%d_%s_h", algType[setNum], subOrLead, jtVarIn, centLow, centHi, fileTag));
  canv_p->cd(pos);

  if(logY)
    gPad->SetLogy();

  if(pos == 4 || pos == 6 || pos == 3 || pos == 1)
    hist_p->SetXTitle("");

  if(pos == 2 || pos == 3 || pos == 5 || pos == 6)
    hist_p->SetYTitle("");

  hist_p->Draw("E1");

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->DrawLatex(.7, .3, Form("%d-%d%%", centLow, centHi));

  if(pos == 1)
    label_p->DrawLatex(.20, .85, Form("%s %s", algType[setNum], jtVarIn));

  if(pos == 2)
    label_p->DrawLatex(.20, .875, Form("Anti-k_{T} %s, R = 0.3", algType[setNum]));

  label_p->DrawLatex(.20, .80, "p_{T,1} > 120, p_{T,2} > 50 GeV/c");
  label_p->DrawLatex(.20, .725, "#Delta #phi_{1,2} > 2#pi/3");
  label_p->DrawLatex(.20, .65, "|#eta|_{1,2} < 1.6");

  delete label_p;
}


void makeJtVarPanel(const char* fileName1, Int_t setNum, const char* jtVarIn, Bool_t logY = false, const char* subOrLead = "", const char* fileName2 = "")
{
  TCanvas* jtVarPanel_p = new TCanvas(Form("%s%sPanel_%s", algType[setNum], jtVarIn, fileTag), Form("%s%sPanel_%s", algType[setNum], jtVarIn, fileTag), 1);
  jtVarPanel_p->Divide(3,2,0,0);

  TFile* panelFile_p = new TFile(fileName1, "UPDATE");

  TLegend* leg = new TLegend(.20, .70, .40, .80);

  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(.06);
  leg->SetBorderSize(0);

  addJtHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 70, 100, 1, logY, subOrLead);
  addJtHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 50, 70, 2, logY, subOrLead);
  addJtHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 30, 50, 3, logY, subOrLead);
  addJtHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 20, 30, 4, logY, subOrLead);
  addJtHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 10, 20, 5, logY, subOrLead);
  addJtHistToPanel(panelFile_p, jtVarPanel_p, setNum, jtVarIn, 0, 10, 6, logY, subOrLead);

  claverCanvasSaving(jtVarPanel_p, Form("../pngDir/%s%sPanel_%s_c", algType[setNum], jtVarIn, fileTag), "png");
  jtVarPanel_p->Write();
  delete leg;

  panelFile_p->Close();
  delete panelFile_p;
  delete jtVarPanel_p;
}


void cfmDiJetAsymm(const char* inName = "inFile_CFMHIST_GAMMA.root", bool montecarlo = 0, const char* outName = "outFile_CFMHIST_GAMMA.root")
{
  //Some personal shorthand for myself, can be ignored as long as you feed it something for filetag

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
  else if(!strcmp(inName, Di80e)){
    std::cout << Di80e << std::endl;
    fileTag = "Di80e";
  }
  else if(!strcmp(inName, Di80f)){
    std::cout << Di80f << std::endl;
    fileTag = "Di80f";
  }
  else if(!strcmp(inName, Di80g)){
    std::cout << Di80g << std::endl;
    fileTag = "Di80g";
  }
  else if(!strcmp(inName, DataA)){
    std::cout << DataA << std::endl;
    fileTag = "DataA";
  }
  else if(!strcmp(inName, DataB)){
    std::cout << DataB << std::endl;
    fileTag = "DataB";
  }
  else if(!strcmp(inName, DataC)){
    std::cout << DataC << std::endl;
    fileTag = "DataC";
  }

  std::cout << "Filetag is: " << fileTag << std::endl;

  inFile1_p = new TFile(inName, "READ");
  inTree_p = (TTree*)inFile1_p->Get("jetTree");
  inTree_p->AddFriend("trackTree");

  Int_t jetAlgMax = 4;

  if(montecarlo){
    inTree_p->AddFriend("genTree");
    jetAlgMax = 5;
  }

  Int_t centLow[8] = {0, 0, 20, 40, 60, 100, 100, 140};
  Int_t centHi[8] = {199, 19, 39, 59, 99, 199, 139, 199};

  const char* subOrLead[2] = {"", "Sub"};
  Int_t leadOrSubBins[6] = {18, 120, 300, 25, 50, 300};

  const char* jtVar[4] = {"Asymm", "DelPhi", "LeadJtPt", "LeadJtEta"};

  for(Int_t algIter = 0; algIter < jetAlgMax; algIter++){

    for(Int_t centIter = 0; centIter < 8; centIter++){
      makeAsymmHist(inFile1_p, inTree_p, outName, algIter, 10, 0, 1, centLow[centIter], centHi[centIter]);
      makeDelPhiHist(inFile1_p, inTree_p, outName, algIter, 10, 2*TMath::Pi()/3, TMath::Pi(), centLow[centIter], centHi[centIter]);

      for(Int_t subIter = 0; subIter < 2; subIter++){
	makeJtPtHist(inFile1_p, inTree_p, outName, algIter, leadOrSubBins[subIter*3], leadOrSubBins[subIter*3 + 1], leadOrSubBins[subIter*3 + 2], centLow[centIter], centHi[centIter], subOrLead[subIter]);
	makeJtEtaHist(inFile1_p, inTree_p, outName, algIter, 10, -1.6, 1.6, centLow[centIter], centHi[centIter], subOrLead[subIter]);
      }
      
    }

    makeJtVarPanel(outName, algIter, jtVar[0]);
  }

  inFile1_p->Close();
  delete inFile1_p;
}
