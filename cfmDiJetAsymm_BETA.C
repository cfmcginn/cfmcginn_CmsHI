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

//auto ln -s title of address /net/hisrv0001/home/yenjie/scratch/tmp/Pythia80_HydjetDrum_mix01_HiForest2_v20.root                           
const char* Di80a = "Pythia80_HydjetDrum_mix01_HiForest2_v20_CFMSKIM.root";

//auto ln -s title of address /mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet80_HydjetDrum_v27_mergedV1.root                             
const char* Di80b = "Dijet80_HydjetDrum_v27_mergedV1_CFMSKIM.root";

//auto ln -s title of address /mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet100_HydjetDrum_v27_mergedV1.root                           
const char* Di100a = "Dijet100_HydjetDrum_v27_mergedV1_CFMSKIM.root";

//auto ln -s title of address /mnt/hadoop/cms/store/user/luck/PbPb_pythiaHYDJET_forest_EmEnrichedDijet/PbPb_pythiaHYDJET_forest_EmEnrichedDijet80.root                                                                                                                                   
const char* EmDi80a = "PbPb_pythiaHYDJET_forest_EmEnrichedDijet80_CFMSKIM.root";

const char* Di80c = "HydjetDrum_Pyquen_Dijet80_Embedded_d20140122_Track7_v2_CFMSKIM.root";

const char* Di80d = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";

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

Int_t checkLORS(const char* lors){
  if(strcmp(lors,"Lead") == 0 || strcmp(lors, "Sublead") == 0)    return 1;

  return 0;
}


Int_t checkPTPHIETA(const char* ptPhiEta){
  if(strcmp(ptPhiEta, "Pt") == 0 || strcmp(ptPhiEta, "Phi") == 0 || strcmp(ptPhiEta, "Eta") == 0)    return 1;

  return 0;
}


Int_t checkGORR(const char* gorr){
  if(strcmp(gorr, "g") == 0 || strcmp(gorr, "gR") == 0 || strcmp(gorr, "r") == 0)    return 1;

    return 0;
}


Int_t checkPFCalo(const char* PFCalo){
  if(strcmp(PFCalo, "PF") == 0 || strcmp(PFCalo, "Calo") == 0 || strcmp(PFCalo, "") == 0) return 1;

  return 0;
}


Int_t checkPERPPROJ(const char* perpProj){
  if(strcmp(perpProj, "Proj") == 0 || strcmp(perpProj, "Perp") == 0)    return 1;

  return 0;
}


Int_t checkFHL(const char* FHL){
  if(strcmp(FHL, "F") == 0 || strcmp(FHL, "H") == 0 || strcmp(FHL, "L") == 0)    return 1;

  return 0;
}


Int_t checkGLN(const char* GLN){
  if(strcmp(GLN, "G") == 0 || strcmp(GLN, "L") == 0 || strcmp(GLN, "N") == 0)    return 1;

  return 0;
}


Int_t checkChar(const char* lors = "Lead", const char* ptPhiEta = "Pt", const char* gorr = "r", const char* perpProj = "Proj", const char* FHL = "F", const char* GLN = "N"){
  if(!checkLORS(lors)){
    std::cout << "LORS not specified; Return" << std::endl;
    return 0;
  }

  if(!checkPTPHIETA(ptPhiEta)){
    std::cout << "PTPHIETA not specified; Return" << std::endl;
    return 0;
  }

  if(!checkGORR(gorr)){
    std::cout << "GORR not specified; Return" << std::endl;
    return 0;
  }

  if(!checkPERPPROJ(perpProj)){
    std::cout << "PERPPROJ not specified; Return" << std::endl;
    return 0;
  }

  if(!checkFHL(FHL)){
    std::cout << "FHL not specified; Return" << std::endl;
    return 0;
  }

  if(!checkGLN(GLN)){
    std::cout << "GLN not specified; Return" << std::endl;
    return 0;
  }

  return 1;
}


TCut makeSetCut(const char* gorr, const char* PFCalo)
{
  TCut setCut = "";
  if(strcmp(gorr, "g") == 0 || strcmp(gorr, "gR") == 0)
    setCut = Form("truthSet == 1");
  else if(strcmp(gorr, "r") == 0){
    if(strcmp(PFCalo, "PF") == 0)
      setCut = Form("recoPFSet == 1");
    else if(strcmp(PFCalo, "Calo") == 0)
      setCut = Form("recoCaloSet == 1");
  }
  else
    std::cout << "setCut unspecified, returning blank cut" << std::endl;

  return setCut;
}


TCut makeCentCut(Int_t centLow, Int_t centHi)
{
  TCut centCut = "";
  if(centLow >= 0 && centHi >= centLow && centHi <= 39)
    centCut = Form("hiBin >= %d && hiBin <= %d", centLow, centHi);
  else
    std::cout << "Warning: Cent Cut empty, centLow/centHi incorrectly specified" << std::endl;

  return centCut;
}


TCut makeAsymmCut(const char* grpfcalo, Float_t asymmLow, Float_t asymmHi)
{
  TCut asymmCut = "";
  const char* asymmJetType = Form("(%sLeadJtPt - %sSubLeadJtPt)/(%sLeadJtPt + %sSubLeadJtPt)", grpfcalo, grpfcalo, grpfcalo, grpfcalo);

  //For backwards compat, will change
  if(strcmp(grpfcalo, "gRPF") == 0)
    asymmJetType = Form("(%sLeadJtPt - %sSubLeadJtPt)/(%sLeadJtPt + %sSubLeadJtPt)", "g", "g", "g", "g");

  if(asymmLow >= 0 && asymmHi >= asymmLow && asymmHi <= 1)
    asymmCut = Form("%s > %f && %s < %f ", asymmJetType, asymmLow, asymmJetType, asymmHi);
  else
    std::cout << "Warning: Asymm Cut empty, asymmLow/asymmHi incorrectly specified" << std::endl;

  return asymmCut;
}


void makeJtCompHist(TTree* getTree_p, const char* outName, const char* lors, const char* ptPhiEta, Int_t xBins, Float_t xLow, Float_t xHi, Int_t yBins, Float_t yLow, Float_t yHi)
{
  if(!checkChar(lors, ptPhiEta)) return;

  inFile_p->cd();

  const char* title = Form("rVGJt%s_%s_%s", ptPhiEta, lors, fileTag);

  TString name = Form("%s_h", title);
  TH2F* jtCompHist_p = new TH2F(name, name, xBins, xLow, xHi, yBins, yLow, yHi);
  TCut histCut = Form("rPF%sJt%s > %f && rPF%sJt%s < %f && g%sJt%s > %f && g%sJt%s < %f", lors, ptPhiEta, xLow, lors, ptPhiEta, xHi, lors, ptPhiEta, yLow, lors, ptPhiEta, yHi);
  getTree_p->Project(name, Form("rPF%sJt%s:g%sJt%s", lors, ptPhiEta, lors, ptPhiEta), histCut);

  outFile_p = new TFile(outName, "UPDATE");
  jtCompHist_p->Write();
  outFile_p->Close();
  delete outFile_p;
}


void makeAsymmHist(TTree* getTree_p, const char* outName, const char* gorr, const char* PFCalo, Int_t nBins, Int_t histLow, Int_t histHi, Int_t centLow, Int_t centHi)
{
  if(!checkChar("Lead", "Pt", gorr)) return;

  const char* grpfcalo = Form("%s%s", gorr, PFCalo);

  inFile_p->cd();

  const char* title = Form("%sAsymm_%d%d_%s", grpfcalo, (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5), fileTag);

  //  TCanvas* asymmCanvas_p = new TCanvas(Form("%s_c", title), Form("%s_c", title), 1);
  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);
  TCut setCut = makeSetCut(gorr, PFCalo);
  TCut centCut = makeCentCut(centLow, centHi);
  getTree_p->Project(name, Form("(%sLeadJtPt - %sSubLeadJtPt)/(%sLeadJtPt + %sSubLeadJtPt)", grpfcalo, grpfcalo, grpfcalo, grpfcalo), centCut && setCut);
  asymmHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  niceTH1(asymmHist_p, .6, 0., 405, 506);

  asymmHist_p->SetYTitle("Event Fraction");
  asymmHist_p->SetXTitle("A_{J} = (p_{T,1} - p_{T,2})/(p_{T,1} + p_{T,2})");
  //  asymmHist_p->Draw();
  //  cent_p->SetNDC();
  //  cent_p->Draw();

  outFile_p = new TFile(outName, "UPDATE");
  asymmHist_p->Write(Form("%s_h", title));
  //  asymmCanvas_p->Write();
  outFile_p->Close();

  delete outFile_p;
  //  delete cent_p;
  //  delete asymmCanvas_p;
}


void addHistToPanel(TFile* file_p, TH1F* hist_p, TCanvas* canv_p, const char* gorr, const char* PFCalo, Int_t centLow, Int_t centHi, Int_t pos)
{
  if(!checkChar("Lead", "Pt" , gorr)) return;

  const char* grpfcalo = Form("%s%s", gorr, PFCalo);

  hist_p = (TH1F*)file_p->Get(Form("%sAsymm_%d%d_%s_h", grpfcalo, centLow, centHi, fileTag));
  canv_p->cd(pos);
  hist_p->Draw();

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->DrawLatex(.6, .3, Form("%d-%d%%", centLow, centHi));

  if(strcmp(gorr, "g") == 0 && pos == 1){
    label_p->DrawLatex(.3, .85, Form("Truth DiJet Asymmetry, Truth Set"));
    label_p->DrawLatex(.3, .8, Form("%s", fileTag));
  }
  else if(strcmp(gorr, "r") == 0 && pos == 1){
    label_p->DrawLatex(.3, .85, Form("Reco DiJet Asymmetry, Reco Set"));
    label_p->DrawLatex(.3, .8, Form("%s", fileTag));
  }
  else if(strcmp(gorr, "gR") == 0 && pos == 1){
    label_p->DrawLatex(.3, .85, Form("Reco DiJet Asymmetry, Truth Set"));
    label_p->DrawLatex(.3, .8, Form("%s", fileTag));
  }


  delete label_p;
}


void makeAsymmPanel(const char* fileName, const char* gorr, const char* PFCalo)
{
  if(!checkChar("Lead", "Pt",gorr)) return;

  const char* grpfcalo = Form("%s%s", gorr, PFCalo);

  TFile* panelFile_p = new TFile(fileName, "UPDATE");
  TH1F* asymmHist_p;

  TCanvas* asymmPanel_p = new TCanvas(Form("%sAsymmPanel_%s_c", grpfcalo, fileTag), Form("%sAsymmPanel_%s_c", grpfcalo, fileTag), 1);
  asymmPanel_p->Divide(3, 2, 0, 0);

  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, PFCalo, 0, 100, 1);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, PFCalo, 50, 100, 2);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, PFCalo, 30, 50, 3);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, PFCalo, 20, 30, 4);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, PFCalo, 10, 20, 5);
  addHistToPanel(panelFile_p, asymmHist_p, asymmPanel_p, gorr, PFCalo, 0, 10, 6);

  asymmPanel_p->Write();
  claverCanvasSaving(asymmPanel_p, Form("../pngDir/%sAsymmPanel_%s", grpfcalo, fileTag), "png");                                                 
  panelFile_p->Close();
  delete panelFile_p;
  delete asymmPanel_p;
}

//Test

void makePtProjHist(TTree* getTree_p, const char* outName)
{
  inFile_p->cd();

  //  TCanvas* asymmCanvas_p = new TCanvas(Form("%s_c", title), Form("%s_c", title), 1);
  gStyle->SetOptStat(1110);
  TH1F* trkPtRPFHist_0_1_030_Symm_p;
  TH1F* trkPtRPFHist_1_8_030_Symm_p;
  TH1F* trkPtRPFHist_0_1_30100_Symm_p;
  TH1F* trkPtRPFHist_0_1_030_Asymm_p;

  TH1F* trkPtGRPFHist_0_1_030_Symm_p;
  TH1F* trkPtGRPFHist_1_8_030_Symm_p;
  TH1F* trkPtGRPFHist_0_1_30100_Symm_p;
  TH1F* trkPtGRPFHist_0_1_030_Asymm_p;


  TH1F* genPtJTHist_0_1_030_Symm_p;
  TH1F* genPtJTHist_1_8_030_Symm_p;
  TH1F* genPtJTHist_0_1_30100_Symm_p;
  TH1F* genPtJTHist_0_1_030_Asymm_p;



  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextSize(.03);

  TString name = Form("trkPtRPF_0_1_030_Symm");
  TCut setCut = makeSetCut("r", "PF");
  TCut centCut = makeCentCut(0, 11);
  TCut asymmCut = makeAsymmCut("rPF", 0, .1);
  TCut ptCut = "trkPt < 1";
  getTree_p->Project(name, "trkPtRPF", centCut && setCut && ptCut && asymmCut);
  trkPtRPFHist_0_1_030_Symm_p = (TH1F*)inFile_p->Get(Form("trkPtRPF_0_1_030_Symm"));
  handsomeTH1(trkPtRPFHist_0_1_030_Symm_p);
  trkPtRPFHist_0_1_030_Symm_p->SetXTitle("trkPt Proj. onto Lead akPu3PF (GeV/c)");
  trkPtRPFHist_0_1_030_Symm_p->SetYTitle("N trks");

  name = Form("trkPtRPF_1_8_030_Symm");
  setCut = makeSetCut("r", "PF");
  centCut = makeCentCut(0, 11);
  asymmCut = makeAsymmCut("rPF", 0, .1);
  ptCut = "trkPt > 1 && trkPt < 8";
  getTree_p->Project(name, "trkPtRPF", centCut && setCut && ptCut && asymmCut);
  trkPtRPFHist_1_8_030_Symm_p = (TH1F*)inFile_p->Get(Form("trkPtRPF_1_8_030_Symm"));
  handsomeTH1(trkPtRPFHist_1_8_030_Symm_p);
  trkPtRPFHist_1_8_030_Symm_p->SetXTitle("trkPt Proj. onto Lead akPu3PF (GeV/c)");
  trkPtRPFHist_1_8_030_Symm_p->SetYTitle("N trks");

  name = Form("trkPtRPF_0_1_30100_Symm");
  setCut = makeSetCut("r", "PF");
  centCut = makeCentCut(12, 39);
  asymmCut = makeAsymmCut("rPF", 0, .1);
  ptCut = "trkPt < 1";
  getTree_p->Project(name, "trkPtRPF", centCut && setCut && ptCut && asymmCut);
  trkPtRPFHist_0_1_30100_Symm_p = (TH1F*)inFile_p->Get(Form("trkPtRPF_0_1_30100_Symm"));
  handsomeTH1(trkPtRPFHist_0_1_30100_Symm_p);
  trkPtRPFHist_0_1_30100_Symm_p->SetXTitle("trkPt Proj. onto Lead akPu3PF (GeV/c)");
  trkPtRPFHist_0_1_30100_Symm_p->SetYTitle("N trks");


  name = Form("trkPtRPF_0_1_030_Asymm");
  setCut = makeSetCut("r", "PF");
  centCut = makeCentCut(0, 11);
  asymmCut = makeAsymmCut("rPF", .3, 1.);
  ptCut = "trkPt < 1";
  getTree_p->Project(name, "trkPtRPF", centCut && setCut && ptCut && asymmCut);
  trkPtRPFHist_0_1_030_Asymm_p = (TH1F*)inFile_p->Get(Form("trkPtRPF_0_1_030_Asymm"));
  handsomeTH1(trkPtRPFHist_0_1_030_Asymm_p);
  trkPtRPFHist_0_1_030_Asymm_p->SetXTitle("trkPt Proj. onto Lead akPu3PF (GeV/c)");
  trkPtRPFHist_0_1_030_Asymm_p->SetYTitle("N trks");


  name = Form("trkPtGRPF_0_1_030_Symm");
  setCut = makeSetCut("gR", "PF");
  centCut = makeCentCut(0, 11);
  asymmCut = makeAsymmCut("gRPF", 0, .1);
  ptCut = "trkPt < 1";
  getTree_p->Project(name, "trkPtGRPF", centCut && setCut && ptCut && asymmCut);
  trkPtGRPFHist_0_1_030_Symm_p = (TH1F*)inFile_p->Get(Form("trkPtGRPF_0_1_030_Symm"));
  handsomeTH1(trkPtGRPFHist_0_1_030_Symm_p);
  trkPtGRPFHist_0_1_030_Symm_p->SetXTitle("trkPt Proj. onto Lead genJet (GeV/c)");
  trkPtGRPFHist_0_1_030_Symm_p->SetYTitle("N trks");

  name = Form("trkPtGRPF_1_8_030_Symm");
  setCut = makeSetCut("gR", "PF");
  centCut = makeCentCut(0, 11);
  asymmCut = makeAsymmCut("gRPF", 0, .1);
  ptCut = "trkPt > 1 && trkPt < 8";
  getTree_p->Project(name, "trkPtGRPF", centCut && setCut && ptCut && asymmCut);
  trkPtGRPFHist_1_8_030_Symm_p = (TH1F*)inFile_p->Get(Form("trkPtGRPF_1_8_030_Symm"));
  handsomeTH1(trkPtGRPFHist_1_8_030_Symm_p);
  trkPtGRPFHist_1_8_030_Symm_p->SetXTitle("trkPt Proj. onto Lead genJet (GeV/c)");
  trkPtGRPFHist_1_8_030_Symm_p->SetYTitle("N trks");

  name = Form("trkPtGRPF_0_1_30100_Symm");
  setCut = makeSetCut("gR", "PF");
  centCut = makeCentCut(12, 39);
  asymmCut = makeAsymmCut("gRPF", 0, .1);
  ptCut = "trkPt < 1";
  getTree_p->Project(name, "trkPtGRPF", centCut && setCut && ptCut && asymmCut);
  trkPtGRPFHist_0_1_30100_Symm_p = (TH1F*)inFile_p->Get(Form("trkPtGRPF_0_1_30100_Symm"));
  handsomeTH1(trkPtGRPFHist_0_1_30100_Symm_p);
  trkPtGRPFHist_0_1_30100_Symm_p->SetXTitle("trkPt Proj. onto Lead genJet (GeV/c)");
  trkPtGRPFHist_0_1_30100_Symm_p->SetYTitle("N trks");

  name = Form("trkPtGRPF_0_1_030_Asymm");
  setCut = makeSetCut("gR", "PF");
  centCut = makeCentCut(0, 11);
  asymmCut = makeAsymmCut("gRPF", .3, 1.);
  ptCut = "trkPt < 1";
  getTree_p->Project(name, "trkPtGRPF", centCut && setCut && ptCut && asymmCut);
  trkPtGRPFHist_0_1_030_Asymm_p = (TH1F*)inFile_p->Get(Form("trkPtGRPF_0_1_030_Asymm"));
  handsomeTH1(trkPtGRPFHist_0_1_030_Asymm_p);
  trkPtGRPFHist_0_1_030_Asymm_p->SetXTitle("trkPt Proj. onto Lead genJet (GeV/c)");
  trkPtGRPFHist_0_1_030_Asymm_p->SetYTitle("N trks");



  name = Form("genPtJT_0_1_030_Symm");
  setCut = makeSetCut("g", "");
  centCut = makeCentCut(0, 11);
  asymmCut = makeAsymmCut("g", 0, .1);
  ptCut = "genPt < 1";
  getTree_p->Project(name, "genPtJT", centCut && setCut && ptCut && asymmCut);
  genPtJTHist_0_1_030_Symm_p = (TH1F*)inFile_p->Get(Form("genPtJT_0_1_030_Symm"));
  handsomeTH1(genPtJTHist_0_1_030_Symm_p);
  genPtJTHist_0_1_030_Symm_p->SetXTitle("genPt Proj. onto Lead genJet (GeV/c)");
  genPtJTHist_0_1_030_Symm_p->SetYTitle("N Particles");


  name = Form("genPtJT_1_8_030_Symm");
  setCut = makeSetCut("g", "");
  centCut = makeCentCut(0, 11);
  asymmCut = makeAsymmCut("g", 0, .1);
  ptCut = "genPt < 8 && genPt > 1";
  getTree_p->Project(name, "genPtJT", centCut && setCut && ptCut && asymmCut);
  genPtJTHist_1_8_030_Symm_p = (TH1F*)inFile_p->Get(Form("genPtJT_1_8_030_Symm"));
  handsomeTH1(genPtJTHist_1_8_030_Symm_p);
  genPtJTHist_1_8_030_Symm_p->SetXTitle("genPt Proj. onto Lead genJet (GeV/c)");
  genPtJTHist_1_8_030_Symm_p->SetYTitle("N Particles");


  name = Form("genPtJT_0_1_30100_Symm");
  setCut = makeSetCut("g", "");
  centCut = makeCentCut(12, 39);
  asymmCut = makeAsymmCut("g", 0, .1);
  ptCut = "genPt < 1";
  getTree_p->Project(name, "genPtJT", centCut && setCut && ptCut && asymmCut);
  genPtJTHist_0_1_30100_Symm_p = (TH1F*)inFile_p->Get(Form("genPtJT_0_1_30100_Symm"));
  handsomeTH1(genPtJTHist_0_1_30100_Symm_p);
  genPtJTHist_0_1_30100_Symm_p->SetXTitle("genPt Proj. onto Lead genJet (GeV/c)");
  genPtJTHist_0_1_30100_Symm_p->SetYTitle("N Particles");


  name = Form("genPtJT_0_1_030_Asymm");
  setCut = makeSetCut("g", "");
  centCut = makeCentCut(0, 11);
  asymmCut = makeAsymmCut("g", .3, 1.);
  ptCut = "genPt < 1";
  getTree_p->Project(name, "genPtJT", centCut && setCut && ptCut && asymmCut);
  genPtJTHist_0_1_030_Asymm_p = (TH1F*)inFile_p->Get(Form("genPtJT_0_1_030_Asymm"));
  handsomeTH1(genPtJTHist_0_1_030_Asymm_p);
  genPtJTHist_0_1_030_Asymm_p->SetXTitle("genPt Proj. onto Lead genJet (GeV/c)");
  genPtJTHist_0_1_030_Asymm_p->SetYTitle("N Particles");


  outFile_p = new TFile(outName, "UPDATE");
  //  trkPtRPFHist_0_1_030_Symm_p->Write(Form("trkPtRPF_0_1_030_Symm_h"));
  //  trkPtRPFHist_1_8_030_Symm_p->Write(Form("trkPtRPF_1_8_030_Symm_h"));
  //  trkPtRPFHist_0_1_30100_Symm_p->Write(Form("trkPtRPF_0_1_30100_Symm_h"));
  //  trkPtRPFHist_0_1_030_Asymm_p->Write(Form("trkPtRPF_0_1_030_Asymm_h"));


  TCanvas* trkPtProjCanvas_p = new TCanvas("trkPtProjCanvas", "trkPtProjCanvas", 1);
  trkPtRPFHist_0_1_030_Symm_p->Draw();
  label_p->DrawLatex(.25, .35, "Lead akPu3PF p_{T} > 120 GeV/c");
  label_p->DrawLatex(.25, .32, "Sublead akPu3PF p_{T} > 50 GeV/c");
  label_p->DrawLatex(.25, .29, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.25, .26, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 0-30%; A_{J} < .1; .5 < trkPt < 1");
  trkPtProjCanvas_p->Write("trkPtRPF_0_1_030_Symm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/trkPtRPF_0_1_030_Symm_%s", fileTag), "png");

  trkPtRPFHist_1_8_030_Symm_p->Draw();
  label_p->DrawLatex(.19, .82, "Lead akPu3PF p_{T} > 120 GeV/c");
  label_p->DrawLatex(.19, .79, "Sublead akPu3PF p_{T} > 50 GeV/c");
  label_p->DrawLatex(.19, .76, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.19, .73, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 0-30%; A_{J} < .1; 1 < trkPt < 8");
  trkPtProjCanvas_p->Write("trkPtRPF_1_8_030_Symm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/trkPtRPF_1_8_030_Symm_%s", fileTag), "png");

  trkPtRPFHist_0_1_30100_Symm_p->Draw();
  label_p->DrawLatex(.25, .35, "Lead akPu3PF p_{T} > 120 GeV/c");
  label_p->DrawLatex(.25, .32, "Sublead akPu3PF p_{T} > 50 GeV/c");
  label_p->DrawLatex(.25, .29, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.25, .26, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 30-100%; A_{J} < .1; .5 < trkPt < 1");
  trkPtProjCanvas_p->Write("trkPtRPF_0_1_30100_Symm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/trkPtRPF_0_1_30100_Symm_%s", fileTag), "png");

  trkPtRPFHist_0_1_030_Asymm_p->Draw();
  label_p->DrawLatex(.25, .35, "Lead akPu3PF p_{T} > 120 GeV/c");
  label_p->DrawLatex(.25, .32, "Sublead akPu3PF p_{T} > 50 GeV/c");
  label_p->DrawLatex(.25, .29, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.25, .26, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 0-30%; A_{J} > .3; .5 < trkPt < 1");
  trkPtProjCanvas_p->Write("trkPtRPF_0_1_030_Asymm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/trkPtRPF_0_1_030_Asymm_%s", fileTag), "png");



  trkPtGRPFHist_0_1_030_Symm_p->Draw();
  label_p->DrawLatex(.25, .35, "Lead genJet p_{T} > 120 GeV/c");
  label_p->DrawLatex(.25, .32, "Sublead genJet p_{T} > 50 GeV/c");
  label_p->DrawLatex(.25, .29, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.25, .26, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 0-30%; A_{J} < .1; .5 < trkPt < 1");
  trkPtProjCanvas_p->Write("trkPtGRPF_0_1_030_Symm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/trkPtGRPF_0_1_030_Symm_%s", fileTag), "png");

  trkPtGRPFHist_1_8_030_Symm_p->Draw();
  label_p->DrawLatex(.19, .82, "Lead genJet p_{T} > 120 GeV/c");
  label_p->DrawLatex(.19, .79, "Sublead genJet p_{T} > 50 GeV/c");
  label_p->DrawLatex(.19, .76, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.19, .73, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 0-30%; A_{J} < .1; 1 < trkPt < 8");
  trkPtProjCanvas_p->Write("trkPtGRPF_1_8_030_Symm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/trkPtGRPF_1_8_030_Symm_%s", fileTag), "png");

  trkPtGRPFHist_0_1_30100_Symm_p->Draw();
  label_p->DrawLatex(.25, .35, "Lead genJet p_{T} > 120 GeV/c");
  label_p->DrawLatex(.25, .32, "Sublead genJet p_{T} > 50 GeV/c");
  label_p->DrawLatex(.25, .29, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.25, .26, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 30-100%; A_{J} < .1; .5 < trkPt < 1");
  trkPtProjCanvas_p->Write("trkPtGRPF_0_1_30100_Symm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/trkPtGRPF_0_1_30100_Symm_%s", fileTag), "png");

  trkPtGRPFHist_0_1_030_Asymm_p->Draw();
  label_p->DrawLatex(.25, .35, "Lead genJet p_{T} > 120 GeV/c");
  label_p->DrawLatex(.25, .32, "Sublead genJet p_{T} > 50 GeV/c");
  label_p->DrawLatex(.25, .29, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.25, .26, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 0-30%; A_{J} > .3; .5 < trkPt < 1");
  trkPtProjCanvas_p->Write("trkPtGRPF_0_1_030_Asymm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/trkPtGRPF_0_1_030_Asymm_%s", fileTag), "png");


  genPtJTHist_0_1_030_Symm_p->Draw();
  label_p->DrawLatex(.25, .35, "Lead genJet p_{T} > 120 GeV/c");
  label_p->DrawLatex(.25, .32, "Sublead genJet p_{T} > 50 GeV/c");
  label_p->DrawLatex(.25, .29, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.25, .26, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 0-30%; A_{J} < .1; .5 < genPt < 1");
  trkPtProjCanvas_p->Write("genPtJT_0_1_030_Symm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/genPtJT_0_1_030_Symm_%s", fileTag), "png");

  genPtJTHist_1_8_030_Symm_p->Draw();
  label_p->DrawLatex(.25, .35, "Lead genJet p_{T} > 120 GeV/c");
  label_p->DrawLatex(.25, .32, "Sublead genJet p_{T} > 50 GeV/c");
  label_p->DrawLatex(.25, .29, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.25, .26, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 0-30%; A_{J} < .1; 1 < genPt < 8");
  trkPtProjCanvas_p->Write("genPtJT_1_8_030_Symm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/genPtJT_1_8_030_Symm_%s", fileTag), "png");

  genPtJTHist_0_1_30100_Symm_p->Draw();
  label_p->DrawLatex(.25, .35, "Lead genJet p_{T} > 120 GeV/c");
  label_p->DrawLatex(.25, .32, "Sublead genJet p_{T} > 50 GeV/c");
  label_p->DrawLatex(.25, .29, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.25, .26, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 30-100%; A_{J} < .1; .5 < genPt < 1");
  trkPtProjCanvas_p->Write("genPtJT_0_1_30100_Symm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/genPtJT_0_1_30100_Symm_%s", fileTag), "png");

  genPtJTHist_0_1_030_Asymm_p->Draw();
  label_p->DrawLatex(.25, .35, "Lead genJet p_{T} > 120 GeV/c");
  label_p->DrawLatex(.25, .32, "Sublead genJet p_{T} > 50 GeV/c");
  label_p->DrawLatex(.25, .29, "Lead-SubLead #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.25, .26, "Lead/Sublead Jet abs(#eta) < 1.6");
  label_p->DrawLatex(.18, .91, "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0.root");
  label_p->DrawLatex(.18, .86, "Bin: Cent 0-30%; A_{J} > .3; .5 < genPt < 1");
  trkPtProjCanvas_p->Write("genPtJT_0_1_030_Asymm_c");
  claverCanvasSaving(trkPtProjCanvas_p, Form("../pngDir/genPtJT_0_1_030_Asymm_%s", fileTag), "png");

  outFile_p->Close();
  delete outFile_p;
  gStyle->SetOptStat(0);

  //  delete cent_p;
  //  delete asymmCanvas_p;
}






//Test


void makeAsymmImbProf(TTree* getTree_p, const char* outName, const char* gorr, const char* PFCalo, const char* perpProj, Int_t centLow, Int_t centHi, Int_t profLow, Int_t profHi, const char* GLN = "N", const char* Corr = "")
{
  if(!checkChar("Lead", "Pt", gorr, perpProj, "F", GLN)) return;

  const char* grpfcalo = Form("%s%s", gorr, PFCalo);

  inFile_p->cd();

  const char* title = Form("%sAsymmImb%s%s%s_%d%d_%s_%s", grpfcalo, perpProj, "F", Corr, (Int_t)(centLow*2.5), (Int_t)((centHi + 1)*2.5), GLN, fileTag);

  TString name = Form("%s_g", title);

  TGraphErrors* asymmImbGraph_p = new TGraphErrors(4);
  asymmImbGraph_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTGraphErrors(asymmImbGraph_p, profHi, profLow);
  TH1F* getHist_p;

  TString var = Form("%sImb%s%s%s", grpfcalo, perpProj, "F", Corr);

  TCut setCut = makeSetCut(gorr, PFCalo);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = "";
  TCut asymmCut = makeAsymmCut(grpfcalo, .00, .10);

  if(strcmp(gorr, "gR") == 0)    asymmCut = makeAsymmCut("g", .00, .10);

  //for backwards compat. temporary                                                                      
  if(strcmp(gorr, "r") == 0){
    if(*GLN == 71)
      etaCut = Form("TMath::Abs(rPFLeadJtEta) > 1.0 || TMath::Abs(rPFSubLeadJtEta) > 1.0");
    else if(*GLN == 76)
      etaCut = Form("TMath::Abs(rPFLeadJtEta) < 1.0 && TMath::Abs(rPFSubLeadJtEta) < 1.0");
  }
  else
    if(*GLN == 71)
      etaCut = Form("TMath::Abs(gLeadJtEta) > 1.0 || TMath::Abs(gSubLeadJtEta) > 1.0");
    else if(*GLN == 76)
      etaCut = Form("TMath::Abs(gLeadJtEta) < 1.0 && TMath::Abs(gSubLeadJtEta) < 1.0");

  getTree_p->Project(name, var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get(name);

  asymmImbGraph_p->SetPoint(0, 0.05, getHist_p->GetMean());
  asymmImbGraph_p->SetPointError(0, 0.05, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));

  asymmCut = makeAsymmCut(grpfcalo, .10, .20);
  getTree_p->Project(name, var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get(name);
  asymmImbGraph_p->SetPoint(1, .15, getHist_p->GetMean());
  asymmImbGraph_p->SetPointError(1, 0.05, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));

  asymmCut = makeAsymmCut(grpfcalo, .20, .30);
  getTree_p->Project(name, var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get(name);
  asymmImbGraph_p->SetPoint(2, .25, getHist_p->GetMean());
  asymmImbGraph_p->SetPointError(2, 0.05, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));

  asymmCut = makeAsymmCut(grpfcalo, .30, .50);
  getTree_p->Project(name, var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get(name);
  asymmImbGraph_p->SetPoint(3, .40, getHist_p->GetMean());
  asymmImbGraph_p->SetPointError(3, .10, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));


  asymmImbGraph_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTGraphErrors(asymmImbGraph_p, profHi, profLow);

  outFile_p = new TFile(outName, "UPDATE");
  asymmImbGraph_p->Write(Form("%s_g", title));
  outFile_p->Close();

  delete outFile_p;
  delete asymmImbGraph_p;
}

void makeAsymmImbProf_PT(TTree* getTree_p, const char*outName, const char* gorr, const char* PFCalo, const char* perpProj, const char* PT, Int_t centLow, Int_t centHi, Int_t profLow, Int_t profHi, const char* GLN = "N", const char* Corr = "")
{
  if(!checkChar("Lead", "Pt", gorr, perpProj, "F", GLN)) return;

  const char* grpfcalo = Form("%s%s", gorr, PFCalo);

  inFile_p->cd();

  const char* title = Form("%sAsymmImb%s%s%s_%d%d_%s_%s", grpfcalo, perpProj, Corr, PT, (Int_t)(centLow*2.5), (Int_t)((centHi + 1)*2.5), GLN, fileTag);

  TString name = Form("%s_g", title);

  TGraphErrors* asymmImbGraph_p = new TGraphErrors(4);
  asymmImbGraph_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTGraphErrors(asymmImbGraph_p, profHi, profLow);
  TH1F* getHist_p;

  TString var = Form("%sImb%s%s%s", grpfcalo, perpProj, Corr, PT);

  TCut setCut = makeSetCut(gorr, PFCalo);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = "";
  TCut asymmCut = makeAsymmCut(grpfcalo, .00, .10);

  if(strcmp(gorr, "gR") == 0)    asymmCut = makeAsymmCut("g", .00, .10);

  //for backwards compat. temporary                                                                                               
  if(strcmp(gorr, "r") == 0){
    if(*GLN == 71)
      etaCut = Form("TMath::Abs(rPFLeadJtEta) > 1.0 || TMath::Abs(rPFSubLeadJtEta) > 1.0");
    else if(*GLN == 76)
      etaCut = Form("TMath::Abs(rPFLeadJtEta) < 1.0 && TMath::Abs(rPFSubLeadJtEta) < 1.0");
  }
  else
    if(*GLN == 71)
      etaCut = Form("TMath::Abs(gLeadJtEta) > 1.0 || TMath::Abs(gSubLeadJtEta) > 1.0");
    else if(*GLN == 76)
      etaCut = Form("TMath::Abs(gLeadJtEta) < 1.0 && TMath::Abs(gSubLeadJtEta) < 1.0");

  getTree_p->Project(name, var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get(name);

  asymmImbGraph_p->SetPoint(0, 0.05, getHist_p->GetMean());
  asymmImbGraph_p->SetPointError(0, 0.05, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));

  asymmCut = makeAsymmCut(grpfcalo, .10, .20);
  getTree_p->Project(name, var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get(name);
  asymmImbGraph_p->SetPoint(1, .15, getHist_p->GetMean());
  asymmImbGraph_p->SetPointError(1, 0.05, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));

  asymmCut = makeAsymmCut(grpfcalo, .20, .30);
  getTree_p->Project(name, var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get(name);
  asymmImbGraph_p->SetPoint(2, .25, getHist_p->GetMean());
  asymmImbGraph_p->SetPointError(2, 0.05, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));

  asymmCut = makeAsymmCut(grpfcalo, .30, .50);
  getTree_p->Project(name, var, setCut && centCut && etaCut && asymmCut);
  getHist_p = (TH1F*)inFile_p->Get(name);
  asymmImbGraph_p->SetPoint(3, .40, getHist_p->GetMean());
  asymmImbGraph_p->SetPointError(3, .10, getHist_p->GetRMS()/(TMath::Sqrt(getHist_p->GetEntries())));

  asymmImbGraph_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTGraphErrors(asymmImbGraph_p, profHi, profLow);

  outFile_p = new TFile(outName, "UPDATE");
  asymmImbGraph_p->Write(Form("%s_g", title));
  outFile_p->Close();

  delete outFile_p;
  delete asymmImbGraph_p;

}


Double_t sumYForStack_PT(Double_t dIn = 0, Double_t comp1 = 0, Double_t comp2 = 0, Double_t comp3 = 0, Double_t comp4 = 0)
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


void makeHistForStack_PT(TGraph* g5_1_p, TGraph* g1_2_p, TGraph* g2_4_p, TGraph* g4_8_p, TGraph* g8_100_p, TGraph* gF_p, TH1F* h5_1_p, TH1F* h1_2_p, TH1F* h2_4_p, TH1F* h4_8_p, TH1F* h8_100_p, TH1F* hF_p)
{
  Double_t x5_1[4] = {0, 0, 0, 0};
  Double_t y5_1[4] = {0, 0, 0, 0};
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
    g5_1_p->GetPoint(iter, x5_1[iter], y5_1[iter]);
    g1_2_p->GetPoint(iter, x1_2[iter], y1_2[iter]);
    g2_4_p->GetPoint(iter, x2_4[iter], y2_4[iter]);
    g4_8_p->GetPoint(iter, x4_8[iter], y4_8[iter]);
    g8_100_p->GetPoint(iter, x8_100[iter], y8_100[iter]);
    gF_p->GetPoint(iter, xF[iter], yF[iter]);

    h8_100_p->SetBinContent(iter + 1, y8_100[iter]);
    h8_100_p->SetBinError(iter + 1, g8_100_p->GetErrorY(iter));

    h4_8_p->SetBinContent(iter + 1, sumYForStack_PT(y4_8[iter], y8_100[iter]));
    h4_8_p->SetBinError(iter + 1, g4_8_p->GetErrorY(iter));

    h2_4_p->SetBinContent(iter + 1, sumYForStack_PT(y2_4[iter], y4_8[iter], y8_100[iter]));
    h2_4_p->SetBinError(iter + 1, g2_4_p->GetErrorY(iter));

    h1_2_p->SetBinContent(iter + 1, sumYForStack_PT(y1_2[iter], y2_4[iter], y4_8[iter], y8_100[iter]));
    h1_2_p->SetBinError(iter + 1, g1_2_p->GetErrorY(iter));

    h5_1_p->SetBinContent(iter + 1, sumYForStack_PT(y5_1[iter], y1_2[iter], y2_4[iter], y4_8[iter], y8_100[iter]));
    h5_1_p->SetBinError(iter + 1, g5_1_p->GetErrorY(iter));

    hF_p->SetBinContent(iter + 1, yF[iter]);
    hF_p->SetBinError(iter + 1, gF_p->GetErrorY(iter));
  }

  h8_100_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  h8_100_p->SetXTitle("A_{J}");
  h4_8_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  h4_8_p->SetXTitle("A_{J}");
  h2_4_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  h2_4_p->SetXTitle("A_{J}");
  h1_2_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  h1_2_p->SetXTitle("A_{J}");
  h5_1_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  h5_1_p->SetXTitle("A_{J}");
  hF_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  hF_p->SetXTitle("A_{J}");

}

//Above, unchanged in case
/*
void makeHistForStack_PT(TGraph* g5_1_p, TGraph* g1_2_p, TGraph* g2_4_p, TGraph* g4_8_p, TGraph* g8_100_p, TGraph* gF_p, TH1F* h5_1_p, TH1F* h1_2_p, TH1F* h2_4_p, TH1F* h4_8_p, TH1F* h8_100_p, TH1F* hF_p)
{
  Double_t x5_1[4] = {0, 0, 0, 0};
  Double_t y5_1[4] = {0, 0, 0, 0};
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
    g5_1_p->GetPoint(iter, x5_1[iter], y5_1[iter]);
    g1_2_p->GetPoint(iter, x1_2[iter], y1_2[iter]);
    g2_4_p->GetPoint(iter, x2_4[iter], y2_4[iter]);
    g4_8_p->GetPoint(iter, x4_8[iter], y4_8[iter]);
    g8_100_p->GetPoint(iter, x8_100[iter], y8_100[iter]);
    gF_p->GetPoint(iter, xF[iter], yF[iter]);

    h8_100_p->SetBinContent(iter + 1, y8_100[iter]);
    h8_100_p->SetBinError(iter + 1, g8_100_p->GetErrorY(iter));

    h4_8_p->SetBinContent(iter + 1, y4_8[iter]);
    h4_8_p->SetBinError(iter + 1, g4_8_p->GetErrorY(iter));

    h2_4_p->SetBinContent(iter + 1, y4_8[iter] + y2_4[iter]);
    h2_4_p->SetBinError(iter + 1, g2_4_p->GetErrorY(iter));

    h1_2_p->SetBinContent(iter + 1, y4_8[iter] + y2_4[iter] + y1_2[iter]);
    h1_2_p->SetBinError(iter + 1, g1_2_p->GetErrorY(iter));

    h5_1_p->SetBinContent(iter + 1, y4_8[iter] + y2_4[iter] + y1_2[iter] + y5_1[iter]);
    h5_1_p->SetBinError(iter + 1, g5_1_p->GetErrorY(iter));

    hF_p->SetBinContent(iter + 1, yF[iter]);
    hF_p->SetBinError(iter + 1, gF_p->GetErrorY(iter));
  }
}
*/

void drawHistToStack_PT(TH1F* drawHist_p, Int_t color, const char* drawOpt)
{
  drawHist_p->SetFillColor(color);
  drawHist_p->SetMarkerStyle(6);
  drawHist_p->SetMarkerSize(.5);
  drawHist_p->Draw(drawOpt);
  drawHist_p->Draw("E1 SAME");
}


void makeAsymmImbStack_PT(const char* fileName, const char* gorr, const char* PFCalo, const char* perpProj, const char* GLN, const char* Corr = "")
{
  if(!checkChar("Lead", "Pt", gorr, perpProj, "F", GLN)) return;

  const char* grpfcalo = Form("%s%s", gorr, PFCalo);

  TFile* panelFile_p = new TFile(fileName, "UPDATE");
  TGraphErrors* getGraphC5_1_p;
  TGraphErrors* getGraphC1_2_p;
  TGraphErrors* getGraphC2_4_p;
  TGraphErrors* getGraphC4_8_p;
  TGraphErrors* getGraphC8_100_p;  
  TGraphErrors* getGraphCF_p;  

  TGraphErrors* getGraphP5_1_p;
  TGraphErrors* getGraphP1_2_p;
  TGraphErrors* getGraphP2_4_p;
  TGraphErrors* getGraphP4_8_p;
  TGraphErrors* getGraphP8_100_p;  
  TGraphErrors* getGraphPF_p;  

  Float_t binArrayX[5] = {.00, .10, .20, .30, .50};

  TH1F* histC_5_1_p = new TH1F("histC_5_1_p", "histC_5_1_p", 4, binArrayX);
  TH1F* histC_1_2_p = new TH1F("histC_1_2_p", "histC_1_2_p", 4, binArrayX);
  TH1F* histC_2_4_p = new TH1F("histC_2_4_p", "histC_2_4_p", 4, binArrayX);
  TH1F* histC_4_8_p = new TH1F("histC_4_8_p", "histC_4_8_p", 4, binArrayX);
  TH1F* histC_8_100_p = new TH1F("histC_8_100_p", "histC_8_100_p", 4, binArrayX);
  TH1F* histC_F_p = new TH1F("histC_F_p", "histC_F_p", 4, binArrayX);

  TH1F* histP_5_1_p = new TH1F("histP_5_1_p", "histP_5_1_p", 4, binArrayX);
  TH1F* histP_1_2_p = new TH1F("histP_1_2_p", "histP_1_2_p", 4, binArrayX);
  TH1F* histP_2_4_p = new TH1F("histP_2_4_p", "histP_2_4_p", 4, binArrayX);
  TH1F* histP_4_8_p = new TH1F("histP_4_8_p", "histP_4_8_p", 4, binArrayX);
  TH1F* histP_8_100_p = new TH1F("histP_8_100_p", "histP_8_100_p", 4, binArrayX);
  TH1F* histP_F_p = new TH1F("histP_F_p", "histP_F_p", 4, binArrayX);

  niceTH1(histC_5_1_p, 60, -60, -505, 406);
  niceTH1(histC_1_2_p, 60, -60, -505, 406);
  niceTH1(histC_2_4_p, 60, -60, -505, 406);
  niceTH1(histC_4_8_p, 60, -60, -505, 406);
  niceTH1(histC_8_100_p, 60, -60, -505, 406);
  niceTH1(histC_F_p, 60, -60, -505, 406);

  niceTH1(histP_5_1_p, 60, -60, -505, 406);
  niceTH1(histP_1_2_p, 60, -60, -505, 406);
  niceTH1(histP_2_4_p, 60, -60, -505, 406);
  niceTH1(histP_4_8_p, 60, -60, -505, 406);
  niceTH1(histP_8_100_p, 60, -60, -505, 406);
  niceTH1(histP_F_p, 60, -60, -505, 406);
  
  getGraphC5_1_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%s%s%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, "5_1", 30, 100, GLN, fileTag));
  getGraphC1_2_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%s%s%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, "1_2", 30, 100, GLN, fileTag));
  getGraphC2_4_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%s%s%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, "2_4", 30, 100, GLN, fileTag));
  getGraphC4_8_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%s%s%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, "4_8", 30, 100, GLN, fileTag));
  getGraphC8_100_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%s%s%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, "8_100", 30, 100, GLN, fileTag));
  getGraphCF_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%sF%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, 30, 100, GLN, fileTag));

  makeHistForStack_PT(getGraphC5_1_p, getGraphC1_2_p, getGraphC2_4_p, getGraphC4_8_p, getGraphC8_100_p, getGraphCF_p, histC_5_1_p, histC_1_2_p, histC_2_4_p, histC_4_8_p, histC_8_100_p, histC_F_p);

  TCanvas* profPanel_p = new TCanvas(Form("%sAsymmImb%s%sPTStack_%s_%s_c", grpfcalo, perpProj, Corr, GLN, fileTag), Form("%sAsymmImb%s%sPTStack_%s_%s_c", grpfcalo, perpProj, Corr, GLN, fileTag), 1);
  profPanel_p->Divide(2, 1, 0, 0);

  TLegend* leg;
  if(strcmp(gorr, "g") == 0)
    leg = new TLegend(0.15, 0.75, 0.95, 0.95, Form("Truth #slash{p}_{T}^{||} v. Truth A_{J}, %s, Truth Set", fileTag));
  else if(strcmp(gorr, "gR") == 0){
    if(strcmp(Corr, "Corr") == 0)
      leg = new TLegend(0.15, 0.75, 0.95, 0.95, Form("Corr. Reco #slash{p}_{T}^{||} v. Truth A_{J}, %s, Truth Set", fileTag));
    else if(strcmp(Corr, "") == 0)
      leg = new TLegend(0.15, 0.75, 0.95, 0.95, Form("Uncorr. Reco #slash{p}_{T}^{||} v. Truth A_{J}, %s, Truth Set", fileTag));
  }
  else if(strcmp(gorr, "r") == 0){
    if(strcmp(Corr, "Corr") == 0)
      leg = new TLegend(0.15, 0.75, 0.95, 0.95, Form("Corr. Reco #slash{p}_{T}^{||} v. Reco A_{J}, %s, Reco Set", fileTag));
    else if(strcmp(Corr, "") == 0)
      leg = new TLegend(0.15, 0.75, 0.95, 0.95, Form("Uncorr. Reco #slash{p}_{T}^{||} v. Reco A_{J}, %s, Reco Set", fileTag));
  }
  else {
    std::cout << "Legend not init: return" << std::endl;
    return;
  }

  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(.04);
  leg->SetBorderSize(0);

  profPanel_p->cd(1);


  drawHistToStack_PT(histC_5_1_p, kBlue - 9, "E1 HIST");
  leg->AddEntry(histC_5_1_p, ".5 < p_{T} < 1", "F");
  drawHistToStack_PT(histC_1_2_p, kYellow - 9, "E1 HIST SAME");
  leg->AddEntry(histC_1_2_p, "1 < p_{T} < 2", "F");
  drawHistToStack_PT(histC_2_4_p, kOrange + 1, "E1 HIST SAME");
  leg->AddEntry(histC_2_4_p, "2 < p_{T} < 4", "F");
  drawHistToStack_PT(histC_4_8_p, kGreen + 3, "E1 HIST SAME");
  leg->AddEntry(histC_4_8_p, "4 < p_{T} < 8", "F");
  drawHistToStack_PT(histC_8_100_p, kRed + 1, "E1 HIST SAME");
  leg->AddEntry(histC_8_100_p, "8 < p_{T}", "F");

  histC_F_p->Draw("SAME E1");

  leg->Draw("SAME");
  TLine* zeroLine_p = new TLine(0., 0., 0.5, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw();

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->DrawLatex(.2, .3, "30-100%");

  profPanel_p->cd(2);

  //Problem here!!!!
  
  getGraphP5_1_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%s%s%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, "5_1", 0, 30, GLN, fileTag));
  getGraphP1_2_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%s%s%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, "1_2", 0, 30, GLN, fileTag));
  getGraphP2_4_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%s%s%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, "2_4", 0, 30, GLN, fileTag));
  getGraphP4_8_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%s%s%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, "4_8", 0, 30, GLN, fileTag));
  getGraphP8_100_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%s%s%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, "8_100", 0, 30, GLN, fileTag));
  getGraphPF_p = (TGraphErrors*)panelFile_p->Get(Form("%sAsymmImb%sF%s_%d%d_%s_%s_g", grpfcalo, perpProj, Corr, 0, 30, GLN, fileTag));

  makeHistForStack_PT(getGraphP5_1_p, getGraphP1_2_p, getGraphP2_4_p, getGraphP4_8_p, getGraphP8_100_p, getGraphPF_p, histP_5_1_p, histP_1_2_p, histP_2_4_p, histP_4_8_p, histP_8_100_p, histP_F_p);

  drawHistToStack_PT(histP_5_1_p, kBlue - 9, "E1 HIST");
  drawHistToStack_PT(histP_1_2_p, kYellow - 9, "E1 HIST SAME");
  drawHistToStack_PT(histP_2_4_p, kOrange + 1, "E1 HIST SAME");
  drawHistToStack_PT(histP_4_8_p, kGreen + 3, "E1 HIST SAME");
  drawHistToStack_PT(histP_8_100_p, kRed + 1, "E1 HIST SAME");  

  histP_F_p->Draw("SAME E1");

  zeroLine_p->Draw();
  label_p->DrawLatex(.2, .3, "0-30%");
  
  label_p->DrawLatex(.1, .92, "Lead Jet p_{T} > 120 GeV/c");
  label_p->DrawLatex(.1, .88, "Sublead Jet p_{T} > 50 GeV/c");
  label_p->DrawLatex(.1, .84, "Jet #Delta #phi > 7#pi/8");
  label_p->DrawLatex(.1, .80, "Lead/Sublead Jet abs(#eta) < 1.6");

  profPanel_p->Write();

  claverCanvasSaving(profPanel_p, Form("../pngDir/%sAsymmImb%s%sPTStack_%s_%s", grpfcalo, perpProj, Corr, GLN, fileTag), "png");

  delete profPanel_p;

  panelFile_p->Close();
  delete panelFile_p;
}


void cfmDiJetAsymm(const char* inName = "inFile_CFMHIST_BETA.root", bool montecarlo = 0, const char* outName = "outFile_CFMHIST_BETA.root")
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
  makeAsymmHist(inTree_p, outName, "r", "PF", 10, 0, 1, 0, 39);                                          
  makeAsymmHist(inTree_p, outName, "r", "PF", 10, 0, 1, 0, 3);                                           
  makeAsymmHist(inTree_p, outName, "r", "PF", 10, 0, 1, 4, 7);                                           
  makeAsymmHist(inTree_p, outName, "r", "PF", 10, 0, 1, 8, 11);                                          
  makeAsymmHist(inTree_p, outName, "r", "PF", 10, 0, 1, 12, 19);                                         
  makeAsymmHist(inTree_p, outName, "r", "PF", 10, 0, 1, 20, 39);                                         
                                                                                                         
  makeAsymmPanel(outName, "r", "PF");                                                                    
  */

  //  makeImbHist(inTree_p, outName, "r", "PF", "Proj", "F", 50, -250, 250);
  //  makeImbHist(inTree_p, outName, "r", "PF", "Proj", "H", 21, -210, 210);                             
  //  makeImbHist(inTree_p, outName, "r", "PF", "Proj", "L", 21, -210, 210);                             
  //  makeImbHist(inTree_p, outName, "r", "PF", "Perp", "F", 50, -250, 250);
  //  makeImbHist(inTree_p, outName, "r", "PF", "Perp", "H", 21, -210, 210);                             
  //  makeImbHist(inTree_p, outName, "r", "PF", "Perp", "L", 21, -210, 210);                             

  //Corrected pt for tracks                                                                              
  //  makeImbHist(inTree_p, outName, "r", "PF", "Proj", "F", 50, -250, 250, "Corr");
  //  makeImbHist(inTree_p, outName, "r", "PF", "Perp", "F", 50, -250, 250, "Corr");

  makePtProjHist(inTree_p, outName);

  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 0, 39, -40, 40, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 0, 11, -40, 40, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 12, 39, -40, 40, "N");

  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 0, 39, -40, 40, "N", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 0, 11, -40, 40, "N", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 12, 39, -40, 40, "N", "Corr");

  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 0, 39, -40, 40, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 0, 11, -40, 40, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 12, 39, -40, 40, "N");

  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 0, 39, -40, 40, "N", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 0, 11, -40, 40, "N", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 12, 39, -40, 40, "N", "Corr");

  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 0, 11, -40, 40, "G");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 12, 39, -40, 40, "G");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 0, 11, -40, 40, "G");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 12, 39, -40, 40, "G");

  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 0, 11, -40, 40, "L");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 12, 39, -40, 40, "L");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 0, 11, -40, 40, "L");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 12, 39, -40, 40, "L");

  //Corrected pt for tracks                                                                                
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 0, 11, -40, 40, "G", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 12, 39, -40, 40, "G", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 0, 11, -40, 40, "G", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 12, 39, -40, 40, "G", "Corr");

  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 0, 11, -40, 40, "L", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Proj", 12, 39, -40, 40, "L", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 0, 11, -40, 40, "L", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "PF", "Perp", 12, 39, -40, 40, "L", "Corr");


  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "5_1", 0, 11, -60, 60, "N");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "1_2", 0, 11, -60, 60, "N");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "2_4", 0, 11, -60, 60, "N");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "4_8", 0, 11, -60, 60, "N");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "8_100", 0, 11, -60, 60, "N");

  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "5_1", 12, 39, -60, 60, "N");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "1_2", 12, 39, -60, 60, "N");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "2_4", 12, 39, -60, 60, "N");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "4_8", 12, 39, -60, 60, "N");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "8_100", 12, 39, -60, 60, "N");

  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "5_1", 0, 11, -60, 60, "N", "Corr");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "1_2", 0, 11, -60, 60, "N", "Corr");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "2_4", 0, 11, -60, 60, "N", "Corr");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "4_8", 0, 11, -60, 60, "N", "Corr");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "8_100", 0, 11, -60, 60, "N", "Corr");

  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "5_1", 12, 39, -60, 60, "N", "Corr");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "1_2", 12, 39, -60, 60, "N", "Corr");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "2_4", 12, 39, -60, 60, "N", "Corr");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "4_8", 12, 39, -60, 60, "N", "Corr");
  makeAsymmImbProf_PT(inTree_p, outName, "r", "PF", "Proj", "8_100", 12, 39, -60, 60, "N", "Corr");

  makeAsymmImbStack_PT(outName, "r", "PF", "Proj", "N", "");
  makeAsymmImbStack_PT(outName, "r", "PF", "Proj", "N", "Corr");

  if(montecarlo){
    makeJtCompHist(inTree_p, outName, "Lead", "Phi", 32, -3.2, 3.2, 32, -3.2, 3.2);

    //The "g" unit

    /*                                                                                                                     
    makeAsymmHist(inTree_p, outName, "g", "", 10, 0, 1, 0, 39);                                                            
    makeAsymmHist(inTree_p, outName, "g", "", 10, 0, 1, 0, 3);                                                             
    makeAsymmHist(inTree_p, outName, "g", "", 10, 0, 1, 4, 7);                                                             
    makeAsymmHist(inTree_p, outName, "g", "", 10, 0, 1, 8, 11);                                                            
    makeAsymmHist(inTree_p, outName, "g", "", 10, 0, 1, 12, 19);                                                           
    makeAsymmHist(inTree_p, outName, "g", "", 10, 0, 1, 20, 39);                                                           
                                                                                                                           
    makeAsymmPanel(outName, "g", "");                                                                                      
    */

    //    makeImbHist(inTree_p, outName, "g", "", "Proj", "F", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "g", "", "Proj", "H", 21, -210, 210);                                             
    //    makeImbHist(inTree_p, outName, "g", "", "Proj", "L", 21, -210, 210);                                             

    //    makeImbHist(inTree_p, outName, "g", "", "Perp", "F", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "g", "", "Perp", "H", 21, -210, 210);                                             
    //    makeImbHist(inTree_p, outName, "g", "", "Perp", "L", 21, -210, 210);                                             

   
    makeAsymmImbProf(inTree_p, outName, "g", "", "Proj", 0, 39, -20, 20, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "", "Proj", 0, 11, -15, 15, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "", "Proj", 12, 39, -15, 15, "N");

    makeAsymmImbProf(inTree_p, outName, "g", "", "Perp", 0, 39, -20, 20, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "", "Perp", 0, 11, -20, 20, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "", "Perp", 12, 39, -20, 20, "N");

    makeAsymmImbProf(inTree_p, outName, "g", "", "Proj", 0, 11, -20, 20, "G");
    makeAsymmImbProf(inTree_p, outName, "g", "", "Proj", 12, 39, -20, 20, "G");
    makeAsymmImbProf(inTree_p, outName, "g", "", "Perp", 0, 11, -20, 20, "G");
    makeAsymmImbProf(inTree_p, outName, "g", "", "Perp", 12, 39, -20, 20, "G");

    makeAsymmImbProf(inTree_p, outName, "g", "", "Proj", 0, 11, -20, 20, "L");
    makeAsymmImbProf(inTree_p, outName, "g", "", "Proj", 12, 39, -20, 20, "L");
    makeAsymmImbProf(inTree_p, outName, "g", "", "Perp", 0, 11, -20, 20, "L");
    makeAsymmImbProf(inTree_p, outName, "g", "", "Perp", 12, 39, -20, 20, "L");


    makeAsymmImbProf_PT(inTree_p, outName, "g", "", "Proj", "5_1", 0, 11, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "g", "", "Proj", "1_2", 0, 11, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "g", "", "Proj", "2_4", 0, 11, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "g", "", "Proj", "4_8", 0, 11, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "g", "", "Proj", "8_100", 0, 11, -60, 60, "N", "");
    
    makeAsymmImbProf_PT(inTree_p, outName, "g", "", "Proj", "5_1", 12, 39, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "g", "", "Proj", "1_2", 12, 39, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "g", "", "Proj", "2_4", 12, 39, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "g", "", "Proj", "4_8", 12, 39, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "g", "", "Proj", "8_100", 12, 39, -60, 60, "N", "");   

  makeAsymmImbStack_PT(outName, "g", "", "Proj", "N", "");

    //gRPF

    /*
    makeAsymmHist(inTree_p, outName, "gR", "PF", 10, 0, 1, 0, 39);
    makeAsymmHist(inTree_p, outName, "gR", "PF", 10, 0, 1, 0, 3);
    makeAsymmHist(inTree_p, outName, "gR", "PF", 10, 0, 1, 7);
    makeAsymmHist(inTree_p, outName, "gR", "PF", 10, 0, 1, 8, 11);
    makeAsymmHist(inTree_p, outName, "gR", "PF", 10, 0, 1, 12, 19);
    makeAsymmHist(inTree_p, outName, "gR", "PF", 10, 0, 1, 20, 39);
    
    makeAsymmPanel(outName, "gR", "PF");
    */

    //    makeImbHist(inTree_p, outName, "gR", "PF", "Proj", "F", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "gR", "PF", "Proj", "H", 21, -210, 210);                                        
    //    makeImbHist(inTree_p, outName, "gR", "PF", "Proj", "L", 21, -210, 210);                                       
    //    makeImbHist(inTree_p, outName, "gR", "PF", "Perp", "F", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "gR", "PF", "Perp", "H", 21, -210, 210);                                     
    //    makeImbHist(inTree_p, outName, "gR", "PF", "Perp", "L", 21, -210, 210);                         


    //    makeImbHist(inTree_p, outName, "gR", "PF", "Proj", "F", 21, -210, 210, "Corr");
    //    makeImbHist(inTree_p, outName, "gR", "PF", "Perp", "F", 21, -210, 210, "Corr");

    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 0, 39, -20, 20, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 0, 11, -15, 15, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 12, 39, -15, 15, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 0, 39, -20, 20, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 0, 11, -20, 20, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 12, 39, -20, 20, "N");

    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 0, 39, -20, 20, "N", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 0, 11, -15, 15, "N", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 12, 39, -15, 15, "N", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 0, 39, -20, 20, "N", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 0, 11, -20, 20, "N", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 12, 39, -20, 20, "N", "Corr");

    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 0, 11, -20, 20, "G");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 12, 39, -20, 20, "G");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 0, 11, -20, 20, "G");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 12, 39, -20, 20, "G");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 0, 11, -20, 20, "L");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 12, 39, -20, 20, "L");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 0, 11, -20, 20, "L");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 12, 39, -20, 20, "L");

    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 0, 11, -20, 20, "G", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 12, 39, -20, 20, "G", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 0, 11, -20, 20, "G", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 12, 39, -20, 20, "G", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 0, 11, -20, 20, "L", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Proj", 12, 39, -20, 20, "L", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 0, 11, -20, 20, "L", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "PF", "Perp", 12, 39, -20, 20, "L", "Corr");

    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "5_1", 0, 11, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "1_2", 0, 11, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "2_4", 0, 11, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "4_8", 0, 11, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "8_100", 0, 11, -60, 60, "N", "");
    
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "5_1", 12, 39, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "1_2", 12, 39, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "2_4", 12, 39, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "4_8", 12, 39, -60, 60, "N", "");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "8_100", 12, 39, -60, 60, "N", "");   

    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "5_1", 0, 11, -60, 60, "N", "Corr");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "1_2", 0, 11, -60, 60, "N", "Corr");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "2_4", 0, 11, -60, 60, "N", "Corr");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "4_8", 0, 11, -60, 60, "N", "Corr");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "8_100", 0, 11, -60, 60, "N", "Corr");
    
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "5_1", 12, 39, -60, 60, "N", "Corr");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "1_2", 12, 39, -60, 60, "N", "Corr");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "2_4", 12, 39, -60, 60, "N", "Corr");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "4_8", 12, 39, -60, 60, "N", "Corr");
    makeAsymmImbProf_PT(inTree_p, outName, "gR", "PF", "Proj", "8_100", 12, 39, -60, 60, "N", "Corr");   

    makeAsymmImbStack_PT(outName, "gR", "PF", "Proj", "N", "");
    makeAsymmImbStack_PT(outName, "gR", "PF", "Proj", "N", "Corr");
  }

  inFile_p->Close();
  delete inFile_p;
}
