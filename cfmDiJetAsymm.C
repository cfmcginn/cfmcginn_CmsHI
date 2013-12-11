#include "commonUtility.h"
#include "TTree.h"
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

  uglyTProf->SetMarkerColor(1);
  uglyTProf->SetMarkerSize(1);
  uglyTProf->SetMarkerStyle(20);
  uglyTProf->SetLineColor(1);
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


TCut makeJetCut(const char* gorr)
{
  TCut jetCut = "";
  if(checkGORR(gorr))
    jetCut =  Form("%sLeadJtPt > 120 && %sSubLeadJtPt > 50 && getAbsDphi(%sLeadJtPhi, %sSubLeadJtPhi) > 7*(TMath::Pi())/8 && TMath::Abs(%sLeadJtEta) < 2.4 && TMath::Abs(%sSubLeadJtEta) < 2.4 ", gorr, gorr, gorr, gorr, gorr, gorr);
  else
    std::cout << "Warning: Jet Cut empty; gorr incorrectly specified" << std::endl;

  //Backcheck
  if(strcmp(gorr, "gR") == 0){
    jetCut =  Form("%sLeadJtPt > 120 && %sSubLeadJtPt > 50 && getAbsDphi(%sLeadJtPhi, %sSubLeadJtPhi) > 7*(TMath::Pi())/8 && TMath::Abs(%sLeadJtEta) < 2.4 && TMath::Abs(%sSubLeadJtEta) < 2.4 ", "g", "g", "g", "g", "g", "g");
    std::cout << "gR sub made" << std::endl;
  }

  return jetCut;
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




void makeJtCompHist(TTree* getTree_p, const char* outName, const char* lors, const char* ptPhiEta, Int_t xBins, Float_t xLow, Float_t xHi, Int_t yBins, Float_t yLow, Float_t yHi)
{
  if(!checkChar(lors, ptPhiEta)) return;

  inFile_p->cd();

  const char* title = Form("rVGJt%s_%s_%s", ptPhiEta, lors, fileTag);

  TString name = Form("%s_h", title);
  TH2F* jtCompHist_p = new TH2F(name, name, xBins, xLow, xHi, yBins, yLow, yHi);
  TCut histCut = Form("r%sJt%s > %f && r%sJt%s < %f && g%sJt%s > %f && g%sJt%s < %f", lors, ptPhiEta, xLow, lors, ptPhiEta, xHi, lors, ptPhiEta, yLow, lors, ptPhiEta, yHi);
  getTree_p->Project(name, Form("r%sJt%s:g%sJt%s", lors, ptPhiEta, lors, ptPhiEta), histCut);

  outFile_p = new TFile(outName, "UPDATE");
  jtCompHist_p->Write();
  outFile_p->Close();
  delete outFile_p;
}


void makeAsymmHist(TTree* getTree_p, const char* outName, const char* gorr, Int_t nBins, Int_t histLow, Int_t histHi, Int_t centLow, Int_t centHi)
{
  if(!checkChar("Lead", "Pt", gorr)) return;

  inFile_p->cd();

  const char* title = Form("%sAsymm_%d%d_%s", gorr, (Int_t)(centLow*2.5), (Int_t)((centHi+1)*2.5), fileTag);

  //  TCanvas* asymmCanvas_p = new TCanvas(Form("%s_c", title), Form("%s_c", title), 1);
  TH1F* asymmHist_p;

  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);
  TCut jetCut = makeJetCut(gorr);
  TCut centCut = makeCentCut(centLow, centHi);
  getTree_p->Project(name, Form("(%sLeadJtPt - %sSubLeadJtPt)/(%sLeadJtPt + %sSubLeadJtPt)", gorr, gorr, gorr, gorr), centCut && jetCut);
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


void addHistToPanel(TFile* file_p, TH1F* hist_p, TCanvas* canv_p, const char* gorr, Int_t centLow, Int_t centHi, Int_t pos)
{
  if(!checkChar("Lead", "Pt" , gorr)) return;

  hist_p = (TH1F*)file_p->Get(Form("%sAsymm_%d%d_%s_h", gorr, centLow, centHi, fileTag));
  canv_p->cd(pos);
  hist_p->Draw();

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->DrawLatex(.6, .3, Form("%d-%d%%", centLow, centHi));

  if(strcmp(gorr, "g") == 0 && pos == 1){
    label_p->DrawLatex(.3, .85, Form("Truth DiJet Asymmetry, Truth Subset"));
    label_p->DrawLatex(.3, .8, Form("%s", fileTag));
  }
  else if(strcmp(gorr, "r") == 0 && pos == 1){
    label_p->DrawLatex(.3, .85, Form("Reco DiJet Asymmetry, Reco Subset"));
    label_p->DrawLatex(.3, .8, Form("%s", fileTag));
  }
  else if(strcmp(gorr, "gR") == 0 && pos == 1){
    label_p->DrawLatex(.3, .85, Form("Reco DiJet Asymmetry, Truth Subset"));
    label_p->DrawLatex(.3, .8, Form("%s", fileTag));
  }


  delete label_p;
}


void makeAsymmPanel(const char* fileName, const char* gorr)
{
  if(!checkChar("Lead", "Pt",gorr)) return;

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
  claverCanvasSaving(asymmPanel_p, Form("../pngDir/%sAsymmPanel_%s", gorr, fileTag), "png");                                                 
  panelFile_p->Close();
  delete panelFile_p;
  delete asymmPanel_p;
}


void makeImbHist(TTree* getTree_p, const char* outName, const char* gorr, const char* perpProj, const char* FHL, Int_t nBins, Int_t histLow, Int_t histHi, const char* Corr = "")
{
  if(!checkChar("Lead", "Pt", gorr, perpProj, FHL)) return;

  inFile_p->cd();

  const char* title = Form("%sImb%s%s%s_%s", gorr, perpProj, FHL, Corr, fileTag);
  TH1F* imbHist_p;

  TString var = Form("%sImb%s%s%s", gorr, perpProj, FHL, Corr);
  TString name = Form("%s_h(%d, %d, %d)", title, nBins, histLow, histHi);
  TCut jetCut = makeJetCut(gorr);

  getTree_p->Project(name, var, jetCut);
  imbHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  imbHist_p->SetYTitle("Events");
  imbHist_p->SetXTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  handsomeTH1(imbHist_p);

  outFile_p = new TFile(outName, "UPDATE");
  imbHist_p->Write();

  outFile_p->Close();

  delete outFile_p;
}


void makeAsymmImbProf(TTree* getTree_p, const char* outName, const char* gorr, const char* perpProj, const char* FHL, Int_t xBins, Int_t yBins, Int_t yLow, Int_t yHi, Int_t centLow, Int_t centHi, Int_t profLow, Int_t profHi, Int_t profXDiv, Int_t profYDiv, const char* GLN = "N", const char* Corr = "")
{
  if(!checkChar("Lead", "Pt", gorr, perpProj, FHL, GLN)) return;
 
  inFile_p->cd();
  
  const char* title = Form("%sAsymmImb%s%s%s_%d%d_%s_%s", gorr, perpProj, FHL, Corr, (Int_t)(centLow*2.5), (Int_t)((centHi + 1)*2.5), GLN, fileTag);

  TH2F* asymmImbHist_p;
  TProfile* asymmImbHistProf_p;

  TString var = Form("%sImb%s%s%s:(gLeadJtPt - gSubLeadJtPt)/(gLeadJtPt + gSubLeadJtPt)", gorr, perpProj, FHL, Corr);

  //for backwards compat. temporary
  if(strcmp(gorr, "r") == 0)
    var = Form("%sImb%s%s%s:(%sLeadJtPt - %sSubLeadJtPt)/(%sLeadJtPt + %sSubLeadJtPt)", gorr, perpProj, FHL, Corr, gorr, gorr, gorr, gorr);

  TString name = Form("%s_h(%d, 0., 0.5, %d, %d, %d)", title, xBins, yBins, yLow, yHi);

  TCut centCut = makeCentCut(centLow, centHi);
  TCut jetCut = makeJetCut(gorr);
  TCut etaCut = "";

  //for backwards compat. temporary
  if(strcmp(gorr, "r") == 0){
    if(*GLN == 71)
      etaCut = Form("TMath::Abs(rLeadJtEta) > 1.0 || TMath::Abs(rSubLeadJtEta) > 1.0");
    else if(*GLN == 76)
      etaCut = Form("TMath::Abs(rLeadJtEta) < 1.0 && TMath::Abs(rSubLeadJtEta) < 1.0");
  }
  else
    if(*GLN == 71)
      etaCut = Form("TMath::Abs(gLeadJtEta) > 1.0 || TMath::Abs(gSubLeadJtEta) > 1.0");
    else if(*GLN == 76)
      etaCut = Form("TMath::Abs(gLeadJtEta) < 1.0 && TMath::Abs(gSubLeadJtEta) < 1.0");

  getTree_p->Project(name, var, centCut && jetCut && etaCut);
  asymmImbHist_p = (TH2F*)inFile_p->Get(Form("%s_h", title));

  asymmImbHist_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  asymmImbHist_p->SetXTitle("A_{J}");
  handsomeTH2(asymmImbHist_p);

  asymmImbHistProf_p = (TProfile*)asymmImbHist_p->ProfileX(Form("%s_prof", title));

  niceTProf(asymmImbHistProf_p, profHi, profLow, profXDiv, profYDiv);

  outFile_p = new TFile(outName, "UPDATE");
  asymmImbHist_p->Write();
  asymmImbHistProf_p->Write();
  outFile_p->Close();

  delete outFile_p;
}


void addProfToPanel(TFile* file_p, TProfile* prof_p, TCanvas* canv_p, const char* gorr, const char* perpProj, const char* FHL, Int_t centLow, Int_t centHi, Int_t pos, const char* GLN, const char* Corr = "", Option_t* drawOpt = "", Int_t color = 1, Int_t style = 20, TLegend* leg = 0)
{
  if(!checkChar("Lead", "Pt", gorr, perpProj, FHL, GLN)) return;

  prof_p = (TProfile*)file_p->Get(Form("%sAsymmImb%s%s%s_%d%d_%s_%s_prof", gorr, perpProj, FHL, Corr, centLow, centHi, GLN, fileTag));
  canv_p->cd(pos);
  prof_p->SetMarkerColor(color);
  prof_p->SetMarkerStyle(style);
  prof_p->SetLineColor(color);
  prof_p->Draw(drawOpt);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->DrawLatex(.6, .3, Form("%d-%d%%", centLow, centHi));

  if(strcmp(gorr, "g") == 0 && pos == 1 && color == 1){
    label_p->DrawLatex(.2, .85, Form("Truth <#slash{p}_{T}^{||}> v. Truth A_{J}, Truth Subset"));
    label_p->DrawLatex(.2, .8, Form("%s, %s", perpProj, fileTag));
  }
  else if(strcmp(gorr, "r") == 0 && pos == 1 && color == 1){
    label_p->DrawLatex(.2, .85, Form("Track <#slash{p}_{T}^{||}> v. Reco A_{J}, Reco Subset"));
    label_p->DrawLatex(.2, .8, Form("%s, %s", perpProj, fileTag));
  }
  else if(strcmp(gorr, "gR") == 0 && pos == 1 && color == 1){
    label_p->DrawLatex(.2, .85, Form("Track <#slash{p}_{T}^{||}> v. Truth A_{J}, Truth Subset"));
    label_p->DrawLatex(.2, .8, Form("%s, %s", perpProj, fileTag));
  }
  

  if(*GLN == 71 && leg != 0) 
    leg->AddEntry(prof_p, "One Jet abs(#eta) > 1.0", "p");
  else if(*GLN == 76 && leg != 0)
    leg->AddEntry(prof_p, "Both Jet abs(#eta) < 1.0", "p");


  delete label_p;
}


void makeAsymmImbPanel(const char* fileName, const char* gorr, const char* perpProj, const char* FHL, const char* GLN, const char* Corr = "")
{
  if(!checkChar("Lead", "Pt", gorr, perpProj, FHL, GLN)) return;

  TFile* panelFile_p = new TFile(fileName, "UPDATE");
  TProfile* getProf_p;

  TCanvas* profPanel_p = new TCanvas(Form("%sAsymmImb%s%s%sPanel_%s_%s_c", gorr, perpProj, FHL, Corr, GLN, fileTag), Form("%sAsymmImb%s%s%sPanel_%s_%s_c", gorr, perpProj, FHL, Corr, GLN, fileTag), 1);
  profPanel_p->Divide(2, 1, 0, 0);

  addProfToPanel(panelFile_p, getProf_p, profPanel_p, gorr, perpProj, FHL, 30, 100, 1, GLN, Corr);
  TLine *zeroLine_p = new TLine(0., 0., .5, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw();
  addProfToPanel(panelFile_p, getProf_p, profPanel_p, gorr, perpProj, FHL, 0, 30, 2, GLN, Corr);
  zeroLine_p->Draw();

  profPanel_p->Write();
  if(*FHL == 70)
    claverCanvasSaving(profPanel_p, Form("../pngDir/%sAsymmImb%s%s%sPanel_%s_%s", gorr, perpProj, FHL, Corr, GLN, fileTag), "png");                      
    
  panelFile_p->Close();
  delete panelFile_p;
  delete profPanel_p;
  delete zeroLine_p;
}


void makeEtaAsymmImbPanel(const char* fileName, const char* gorr, const char* perpProj, const char* FHL, const char* Corr = "")
{
  if(!checkChar("Lead", "Pt", gorr, perpProj, FHL)) return;

  TFile* panelFile_p = new TFile(fileName, "UPDATE");
  TProfile* getProf_p;

  TCanvas* profPanel_p = new TCanvas(Form("%sAsymmImb%s%s%sPanel_B_%s_c", gorr, perpProj, FHL, Corr, fileTag) , Form("%sAsymmImb%s%s%sPanel_B_%s_c", gorr, perpProj, FHL, Corr, fileTag), 1);
  profPanel_p->Divide(2,1,0,0);

  TLegend* leg;
  if(strcmp(perpProj, "Proj") == 0){
       if(strcmp(gorr, "g") == 0)
	 leg = new TLegend(0.15, 0.75, 0.95, .95, Form("Truth #slash{p}_{T}^{||} v. Truth A_{J}, %s, Truth Subset", fileTag));
       else if(strcmp(gorr, "gR") == 0)
	 leg = new TLegend(0.15, 0.75, 0.95, .95, Form("Track #slash{p}_{T}^{||} v. Truth A_{J}, %s, Truth Subset", fileTag));
       else if(strcmp(gorr, "r") == 0)
	 leg = new TLegend(0.15, 0.75, 0.95, .95, Form("Track #slash{p}_{T}^{||} v. Reco A_{J}, %s, Reco Subset", fileTag));
  }
  else if(strcmp(perpProj, "Perp") == 0){
       if(strcmp(gorr, "g") == 0)
	 leg = new TLegend(0.15, 0.75, 0.95, .95, Form("Truth #slash{p}_{T}^{+} v. Truth A_{J}, %s, Truth Subset", fileTag));
       else if(strcmp(gorr, "gR") == 0)
	 leg = new TLegend(0.15, 0.75, 0.95, .95, Form("Track #slash{p}_{T}^{+} v. Truth A_{J}, %s, Truth Subset", fileTag));
       else if(strcmp(gorr, "r") == 0)
	 leg = new TLegend(0.15, 0.75, 0.95, .95, Form("Track #slash{p}_{T}^{+} v. Reco A_{J}, %s, Reco Subset", fileTag));
  }
  else{
    std::cout << "Neither 'g' or 'r' or 'gR' input." << std::endl;
    return;
  }

  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(.04);
  addProfToPanel(panelFile_p, getProf_p, profPanel_p, gorr, perpProj, FHL, 30, 100, 1, "G", Corr, "P E1", 2, 20, leg);
  addProfToPanel(panelFile_p, getProf_p, profPanel_p, gorr, perpProj, FHL, 30, 100, 1, "L", Corr, "P E1 SAME", 4, 20, leg); 
  leg->Draw("SAME");
  TLine *zeroLine_p = new TLine(0., 0., .5, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw();

  addProfToPanel(panelFile_p, getProf_p, profPanel_p, gorr, perpProj, FHL, 0, 30, 2, "G", Corr, "P E1", 2, 20);
  addProfToPanel(panelFile_p, getProf_p, profPanel_p, gorr, perpProj, FHL, 0, 30, 2, "L", Corr, "P E1 SAME", 4, 20);

  zeroLine_p->Draw();

  profPanel_p->Write();

  if(*FHL == 70)
    claverCanvasSaving(profPanel_p, Form("../pngDir/%sAsymmImb%s%s%sPanel_B_%s", gorr, perpProj, FHL, Corr, fileTag), "png");                      

  panelFile_p->Close();
  delete zeroLine_p;
  delete profPanel_p;
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

  std::cout << "Filetag is: " << fileTag << std::endl;

  inFile_p = new TFile(inName, "READ");
  inTree_p = (TTree*)inFile_p->Get("jetTree");
  inTree_p->AddFriend("trackTree");

  if(montecarlo)
    inTree_p->AddFriend("genTree");

  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 0, 39);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 0, 3);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 4, 7);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 8, 11);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 12, 19);
  makeAsymmHist(inTree_p, outName, "r", 10, 0, 1, 20, 39);

  makeAsymmPanel(outName, "r");

  makeImbHist(inTree_p, outName, "r", "Proj", "F", 50, -250, 250);
  //  makeImbHist(inTree_p, outName, "r", "Proj", "H", 21, -210, 210);
  //  makeImbHist(inTree_p, outName, "r", "Proj", "L", 21, -210, 210);
  makeImbHist(inTree_p, outName, "r", "Perp", "F", 50, -250, 250);
  //  makeImbHist(inTree_p, outName, "r", "Perp", "H", 21, -210, 210);
  //  makeImbHist(inTree_p, outName, "r", "Perp", "L", 21, -210, 210);

  //Corrected pt for tracks
  makeImbHist(inTree_p, outName, "r", "Proj", "F", 50, -250, 250, "Corr");
  makeImbHist(inTree_p, outName, "r", "Perp", "F", 50, -250, 250, "Corr");

  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
  /*
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "H", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "H", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "H", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "L", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "L", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "L", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
  */

  //Corrected pt for tracks
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N", "Corr");

  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
  /*
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "H", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "H", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "H", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "L", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "L", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "L", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
  */

  //Corrected pt for tracks
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N", "Corr");

  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "G");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "G");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "G");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "G");

  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "L");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "L");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "L");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "L");

  //Corrected pt for tracks
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "G", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "G", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "G", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "G", "Corr");

  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "L", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "Proj", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "L", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "L", "Corr");
  makeAsymmImbProf(inTree_p, outName, "r", "Perp", "F", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "L", "Corr");


  makeAsymmImbPanel(outName, "r", "Proj", "F", "N");
  //  makeAsymmImbPanel(outName, "r", "Proj", "H", "N");
  //  makeAsymmImbPanel(outName, "r", "Proj", "L", "N");

  makeAsymmImbPanel(outName, "r", "Proj", "F", "G");
  makeAsymmImbPanel(outName, "r", "Proj", "F", "L");

  makeAsymmImbPanel(outName, "r", "Proj", "F", "N", "Corr");
  makeAsymmImbPanel(outName, "r", "Proj", "F", "G", "Corr");
  makeAsymmImbPanel(outName, "r", "Proj", "F", "L", "Corr");

  makeAsymmImbPanel(outName, "r", "Perp", "F", "N");
  //  makeAsymmImbPanel(outName, "r", "Perp", "H", "N");
  //  makeAsymmImbPanel(outName, "r", "Perp", "L", "N");

  makeEtaAsymmImbPanel(outName, "r", "Proj", "F");
  makeEtaAsymmImbPanel(outName, "r", "Perp", "F");

  makeAsymmImbPanel(outName, "r", "Perp", "F", "N", "Corr");
  makeEtaAsymmImbPanel(outName, "r", "Proj", "F", "Corr");
  makeEtaAsymmImbPanel(outName, "r", "Perp", "F", "Corr");

  if(montecarlo){
    makeJtCompHist(inTree_p, outName, "Lead", "Phi", 32, -3.2, 3.2, 32, -3.2, 3.2);

    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 0, 39);
    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 0, 3);
    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 4, 7);
    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 8, 11);
    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 12, 19);
    makeAsymmHist(inTree_p, outName, "g", 10, 0, 1, 20, 39);

    makeAsymmPanel(outName, "g");

    makeImbHist(inTree_p, outName, "g", "Proj", "F", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "g", "Proj", "H", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "g", "Proj", "L", 21, -210, 210);

    makeImbHist(inTree_p, outName, "g", "Perp", "F", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "g", "Perp", "H", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "g", "Perp", "L", 21, -210, 210);

    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "F", 10, 21, -210, 210, 0, 39, -20, 20, 505, 402, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "N");
    /*
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "H", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "H", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "H", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "L", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "L", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "L", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
    */

    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "F", 10, 21, -210, 210, 0, 39, -20, 20, 505, 402, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "N");
    /*
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "H", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "H", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "H", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "L", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "L", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "L", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
    */

    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "G");
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "G");
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "G");
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "G");

    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "L");
    makeAsymmImbProf(inTree_p, outName, "g", "Proj", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "L");
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "L");
    makeAsymmImbProf(inTree_p, outName, "g", "Perp", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "L");

    makeAsymmImbPanel(outName, "g", "Proj", "F", "N");
    //    makeAsymmImbPanel(outName, "g", "Proj", "H", "N");
    //    makeAsymmImbPanel(outName, "g", "Proj", "L", "N");

    makeAsymmImbPanel(outName, "g", "Perp", "F", "N");
    //    makeAsymmImbPanel(outName, "g", "Perp", "H", "N");
    //    makeAsymmImbPanel(outName, "g", "Perp", "L", "N");

    makeAsymmImbPanel(outName, "g", "Proj", "F", "G");
    makeAsymmImbPanel(outName, "g", "Proj", "F", "L");

    makeEtaAsymmImbPanel(outName, "g", "Proj", "F");
    makeEtaAsymmImbPanel(outName, "g", "Perp", "F");




    makeAsymmHist(inTree_p, outName, "gR", 10, 0, 1, 0, 39);
    makeAsymmHist(inTree_p, outName, "gR", 10, 0, 1, 0, 3);
    makeAsymmHist(inTree_p, outName, "gR", 10, 0, 1, 4, 7);
    makeAsymmHist(inTree_p, outName, "gR", 10, 0, 1, 8, 11);
    makeAsymmHist(inTree_p, outName, "gR", 10, 0, 1, 12, 19);
    makeAsymmHist(inTree_p, outName, "gR", 10, 0, 1, 20, 39);

    makeAsymmPanel(outName, "gR");

    makeImbHist(inTree_p, outName, "gR", "Proj", "F", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "gR", "Proj", "H", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "gR", "Proj", "L", 21, -210, 210);

    makeImbHist(inTree_p, outName, "gR", "Perp", "F", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "gR", "Perp", "H", 21, -210, 210);
    //    makeImbHist(inTree_p, outName, "gR", "Perp", "L", 21, -210, 210);

    makeImbHist(inTree_p, outName, "gR", "Proj", "F", 21, -210, 210, "Corr");
    makeImbHist(inTree_p, outName, "gR", "Perp", "F", 21, -210, 210, "Corr");

    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 0, 39, -20, 20, 505, 402, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "N");
    /*
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "H", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "H", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "H", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "L", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "L", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "L", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
    */

    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 0, 39, -20, 20, 505, 402, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "N");
    /*
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "H", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "H", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "H", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "L", 10, 21, -210, 210, 0, 39, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "L", 10, 21, -210, 210, 0, 11, -40, 40, 505, 404, "N");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "L", 10, 21, -210, 210, 12, 39, -40, 40, 505, 404, "N");
    */

    //Pt corrected tracks
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 0, 39, -20, 20, 505, 402, "N", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "N", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "N", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 0, 39, -20, 20, 505, 402, "N", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "N", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "N", "Corr");


    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "G");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "G");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "G");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "G");

    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "L");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "L");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "L");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "L");

    //Pt corrected tracks
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "G", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "G", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "G", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "G", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "L", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Proj", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "L", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 0, 11, -20, 20, 505, 402, "L", "Corr");
    makeAsymmImbProf(inTree_p, outName, "gR", "Perp", "F", 10, 21, -210, 210, 12, 39, -20, 20, 505, 402, "L", "Corr");



    makeAsymmImbPanel(outName, "gR", "Proj", "F", "N");
    //    makeAsymmImbPanel(outName, "gR", "Proj", "H", "N");
    //    makeAsymmImbPanel(outName, "gR", "Proj", "L", "N");

    makeAsymmImbPanel(outName, "gR", "Perp", "F", "N");
    //    makeAsymmImbPanel(outName, "gR", "Perp", "H", "N");
    //    makeAsymmImbPanel(outName, "gR", "Perp", "L", "N");

    makeAsymmImbPanel(outName, "gR", "Proj", "F", "G");
    makeAsymmImbPanel(outName, "gR", "Proj", "F", "L");

    makeEtaAsymmImbPanel(outName, "gR", "Proj", "F");
    makeEtaAsymmImbPanel(outName, "gR", "Perp", "F");

    //Pt corrected Tracks
    makeAsymmImbPanel(outName, "gR", "Proj", "F", "N", "Corr");
    makeAsymmImbPanel(outName, "gR", "Perp", "F", "N", "Corr");
    makeAsymmImbPanel(outName, "gR", "Proj", "F", "G", "Corr");
    makeAsymmImbPanel(outName, "gR", "Proj", "F", "L", "Corr");
    makeEtaAsymmImbPanel(outName, "gR", "Proj", "F", "Corr");
    makeEtaAsymmImbPanel(outName, "gR", "Perp", "F", "Corr");

  }

  inFile_p->Close();
  delete inFile_p;
}
