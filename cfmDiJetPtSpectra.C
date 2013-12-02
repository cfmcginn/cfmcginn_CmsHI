#include "commonUtility.h"
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"

TFile* inFile_p = 0;
TFile* outFile_p = 0;

TTree* inTree_p = 0;

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
  TH1I* ptHist_p = new TH1I(Form("%s_h", title), Form("%s_h", title), nBins, histLow, histHi);

  TCut centCut = Form("%sPt < %f && hiBin >= %d && hiBin <= %d && %f < %sLeadDelPhi && %sLeadDelPhi < %f", genTrk, histHi, centLow, centHi, delPhiLow, genTrk, genTrk, delPhiHi);

  getTree_p->Project(Form("%s_h", title), Form("%sPt", genTrk), centCut);
  handsomeTH1(ptHist_p);

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
  TH1I* rHist_p = (TH1I*)outFile_p->Get(Form("rPt_%d%d_%s_%s_h", centLow, centHi, lAAll, fileTag));
  TH1I* gHist_p = (TH1I*)outFile_p->Get(Form("gPt_%d%d_%s_%s_h", centLow, centHi, lAAll, fileTag));

  gHist_p->SetMarkerColor(kBlue);
  gHist_p->SetMarkerStyle(20);
  rHist_p->SetMarkerColor(kRed);
  rHist_p->SetMarkerStyle(20);

  gHist_p->Draw("P E1");
  rHist_p->Draw("SAME P E1");

  rOnGCanvas_p->SetLogy();


  TLegend *leg = new TLegend(0.55, 0.45, 0.9, 0.65, Form("%d-%d%%, %s, %s", centLow, centHi, lAAll, fileTag));
  leg->SetFillColor(0);
  leg->AddEntry(gHist_p, "Truth", "p");
  leg->AddEntry(rHist_p, "Tracks", "p");
  leg->Draw("SAME");

  rOnGCanvas_p->Write();
  claverCanvasSaving(rOnGCanvas_p, Form("../pngDir/rOnGPt_%d%d_%s_%s", centLow, centHi, lAAll, fileTag), "png");                    
  outFile_p->Close();
  delete leg;
  delete rOnGCanvas_p;
  delete outFile_p;
}


void makeRDivGHist(const char* fileName, Int_t centLow, Int_t centHi, const char* lAAll)
{
  outFile_p = new TFile(fileName, "UPDATE");

  TCanvas* divCanvas_p = new TCanvas(Form("rDivGPt_%d%d_%s_%s_c", centLow, centHi, lAAll, fileTag), Form("rDivGPt_%d%d_%s_%s_c", centLow, centHi, lAAll, fileTag), 1);

  TH1I* numHist_p = (TH1I*)outFile_p->Get(Form("rPt_%d%d_%s_%s_h", centLow, centHi, lAAll, fileTag));
  TH1I* denomHist_p = (TH1I*)outFile_p->Get(Form("gPt_%d%d_%s_%s_h", centLow, centHi, lAAll, fileTag));

  TGraphAsymmErrors* divHist_p = new TGraphAsymmErrors(numHist_p, denomHist_p, "cl=.683 b(1,1) mode");
  handsomeTGraph(divHist_p);
  divHist_p->SetMinimum(0.);
  divHist_p->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  divHist_p->GetXaxis()->CenterTitle();
  divHist_p->GetYaxis()->SetTitle("Efficiency (Track/Truth)");
  divHist_p->GetYaxis()->CenterTitle();
  divHist_p->Draw("SAME AP");

  TLatex* cent_p = new TLatex(.5, .4, Form("%d-%d%%, %s, %s", centLow, centHi, lAAll, fileTag));
  cent_p->SetNDC();
  cent_p->Draw();

  divCanvas_p->Write();
  claverCanvasSaving(divCanvas_p, Form("../pngDir/rDivGPt_%d%d_%s_%s", centLow, centHi, lAAll, fileTag), "png");                    
  divHist_p->Write(Form("rDivGPt_%d%d_%s_%s_h", centLow, centHi, lAAll, fileTag));
  outFile_p->Close();

  delete divHist_p;
  delete divCanvas_p;
  delete outFile_p;
}


void cfmDiJetPtSpectra(const char* inName = "inFile_CFMHIST_BETA.root", bool montecarlo = 0, const char* outName = "outFile_CFMHIST_BETA.root")
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

  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 39, 0., TMath::PiOver2());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 11, 0., TMath::PiOver2());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 12, 39, 0., TMath::PiOver2());

  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 39, TMath::PiOver2(), TMath::Pi());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 11, TMath::PiOver2(), TMath::Pi());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 12, 39, TMath::PiOver2(), TMath::Pi());

  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 39, 0., TMath::Pi());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 0, 11, 0., TMath::Pi());
  makePtSpectra(inTree_p, outName, "r", "trk", 20, -.5, 19.5, 12, 39, 0., TMath::Pi());

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
  }

  inFile_p->Close();
  delete inFile_p;
}
