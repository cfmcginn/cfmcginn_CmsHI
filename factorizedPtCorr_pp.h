//=============================================                                                                         
// Author: Chris McGinn                                                                                                 
//                                                                                                                      
// DiJet Analysis Class (MC)                                                                                            
//                                                                                                                      
//=============================================                                                                         

#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"

#include <iostream>

//Current # of correction histograms

const Int_t nHist = 4;

//PF File and Hist array


//Vs Calo File and Hist array

TFile* CaloFile_p[nHist];
TProfile* Calocent_p[nHist];
TProfile2D* CalophiEta_p[nHist];
TProfile* Calopt_p[nHist];
TProfile* CalodelR_p[nHist];


//Fake Vs Calo File and Hist array

TFile* FakeCaloFile_p[nHist];
TProfile* FakeCalocent_p[nHist];
TProfile2D* FakeCalophiEta_p[nHist];
TProfile* FakeCalopt_p[nHist];
TProfile* FakeCalodelR_p[nHist];

//Substitute in appropriate names, these are currently PF derived from 9k sample and Calo derived from 50k sample

void InitCorrFiles()
{
  //File names w/ various binnings, ordered by pt and then centrality. Each Jet Algorithm gets a file array

  CaloFile_p[0] = new TFile("eff_pt0_1_ak3Calo_dogenjet0.root", "READ");
  CaloFile_p[1] = new TFile("eff_pt1_3_ak3Calo_dogenjet0.root", "READ");
  CaloFile_p[2] = new TFile("eff_pt3_8_ak3Calo_dogenjet0.root", "READ");
  CaloFile_p[3] = new TFile("eff_pt8_300_ak3Calo_dogenjet0.root", "READ");

  //Fakes Calo

  FakeCaloFile_p[0] = new TFile("fake_pt0_1_ak3Calo_dogenjet0.root", "READ");
  FakeCaloFile_p[1] = new TFile("fake_pt1_3_ak3Calo_dogenjet0.root", "READ");
  FakeCaloFile_p[2] = new TFile("fake_pt3_8_ak3Calo_dogenjet0.root", "READ");
  FakeCaloFile_p[3] = new TFile("fake_pt8_300_ak3Calo_dogenjet0.root", "READ");

  return;
}


void InitCorrHists()
{
  for(Int_t hIter = 0; hIter < nHist; hIter++){
    CalophiEta_p[hIter] = (TProfile2D*)CaloFile_p[hIter]->Get("p_eff_acceptance");
    Calopt_p[hIter] = (TProfile*)CaloFile_p[hIter]->Get("p_eff_pt");
    CalodelR_p[hIter] = (TProfile*)CaloFile_p[hIter]->Get("p_eff_rmin");

    FakeCalophiEta_p[hIter] = (TProfile2D*)FakeCaloFile_p[hIter]->Get("p_fake_acceptance");
    FakeCalopt_p[hIter] = (TProfile*)FakeCaloFile_p[hIter]->Get("p_fake_pt");
    FakeCalodelR_p[hIter] = (TProfile*)FakeCaloFile_p[hIter]->Get("p_fake_rmin");
  }
  return;
}


Int_t getPtBin(Float_t pt)
{
  Int_t ptPos = -1;

  if(.5 <= pt && pt < 1)
    ptPos = 0;
  else if(1 <= pt && pt < 3)
    ptPos = 1;
  else if(3 <= pt && pt < 8)
    ptPos = 2;
  else if(8 <= pt)
    ptPos = 3;

  return ptPos;
}


//Feed variables and the histograms along w/ the appropriate rmincut (currently rmin only defined to 3 for PF and 5 for calo)

Float_t factorizedPtCorr_pp(Float_t pt, Float_t phi, Float_t eta, Float_t rmin, TProfile2D* etaPhiProf_p, TProfile* ptProf_p, TProfile* rminProf_p, Int_t rMinCut, Bool_t eff = true)
{
  Float_t corrFactor = 1;

  if(pt < .5){
    if(eff)
      return 1;
    else
      return 0;
  }

  if(pt > 100){
    if(eff)
      return .8;
    else 
      return 0;
  }

  corrFactor = corrFactor*(etaPhiProf_p->GetBinContent(etaPhiProf_p->FindBin(phi, eta)));
  corrFactor = corrFactor*(ptProf_p->GetBinContent(ptProf_p->FindBin(pt)));

  if(rmin < rMinCut)
    corrFactor = corrFactor*(rminProf_p->GetBinContent(rminProf_p->FindBin(rmin)));

  return corrFactor;
}



/*
Above implemented in code as follows:

    Int_t hiBinDiv[5] = {20, 40, 60, 100, 200};
    Int_t hiSet[15] = {0, 5, 10, 1, 6, 11, 2, 7, 12, 3, 8, 12, 4, 9, 12};

    for(Int_t hiBinIter = 0; hiBinIter < 5; hiBinIter++){
      if(hiBin_ < hiBinDiv[hiBinIter]){
        for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
          Int_t ptPos = getPtBin(trkPt_[trkEntry], hiSet[hiBinIter*3], hiSet[hiBinIter*3 + 1], hiSet[hiBinIter*3 + 2]);

          trkPtCorrPF_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent_p[ptPos], PFphiEta_p[ptPos], PFpt_p[ptPos], PFdelR_p[ptPos], 3);
          trkPtFactPF_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPF_[trkEntry], PFcent_p[ptPos], PFphiEta_p[ptPos], PFpt_p[ptPos], PFdelR_p[ptPos], 3);

          trkPtCorrCalo_[trkEntry] = trkPt_[trkEntry]/factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinCalo_[trkEntry], Calocent_p[ptPos], CalophiEta_p[ptPos], Calopt_p[ptPos], CalodelR_p[ptPos], 5);
          trkPtFactCalo_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinCalo_[trkEntry], Calocent_p[ptPos], CalophiEta_p[ptPos], Calopt_p[ptPos], CalodelR_p[ptPos], 5);
        }
        break;
      }
    }

*/
