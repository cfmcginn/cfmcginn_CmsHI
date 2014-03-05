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

const Int_t nHist = 14;

//PF File and Hist array

TFile* PFFile_p[nHist];
TProfile* PFcent_p[nHist];
TProfile2D* PFphiEta_p[nHist];
TProfile* PFpt_p[nHist];
TProfile* PFdelR_p[nHist];

//Calo File and Hist array

TFile* CaloFile_p[nHist];
TProfile* Calocent_p[nHist];
TProfile2D* CalophiEta_p[nHist];
TProfile* Calopt_p[nHist];
TProfile* CalodelR_p[nHist];

//Vs Calo File and Hist array

TFile* VsCaloFile_p[nHist];
TProfile* VsCalocent_p[nHist];
TProfile2D* VsCalophiEta_p[nHist];
TProfile* VsCalopt_p[nHist];
TProfile* VsCalodelR_p[nHist];


//Fake Vs Calo File and Hist array

TFile* FakeVsCaloFile_p[nHist + 1];
TProfile* FakeVsCalocent_p[nHist + 1];
TProfile2D* FakeVsCalophiEta_p[nHist + 1];
TProfile* FakeVsCalopt_p[nHist + 1];
TProfile* FakeVsCalodelR_p[nHist + 1];

//Substitute in appropriate names, these are currently PF derived from 9k sample and Calo derived from 50k sample

void InitCorrFiles()
{
  PFFile_p[0] = new TFile("eff_pt0_1_cent0_10_step_cent4accept4pt4rmin3_PF.root", "READ");
  PFFile_p[1] = new TFile("eff_pt0_1_cent10_20_step_cent4accept4pt4rmin3_PF.root", "READ");
  PFFile_p[2] = new TFile("eff_pt0_1_cent20_30_step_cent4accept4pt4rmin3_PF.root", "READ");
  PFFile_p[3] = new TFile("eff_pt0_1_cent30_50_step_cent4accept4pt4rmin3_PF.root", "READ");
  PFFile_p[4] = new TFile("eff_pt0_1_cent50_100_step_cent4accept4pt4rmin3_PF.root", "READ");

  PFFile_p[5] = new TFile("eff_pt1_3_cent0_10_step_cent3accept3pt3rmin2_PF.root", "READ");
  PFFile_p[6] = new TFile("eff_pt1_3_cent10_20_step_cent3accept3pt3rmin2_PF.root", "READ");
  PFFile_p[7] = new TFile("eff_pt1_3_cent20_30_step_cent3accept3pt3rmin2_PF.root", "READ");
  PFFile_p[8] = new TFile("eff_pt1_3_cent30_50_step_cent3accept3pt3rmin2_PF.root", "READ");
  PFFile_p[9] = new TFile("eff_pt1_3_cent50_100_step_cent3accept3pt3rmin2_PF.root", "READ");

  PFFile_p[10] = new TFile("eff_pt3_8_cent0_10_step_cent3accept3pt3rmin2_PF.root", "READ");
  PFFile_p[11] = new TFile("eff_pt3_8_cent10_20_step_cent3accept3pt3rmin2_PF.root", "READ");
  PFFile_p[12] = new TFile("eff_pt3_8_cent20_100_step_cent3accept3pt3rmin2_PF.root", "READ");
  PFFile_p[13] = new TFile("eff_pt8_300_cent0_100_step_cent3accept3pt3rmin2_PF.root", "READ");


  CaloFile_p[0] = new TFile("eff_pt0_1_cent0_10_step_cent4accept4pt4rmin3_Calo.root", "READ");
  CaloFile_p[1] = new TFile("eff_pt0_1_cent10_20_step_cent4accept4pt4rmin3_Calo.root", "READ");
  CaloFile_p[2] = new TFile("eff_pt0_1_cent20_30_step_cent4accept4pt4rmin3_Calo.root", "READ");
  CaloFile_p[3] = new TFile("eff_pt0_1_cent30_50_step_cent4accept4pt4rmin3_Calo.root", "READ");
  CaloFile_p[4] = new TFile("eff_pt0_1_cent50_100_step_cent4accept4pt4rmin3_Calo.root", "READ");

  CaloFile_p[5] = new TFile("eff_pt1_3_cent0_10_step_cent3accept3pt3rmin2_Calo.root", "READ");
  CaloFile_p[6] = new TFile("eff_pt1_3_cent10_20_step_cent3accept3pt3rmin2_Calo.root", "READ");
  CaloFile_p[7] = new TFile("eff_pt1_3_cent20_30_step_cent3accept3pt3rmin2_Calo.root", "READ");
  CaloFile_p[8] = new TFile("eff_pt1_3_cent30_50_step_cent3accept3pt3rmin2_Calo.root", "READ");
  CaloFile_p[9] = new TFile("eff_pt1_3_cent50_100_step_cent3accept3pt3rmin2_Calo.root", "READ");

  CaloFile_p[10] = new TFile("eff_pt3_8_cent0_10_step_cent3accept3pt3rmin2_Calo.root", "READ");
  CaloFile_p[11] = new TFile("eff_pt3_8_cent10_20_step_cent3accept3pt3rmin2_Calo.root", "READ");
  CaloFile_p[12] = new TFile("eff_pt3_8_cent20_100_step_cent3accept3pt3rmin2_Calo.root", "READ");
  CaloFile_p[13] = new TFile("eff_pt8_300_cent0_100_step_cent3accept3pt3rmin2_Calo.root", "READ");

  VsCaloFile_p[0] = new TFile("eff_pt0_1_cent0_10_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[1] = new TFile("eff_pt0_1_cent10_20_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[2] = new TFile("eff_pt0_1_cent20_30_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[3] = new TFile("eff_pt0_1_cent30_50_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[4] = new TFile("eff_pt0_1_cent50_100_step_cent4accept4pt4rmin3_VsCalo.root", "READ");

  VsCaloFile_p[5] = new TFile("eff_pt1_3_cent0_10_step_cent3accept3pt3rmin2_VsCalo.root", "READ");
  VsCaloFile_p[6] = new TFile("eff_pt1_3_cent10_20_step_cent3accept3pt3rmin2_VsCalo.root", "READ");
  VsCaloFile_p[7] = new TFile("eff_pt1_3_cent20_30_step_cent3accept3pt3rmin2_VsCalo.root", "READ");
  VsCaloFile_p[8] = new TFile("eff_pt1_3_cent30_50_step_cent3accept3pt3rmin2_VsCalo.root", "READ");
  VsCaloFile_p[9] = new TFile("eff_pt1_3_cent50_100_step_cent3accept3pt3rmin2_VsCalo.root", "READ");

  VsCaloFile_p[10] = new TFile("eff_pt3_8_cent0_10_step_cent3accept3pt3rmin2_VsCalo.root", "READ");
  VsCaloFile_p[11] = new TFile("eff_pt3_8_cent10_20_step_cent3accept3pt3rmin2_VsCalo.root", "READ");
  VsCaloFile_p[12] = new TFile("eff_pt3_8_cent20_100_step_cent3accept3pt3rmin2_VsCalo.root", "READ");
  VsCaloFile_p[13] = new TFile("eff_pt8_300_cent0_100_step_cent3accept3pt3rmin2_VsCalo.root", "READ");


  //Fakes separate because of binning for now

  FakeVsCaloFile_p[0] = new TFile("fake_pt0_1_cent0_20_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[1] = new TFile("fake_pt0_1_cent20_40_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[2] = new TFile("fake_pt0_1_cent40_60_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[3] = new TFile("fake_pt0_1_cent60_100_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[4] = new TFile("fake_pt0_1_cent100_200_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");

  FakeVsCaloFile_p[5] = new TFile("fake_pt1_3_cent0_20_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[6] = new TFile("fake_pt1_3_cent20_40_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[7] = new TFile("fake_pt1_3_cent40_60_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[8] = new TFile("fake_pt1_3_cent60_100_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[9] = new TFile("fake_pt1_3_cent100_200_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");

  FakeVsCaloFile_p[10] = new TFile("fake_pt3_8_cent0_20_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[11] = new TFile("fake_pt3_8_cent20_40_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[12] = new TFile("fake_pt3_8_cent40_100_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[13] = new TFile("fake_pt3_8_cent100_200_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");

  FakeVsCaloFile_p[14] = new TFile("fake_pt8_300_cent0_200_step_cent4accept4pt4rmin3_akVs3Calo_dogenjet0.root", "READ");

  return;
}


void InitCorrHists()
{
  for(Int_t hIter = 0; hIter < nHist; hIter++){
    PFcent_p[hIter] = (TProfile*)PFFile_p[hIter]->Get("p_eff_cent");
    PFphiEta_p[hIter] = (TProfile2D*)PFFile_p[hIter]->Get("p_eff_acceptance");
    PFpt_p[hIter] = (TProfile*)PFFile_p[hIter]->Get("p_eff_pt");
    PFdelR_p[hIter] = (TProfile*)PFFile_p[hIter]->Get("p_eff_rmin");

    Calocent_p[hIter] = (TProfile*)CaloFile_p[hIter]->Get("p_eff_cent");
    CalophiEta_p[hIter] = (TProfile2D*)CaloFile_p[hIter]->Get("p_eff_acceptance");
    Calopt_p[hIter] = (TProfile*)CaloFile_p[hIter]->Get("p_eff_pt");
    CalodelR_p[hIter] = (TProfile*)CaloFile_p[hIter]->Get("p_eff_rmin");

    VsCalocent_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_cent");
    VsCalophiEta_p[hIter] = (TProfile2D*)VsCaloFile_p[hIter]->Get("p_eff_acceptance");
    VsCalopt_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_pt");
    VsCalodelR_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_rmin");
  }


  for(Int_t hIter = 0; hIter < nHist + 1; hIter++){
    FakeVsCalocent_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_cent");
    FakeVsCalophiEta_p[hIter] = (TProfile2D*)FakeVsCaloFile_p[hIter]->Get("p_fake_acceptance");
    FakeVsCalopt_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_pt");
    FakeVsCalodelR_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_rmin");
  }

  return;
}


Int_t getPtBin(Float_t pt, Int_t hiSet1, Int_t hiSet2, Int_t hiSet3, Int_t hiSet4)
{
  Int_t ptPos = -1;

  if(.5 <= pt && pt < 1)
    ptPos = hiSet1;
  else if(1 <= pt && pt < 3)
    ptPos = hiSet2;
  else if(3 <= pt && pt < 8)
    ptPos = hiSet3;
  else if(8 <= pt && pt < 100)
    ptPos = hiSet4;

  return ptPos;
}


//Feed variables and the histograms along w/ the appropriate rmincut (currently rmin only defined to 3 for PF and 5 for calo)

Float_t factorizedPtCorr(Int_t hiBin, Float_t pt, Float_t phi, Float_t eta, Float_t rmin, TProfile* centProf_p, TProfile2D* etaPhiProf_p, TProfile* ptProf_p, TProfile* rminProf_p, Int_t rMinCut, Bool_t eff = true)
{
  Float_t corrFactor = 1;

  if(hiBin < 0 || hiBin > 200){
    if(eff)
      return 1;
    else
      return 0;
  }

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

  corrFactor = corrFactor*(centProf_p->GetBinContent(hiBin+1));
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
