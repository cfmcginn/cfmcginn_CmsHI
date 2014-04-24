//=============================================                                                         // Author: Chris McGinn                                                                                 //                                                                                                      // DiJet Analysis Class (MC)                                                                            //                                                                                                     
//=============================================                                                                         

#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"

#include <iostream>

//Current # of correction histograms

const Int_t nHist = 14;

//PuPF File and Hist array

TFile* PuPFFile_p[nHist];
TProfile* PuPFcent_p[nHist];
TProfile2D* PuPFphiEta_p[nHist];
TProfile* PuPFpt_p[nHist];
TProfile* PuPFdelR_p[nHist];

//PuCalo File and Hist array

TFile* PuCaloFile_p[nHist];
TProfile* PuCalocent_p[nHist];
TProfile2D* PuCalophiEta_p[nHist];
TProfile* PuCalopt_p[nHist];
TProfile* PuCalodelR_p[nHist];

//Vs Calo File and Hist array

TFile* VsCaloFile_p[nHist];
TProfile* VsCalocent_p[nHist];
TProfile2D* VsCalophiEta_p[nHist];
TProfile* VsCalopt_p[nHist];
TProfile* VsCalodelR_p[nHist];


//Fake Vs Calo File and Hist array

TFile* FakeVsCaloFile_p[nHist];
TProfile* FakeVsCalocent_p[nHist];
TProfile2D* FakeVsCalophiEta_p[nHist];
TProfile* FakeVsCalopt_p[nHist];
TProfile* FakeVsCalodelR_p[nHist];

//Fake Pu Calo File and Hist array

TFile* FakePuCaloFile_p[nHist];
TProfile* FakePuCalocent_p[nHist];
TProfile2D* FakePuCalophiEta_p[nHist];
TProfile* FakePuCalopt_p[nHist];
TProfile* FakePuCalodelR_p[nHist];

//Substitute in appropriate names, these are currently PuPF derived from 9k sample and Calo derived from 50k sample

void InitCorrFiles()
{
  //File names w/ various binnings, ordered by pt and then centrality. Each Jet Algorithm gets a file array
  
  PuPFFile_p[0] = new TFile("eff_pt0_1_cent0_10_step_cent4accept4pt4rmin3_PF.root", "READ");
  PuPFFile_p[1] = new TFile("eff_pt0_1_cent10_20_step_cent4accept4pt4rmin3_PF.root", "READ");
  PuPFFile_p[2] = new TFile("eff_pt0_1_cent20_30_step_cent4accept4pt4rmin3_PF.root", "READ");
  PuPFFile_p[3] = new TFile("eff_pt0_1_cent30_50_step_cent4accept4pt4rmin3_PF.root", "READ");
  PuPFFile_p[4] = new TFile("eff_pt0_1_cent50_100_step_cent4accept4pt4rmin3_PF.root", "READ");

  PuPFFile_p[5] = new TFile("eff_pt1_3_cent0_10_step_cent3accept3pt3rmin2_PF.root", "READ");
  PuPFFile_p[6] = new TFile("eff_pt1_3_cent10_20_step_cent3accept3pt3rmin2_PF.root", "READ");
  PuPFFile_p[7] = new TFile("eff_pt1_3_cent20_30_step_cent3accept3pt3rmin2_PF.root", "READ");
  PuPFFile_p[8] = new TFile("eff_pt1_3_cent30_50_step_cent3accept3pt3rmin2_PF.root", "READ");
  PuPFFile_p[9] = new TFile("eff_pt1_3_cent50_100_step_cent3accept3pt3rmin2_PF.root", "READ");

  PuPFFile_p[10] = new TFile("eff_pt3_8_cent0_10_step_cent3accept3pt3rmin2_PF.root", "READ");
  PuPFFile_p[11] = new TFile("eff_pt3_8_cent10_20_step_cent3accept3pt3rmin2_PF.root", "READ");
  PuPFFile_p[12] = new TFile("eff_pt3_8_cent20_100_step_cent3accept3pt3rmin2_PF.root", "READ");
  PuPFFile_p[13] = new TFile("eff_pt8_300_cent0_100_step_cent3accept3pt3rmin2_PF.root", "READ");


  PuCaloFile_p[0] = new TFile("eff_pt0_1_cent0_10_step_cent4accept4pt4rmin3_PuCalo.root", "READ");
  PuCaloFile_p[1] = new TFile("eff_pt0_1_cent10_20_step_cent4accept4pt4rmin3_PuCalo.root", "READ");
  PuCaloFile_p[2] = new TFile("eff_pt0_1_cent20_30_step_cent4accept4pt4rmin3_PuCalo.root", "READ");
  PuCaloFile_p[3] = new TFile("eff_pt0_1_cent30_50_step_cent4accept4pt4rmin3_PuCalo.root", "READ");
  PuCaloFile_p[4] = new TFile("eff_pt0_1_cent50_100_step_cent4accept4pt4rmin3_PuCalo.root", "READ");

  PuCaloFile_p[5] = new TFile("eff_pt1_3_cent0_10_step_cent4accept4pt4rmin3_PuCalo.root", "READ");
  PuCaloFile_p[6] = new TFile("eff_pt1_3_cent10_20_step_cent4accept4pt4rmin3_PuCalo.root", "READ");
  PuCaloFile_p[7] = new TFile("eff_pt1_3_cent20_30_step_cent4accept4pt4rmin3_PuCalo.root", "READ");
  PuCaloFile_p[8] = new TFile("eff_pt1_3_cent30_50_step_cent4accept4pt4rmin3_PuCalo.root", "READ");
  PuCaloFile_p[9] = new TFile("eff_pt1_3_cent50_100_step_cent4accept4pt4rmin3_PuCalo.root", "READ");

  PuCaloFile_p[10] = new TFile("eff_pt3_8_cent0_10_step_cent4accept4pt4rmin3_PuCalo.root", "READ");
  PuCaloFile_p[11] = new TFile("eff_pt3_8_cent10_20_step_cent4accept4pt4rmin3_PuCalo.root", "READ");
  PuCaloFile_p[12] = new TFile("eff_pt3_8_cent20_100_step_cent4accept4pt4rmin3_PuCalo.root", "READ");
  PuCaloFile_p[13] = new TFile("eff_pt8_300_cent0_100_step_cent3accept3pt3rmin3_PuCalo.root", "READ");
  
  VsCaloFile_p[0] = new TFile("eff_pt0_1_cent0_10_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[1] = new TFile("eff_pt0_1_cent10_20_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[2] = new TFile("eff_pt0_1_cent20_30_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[3] = new TFile("eff_pt0_1_cent30_50_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[4] = new TFile("eff_pt0_1_cent50_100_step_cent4accept4pt4rmin3_VsCalo.root", "READ");

  VsCaloFile_p[5] = new TFile("eff_pt1_3_cent0_10_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[6] = new TFile("eff_pt1_3_cent10_20_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[7] = new TFile("eff_pt1_3_cent20_30_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[8] = new TFile("eff_pt1_3_cent30_50_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[9] = new TFile("eff_pt1_3_cent50_100_step_cent4accept4pt4rmin3_VsCalo.root", "READ");

  VsCaloFile_p[10] = new TFile("eff_pt3_8_cent0_10_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[11] = new TFile("eff_pt3_8_cent10_20_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[12] = new TFile("eff_pt3_8_cent20_100_step_cent4accept4pt4rmin3_VsCalo.root", "READ");
  VsCaloFile_p[13] = new TFile("eff_pt8_300_cent0_100_step_cent3accept3pt3rmin3_VsCalo.root", "READ");


  //Fakes VsCalo

  FakeVsCaloFile_p[0] = new TFile("fake_pt0_1_cent0_10_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[1] = new TFile("fake_pt0_1_cent10_20_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[2] = new TFile("fake_pt0_1_cent20_30_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[3] = new TFile("fake_pt0_1_cent30_50_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[4] = new TFile("fake_pt0_1_cent50_100_akVs3Calo_dogenjet0.root", "READ");

  FakeVsCaloFile_p[5] = new TFile("fake_pt1_3_cent0_10_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[6] = new TFile("fake_pt1_3_cent10_20_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[7] = new TFile("fake_pt1_3_cent20_30_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[8] = new TFile("fake_pt1_3_cent30_50_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[9] = new TFile("fake_pt1_3_cent50_100_akVs3Calo_dogenjet0.root", "READ");

  FakeVsCaloFile_p[10] = new TFile("fake_pt3_8_cent0_10_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[11] = new TFile("fake_pt3_8_cent10_20_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[12] = new TFile("fake_pt3_8_cent20_100_akVs3Calo_dogenjet0.root", "READ");
  FakeVsCaloFile_p[13] = new TFile("fake_pt8_300_cent0_100_akVs3Calo_dogenjet0.root", "READ");

  //Fakes PuCalo
  
  FakePuCaloFile_p[0] = new TFile("fake_pt0_1_cent0_10_akPu3Calo_dogenjet0.root", "READ");
  FakePuCaloFile_p[1] = new TFile("fake_pt0_1_cent10_20_akPu3Calo_dogenjet0.root", "READ");
  FakePuCaloFile_p[2] = new TFile("fake_pt0_1_cent20_30_akPu3Calo_dogenjet0.root", "READ");
  FakePuCaloFile_p[3] = new TFile("fake_pt0_1_cent30_50_akPu3Calo_dogenjet0.root", "READ");
  FakePuCaloFile_p[4] = new TFile("fake_pt0_1_cent50_100_akPu3Calo_dogenjet0.root", "READ");

  FakePuCaloFile_p[5] = new TFile("fake_pt1_3_cent0_10_akPu3Calo_dogenjet0.root", "READ");
  FakePuCaloFile_p[6] = new TFile("fake_pt1_3_cent10_20_akPu3Calo_dogenjet0.root", "READ");
  FakePuCaloFile_p[7] = new TFile("fake_pt1_3_cent20_30_akPu3Calo_dogenjet0.root", "READ");
  FakePuCaloFile_p[8] = new TFile("fake_pt1_3_cent30_50_akPu3Calo_dogenjet0.root", "READ");
  FakePuCaloFile_p[9] = new TFile("fake_pt1_3_cent50_100_akPu3Calo_dogenjet0.root", "READ");

  FakePuCaloFile_p[10] = new TFile("fake_pt3_8_cent0_10_akPu3Calo_dogenjet0.root", "READ");
  FakePuCaloFile_p[11] = new TFile("fake_pt3_8_cent10_20_akPu3Calo_dogenjet0.root", "READ");
  FakePuCaloFile_p[12] = new TFile("fake_pt3_8_cent20_100_akPu3Calo_dogenjet0.root", "READ");
  FakePuCaloFile_p[13] = new TFile("fake_pt8_300_cent0_100_akPu3Calo_dogenjet0.root", "READ");
  
  return;
}


void InitCorrHists()
{
  for(Int_t hIter = 0; hIter < nHist; hIter++){
    
    PuPFcent_p[hIter] = (TProfile*)PuPFFile_p[hIter]->Get("p_eff_cent");
    PuPFphiEta_p[hIter] = (TProfile2D*)PuPFFile_p[hIter]->Get("p_eff_acceptance");
    PuPFpt_p[hIter] = (TProfile*)PuPFFile_p[hIter]->Get("p_eff_pt");
    PuPFdelR_p[hIter] = (TProfile*)PuPFFile_p[hIter]->Get("p_eff_rmin");

    PuCalocent_p[hIter] = (TProfile*)PuCaloFile_p[hIter]->Get("p_eff_cent");
    PuCalophiEta_p[hIter] = (TProfile2D*)PuCaloFile_p[hIter]->Get("p_eff_acceptance");
    PuCalopt_p[hIter] = (TProfile*)PuCaloFile_p[hIter]->Get("p_eff_pt");
    PuCalodelR_p[hIter] = (TProfile*)PuCaloFile_p[hIter]->Get("p_eff_rmin");
    

    VsCalocent_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_cent");
    VsCalophiEta_p[hIter] = (TProfile2D*)VsCaloFile_p[hIter]->Get("p_eff_acceptance");
    VsCalopt_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_pt");
    VsCalodelR_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_rmin");

    FakeVsCalocent_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_cent");
    FakeVsCalophiEta_p[hIter] = (TProfile2D*)FakeVsCaloFile_p[hIter]->Get("p_fake_acceptance");
    FakeVsCalopt_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_pt");
    FakeVsCalodelR_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_rmin");
    
    FakePuCalocent_p[hIter] = (TProfile*)FakePuCaloFile_p[hIter]->Get("p_fake_cent");
    FakePuCalophiEta_p[hIter] = (TProfile2D*)FakePuCaloFile_p[hIter]->Get("p_fake_acceptance");
    FakePuCalopt_p[hIter] = (TProfile*)FakePuCaloFile_p[hIter]->Get("p_fake_pt");
    FakePuCalodelR_p[hIter] = (TProfile*)FakePuCaloFile_p[hIter]->Get("p_fake_rmin");
    

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
  else if(8 <= pt)
    ptPos = hiSet4;

  return ptPos;
}


//Feed variables and the histograms along w/ the appropriate rmincut (currently rmin only defined to 3 for PuPF and 5 for calo)

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

  if(eff){
    corrFactor = corrFactor*(centProf_p->GetBinContent(centProf_p->FindBin(hiBin)));
    corrFactor = corrFactor*(etaPhiProf_p->GetBinContent(etaPhiProf_p->FindBin(phi, eta)));
    corrFactor = corrFactor*(ptProf_p->GetBinContent(ptProf_p->FindBin(pt)));
    if(rmin < rMinCut)
      corrFactor = corrFactor*(rminProf_p->GetBinContent(rminProf_p->FindBin(rmin)));
  }
  else{
    corrFactor = 0;
    corrFactor = corrFactor + centProf_p->GetBinContent(centProf_p->FindBin(hiBin));
    corrFactor = corrFactor + etaPhiProf_p->GetBinContent(etaPhiProf_p->FindBin(phi, eta));
    corrFactor = corrFactor + ptProf_p->GetBinContent(ptProf_p->FindBin(pt)) ;

    if(rmin < rMinCut)
      corrFactor = corrFactor +  rminProf_p->GetBinContent(rminProf_p->FindBin(rmin));
  }

  return corrFactor;
}



/*
Above implemented in code as follows:

 for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
      trkRMinPuPF_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[0]);
      trkRMinPuCalo_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[1]);
      trkRMinVsPF_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[2]);
      trkRMinVsCalo_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[3]);
    }

    if(montecarlo){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
        trkRMinT_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[3], true);
      }
    }

    Int_t hiBinDiv[5] = {20, 40, 60, 100, 200};
    Int_t hiSetEff[15] = {0, 5, 10, 1, 6, 11, 2, 7, 12, 3, 8, 12, 4, 9, 12};

    for(Int_t hiBinIter = 0; hiBinIter < 5; hiBinIter++){
      if(hiBin_ < hiBinDiv[hiBinIter]){
        for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
          Int_t ptPos = getPtBin(trkPt_[trkEntry], hiSetEff[hiBinIter*3], hiSetEff[hiBinIter*3 + 1], hiSetEff[hiBinIter*3 + 2], 13);

          trkPtFactPuPF_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPuPF_[trkEntry], PuPFcent_p[ptPos], PuPFphiEta_p[ptPos], PuPFpt_p[ptPos], PuPFdelR_p[ptPos], 3);
          trkPtCorrPuPF_[trkEntry] = trkPt_[trkEntry]/trkPtFactPuPF_[trkEntry];

          trkPtFactPuCalo_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPuCalo_[trkEntry], PuCalocent_p[ptPos], PuCalophiEta_p[ptPos], PuCalopt_p[ptPos], PuCalodelR_p[ptPos], 5);
          trkPtCorrPuCalo_[trkEntry] = trkPt_[trkEntry]*(1 - factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPuCalo_[trkEntry], FakePuCalocent_p[ptPos], FakePuCalophiEta_p[ptPos], FakePuCalopt_p[ptPos], FakePuCalodelR_p[ptPos], 5, false))/trkPtFactPuCalo_[trkEntry];

          trkPtFactVsCalo_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinVsCalo_[trkEntry], VsCalocent_p[ptPos], VsCalophiEta_p[ptPos], VsCalopt_p[ptPos], VsCalodelR_p[ptPos], 5);
          trkPtCorrVsCalo_[trkEntry] = trkPt_[trkEntry]*(1 - factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinVsCalo_[trkEntry], FakeVsCalocent_p[ptPos], FakeVsCalophiEta_p[ptPos], FakeVsCalopt_p[ptPos], FakeVsCalodelR_p[ptPos], 5, false))/trkPtFactVsCalo_[trkEntry];


	}
	break;
      }
    }
*/
