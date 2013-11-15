#include <string>

void runCfmDiJetHist()
{
  gROOT->ProcessLine(".L ../gammaJetAnalysis/commonUtility.h");

  gROOT->ProcessLine(".L cfmDiJetHist.C++");
    
  gROOT->ProcessLine("cfmDiJetHist(\"Pythia80_HydjetDrum_mix01_HiForest2_v20_CFMSKIM.root\", 1, \"Pythia80_HydjetDrum_mix01_HiForest2_v20_CFMHIST.root\")");
}
