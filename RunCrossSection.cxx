#include "ThinSlice.h"
#include "Unfold.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TChain.h"
#include <iostream>

int main(){

  int nevents[nParTypes] = {0};

  TChain *chain = new TChain();

  //chain->Add("/data/tjyang/dune/pduneana_Prod4.1_5_11_21.root/pduneana/beamana");
  chain->Add("/data/tjyang/dune/pduneana_Prod4_1GeV_5_8_21.root/pduneana/beamana");

  Unfold uf(nthinslices+2, -1, nthinslices+1);

  HadAna evt(chain);
  evt.AddTruePDG(-13);
  evt.AddTruePDG(211);
  evt.SetPandoraSlicePDG(13);

  ThinSlice ths;
  ths.SetOutputFileName("mcprod4.root");
  ths.Run(evt, uf);
  return 0;

}
