#include "ThinSlice.h"
#include "Unfold.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TChain.h"
#include <iostream>

int main(){


  TChain *mcchain = new TChain();

  //chain->Add("/data/tjyang/dune/pduneana_Prod4.1_5_11_21.root/pduneana/beamana");
  //chain->Add("/data/tjyang/dune/pduneana_Prod4_1GeV_5_8_21.root/pduneana/beamana");

  mcchain->Add("/data/tjyang/dune/pduneana_Prod4a_1GeV_5_14_21.root/pduneana/beamana");


  TChain *datachain = new TChain();
  datachain->Add("/data/tjyang/dune/pduneana_Prod4_1GeV_5387_5_12_21.root/pduneana/beamana");

  Unfold uf(nthinslices+2, -1, nthinslices+1);

  HadAna mcevt(mcchain);
  mcevt.AddTruePDG(-13);
  mcevt.AddTruePDG(211);
  mcevt.SetPandoraSlicePDG(13);

  ThinSlice mcths;
  mcths.SetOutputFileName("mcprod4a.root");
  mcths.Run(mcevt, uf);

  HadAna dataevt(datachain);
  dataevt.AddTruePDG(13);
  dataevt.AddTruePDG(211);
  dataevt.SetPandoraSlicePDG(13);

  ThinSlice dataths;
  dataths.SetOutputFileName("data.root");
  dataths.Run(dataevt, uf);


  return 0;

}
