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

  mcchain->Add("/dune/data/users/calcuttj/pduneana_Prod4a_1GeV_5_14_21.root/pduneana/beamana");


  TChain *datachain = new TChain();
  datachain->Add("/dune/data/users/calcuttj/pduneana_Prod4_1GeV_5387_5_12_21.root/pduneana/beamana");

  Unfold uf(pi::nthinslices+2, -1, pi::nthinslices+1);

  anavar mcevt(mcchain);

  ThinSlice mcths;
  mcths.SetOutputFileName("mcprod4a.root");
  mcths.Run(mcevt, uf);

  /*anavar dataevt(datachain);

  ThinSlice dataths;
  dataths.SetOutputFileName("data.root");
  dataths.Run(dataevt, uf);

  ThinSlice cosmicsths;
  cosmicsths.SetSelectCosmics(true);
  cosmicsths.SetOutputFileName("cosmics.root");
  cosmicsths.Run(dataevt, uf);*/

  return 0;

}
