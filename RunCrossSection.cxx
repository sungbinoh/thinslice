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

  //mcchain->Add("pduneana_mc.root/pduneana/beamana"); // test
  mcchain->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/archive/sam_managed_users/tjyang/data/e/1/4/6/f52abbdc-f3f1-4b1e-9b05-0e25fd4bd232-whole_mc.root/pduneana/beamana");

  TChain *datachain = new TChain();
  //datachain->Add("pduneana.root/pduneana/beamana"); // test
  datachain->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/archive/sam_managed_users/tjyang/data/9/f/5/1/102e32ac-8d8d-44c6-8df0-6a0eb76ba1b4-whole_data.root/pduneana/beamana");

  Unfold uf(pi::nthinslices+2, -1, pi::nthinslices+1);

  anavar mcevt(mcchain);

  ThinSlice mcths;
  mcths.SetOutputFileName("mcprod4a.root");
  mcths.Run(mcevt, uf);

  anavar dataevt(datachain);

  ThinSlice dataths;
  dataths.SetOutputFileName("data.root");
  dataths.Run(dataevt, uf);

  ThinSlice cosmicsths;
  cosmicsths.SetSelectCosmics(true);
  cosmicsths.SetOutputFileName("cosmics.root");
  cosmicsths.Run(dataevt, uf);

  return 0;

}
