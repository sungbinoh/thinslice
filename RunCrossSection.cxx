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

  //mcchain->Add("/pnfs/dune/scratch/users/yinrui/hhr_data_50/pduneana_data.root/pduneana/beamana"); //test
  mcchain->Add("/pnfs/dune/scratch/users/yinrui/MC_211121/pduneana_MC.root/pduneana/beamana");

  TChain *datachain = new TChain();
  //datachain->Add("/pnfs/dune/scratch/users/yinrui/pduneana_Prod4_emfind_resel_allTrack/9_26_21_5387/output_sce_1GeV/pduneana_data.root/pduneana/beamana");
  datachain->Add("/pnfs/dune/scratch/users/yinrui/data_211121/pduneana_data.root/pduneana/beamana");

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
