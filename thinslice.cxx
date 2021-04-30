#include "anavar.h"
#include "TChain.h"
#include <iostream>

int main(){

  TChain *chain = new TChain();

  chain->Add("/data/tjyang/dune/pduneana_Prod4_1GeV_2_9_21.root/pduneana/beamana");

  anavar ana(chain);
  
  Long64_t nentries = ana.fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = ana.LoadTree(jentry);
    if (ientry < 0) break;
    nb = ana.fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    std::cout<<ana.run<<" "<<ana.event<<" "<<ana.MC<<std::endl;
  }

  return 0;

}
