#include "HadAna.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TChain.h"
#include <iostream>

int main(){

  int nevents[nParTypes] = {0};

  TChain *chain = new TChain();

  chain->Add("/data/tjyang/dune/pduneana_Prod4_1GeV_2_9_21.root/pduneana/beamana");

  HadAna ana(chain);
  ana.AddTruePDG(-13);
  ana.AddTruePDG(211);
  ana.SetPandoraSlicePDG(13);

  Long64_t nentries = ana.fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = ana.LoadTree(jentry);
    if (ientry < 0) break;
    nb = ana.fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //std::cout<<ana.run<<" "<<ana.event<<" "<<ana.MC<<" "<<ana.reco_beam_true_byE_matched<<" "<<ana.true_beam_PDG<<" "<<(*ana.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    if (!ana.isTrueSelectedPart()) continue;
    int partype = ana.GetParType();
    if (partype<nParTypes){
      ++nevents[partype];
    }
  }

  for (int i = 0; i<nParTypes; ++i){
    std::cout<<i<<" "<<nevents[i]<<std::endl;
  }

  return 0;

}
