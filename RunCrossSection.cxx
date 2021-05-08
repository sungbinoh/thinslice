#include "ThinSlice.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TChain.h"
#include <iostream>

int main(){

  int nevents[nParTypes] = {0};

  TChain *chain = new TChain();

  chain->Add("/data/tjyang/dune/pduneana_Prod4.1_1GeV_3_26_21.root/pduneana/beamana");

  HadAna evt(chain);
  evt.AddTruePDG(-13);
  evt.AddTruePDG(211);
  evt.SetPandoraSlicePDG(13);

  ThinSlice ths;
  TFile f("hadana.root","recreate");
  ths.BookHistograms();

  //Long64_t nentries = evt.fChain->GetEntriesFast();
  Long64_t nentries = evt.fChain->GetEntries();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //std::cout<<evt.run<<" "<<evt.event<<" "<<evt.MC<<" "<<evt.reco_beam_true_byE_matched<<" "<<evt.true_beam_PDG<<" "<<(*evt.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    if (!evt.isTrueSelectedPart()) continue;
    evt.ProcessEvent();
    ths.ProcessEvent(evt);
    ths.FillHistograms(kNocut, evt);
    if (evt.PassPandoraSliceCut()){
      ths.FillHistograms(kPandoraSlice, evt);
      if (evt.PassBeamQualityCut()){
        ths.FillHistograms(kBeamQuality, evt);
        if (evt.PassAPA3Cut()){
          ths.FillHistograms(kAPA3, evt);
          if (evt.PassCaloSizeCut()){
            ths.FillHistograms(kCaloSize, evt);
            if (evt.PassMichelScoreCut()){
              ths.FillHistograms(kMichelScore, evt);
              if (evt.PassMediandEdxCut()){
                ths.FillHistograms(kMediandEdx, evt);
              }
            }
          }
        }
      }
    }
  }

  ths.CalcXS();

  f.Write();
  return 0;

}
