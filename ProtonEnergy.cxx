#include "ProtonEnergy.h"
#include "HadAna.h"
#include "anavar.h"
#include "TH1D.h"
#include "util.h"
#include <iostream>

using namespace std;

ProtonEnergy::ProtonEnergy(){
  hadana.InitP();
}

void ProtonEnergy::BookHistograms(){

  outputFile = TFile::Open(fOutputFileName.c_str(), "recreate");

  for (int i = 0; i < p::nIntTypes+1; ++i){
    htrack_length_ratio[i] = new TH1D(Form("htrack_length_ratio_%d",i),Form("%s;track_length_ratio", p::intTypeName[i]), 100, 0, 2);
  }

}

void ProtonEnergy::ProcessEvent(const anavar & evt){
  hadana.ProcessEvent(evt);
}

void ProtonEnergy::FillHistograms(const anavar & evt){

  double beamKE = sqrt(pow(evt.beam_inst_P*1000, 2) + pow(938.272, 2)) - 938.272;
  double hypoth_length = hadana.map_BB[2212]->RangeFromKESpline(beamKE);
  //cout<<evt.beam_inst_P<<" "<<beamKE<<" "<<hypoth_length<<" "<<evt.reco_beam_alt_len<<endl;
  FillHistVec1D(htrack_length_ratio, evt.reco_beam_alt_len/hypoth_length, hadana.ptype);
}

void ProtonEnergy::SaveHistograms(){
  outputFile->cd();
  outputFile->Write();
}

void ProtonEnergy::Run(anavar & evt, Long64_t nentries=-1){

  BookHistograms();

  //Long64_t nentries = evt.fChain->GetEntriesFast();
  if (nentries == -1) nentries = evt.fChain->GetEntries();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%100000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //std::cout<<evt.run<<" "<<evt.event<<" "<<evt.MC<<" "<<evt.reco_beam_true_byE_matched<<" "<<evt.true_beam_PDG<<" "<<(*evt.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    if (!hadana.isSelectedPart(evt)) continue;
    ProcessEvent(evt);
    if (!hadana.PassPCuts(evt)) continue;
    FillHistograms(evt);
  }
  SaveHistograms();
}

