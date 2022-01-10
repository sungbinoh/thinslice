#include "EffEval.h"
#include "HadAna.h"
#include "anavar.h"
#include "TGraphErrors.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TEfficiency.h"
#include "util.h"

#include <iostream>

EffEval::EffEval(){
  hadana.InitPi();
}

void EffEval::BookHistograms(){

  outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); // why c_str?
  h_true_Ppi_all = new TH1D("h_true_Ppi_all","All;p_{#pi} (MeV/c);Events",120,0,1200);
  h_true_Ppi_all->Sumw2();
  h_true_Ppi_sel = new TH1D("h_true_Ppi_sel","Selected;p_{#pi} (MeV/c);Events",120,0,1200);
  h_true_Ppi_sel->Sumw2();

}

void EffEval::ProcessEvent(const anavar & evt){
  hadana.ProcessEvent(evt);
}

void EffEval::FillHistograms(const anavar & evt){
  for (size_t i = 0; i<evt.true_beam_daughter_ID->size(); ++i){
    if (std::abs((*evt.true_beam_daughter_PDG)[i]) == 211){
      h_true_Ppi_all->Fill((*evt.true_beam_daughter_startP)[i]*1000);
      bool foundreco = false;
      for (size_t j = 0; j<evt.reco_daughter_PFP_ID->size(); ++j){
        if ((*evt.reco_daughter_PFP_true_byHits_ID)[j] == (*evt.true_beam_daughter_ID)[i]){
          foundreco = true;
        }
      }
      if (foundreco) h_true_Ppi_sel->Fill((*evt.true_beam_daughter_startP)[i]*1000);

    }
  }
}

void EffEval::SaveHistograms(){
  outputFile->cd();
  outputFile->Write();
//  for (int i = 0; i<h_true_Ppi_sel->GetNbinsX()+2; ++i){
//    cout<<i<<" "<<h_true_Ppi_sel->GetBinContent(i)<<" "<<h_true_Ppi_all->GetBinContent(i)<<endl;
//    if (h_true_Ppi_all->GetBinContent(i)<h_true_Ppi_sel->GetBinContent(i)) cout<<"Inconsistent"<<endl;
//  }
  if (TEfficiency::CheckConsistency(*h_true_Ppi_sel, *h_true_Ppi_all)){
    h_Eff_Ppi = new TEfficiency(*h_true_Ppi_sel, *h_true_Ppi_all);
    h_Eff_Ppi->Write("h_Eff_Ppi");
  }
}

void EffEval::Run(anavar & evt, Long64_t nentries=-1){

  BookHistograms();

  //Long64_t nentries = evt.fChain->GetEntriesFast();
  if (nentries == -1) nentries = evt.fChain->GetEntries();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //std::cout<<evt.run<<" "<<evt.event<<" "<<evt.MC<<" "<<evt.reco_beam_true_byE_matched<<" "<<evt.true_beam_PDG<<" "<<(*evt.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    if (!hadana.isSelectedPart(evt)) continue;
    ProcessEvent(evt);
    if (!hadana.PassPiCuts(evt)) continue;
    FillHistograms(evt);
  }
  SaveHistograms();
}
