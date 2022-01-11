#include "EffEval.h"
#include "HadAna.h"
#include "anavar.h"
#include "TGraphErrors.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "util.h"

#include <iostream>

EffEval::EffEval(){
  hadana.InitPi();
}

void EffEval::BookHistograms(){

  outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); // why c_str?
  h_true_Ppi_all = new TH1D("h_true_Ppi_all","All;True p_{#pi} (MeV/c);Events",120,0,1200);
  h_true_Ppi_all->Sumw2();

  h_true_Ppi_sel = new TH1D("h_true_Ppi_sel","Selected;True p_{#pi} (MeV/c);Events",120,0,1200);
  h_true_Ppi_sel->Sumw2();
  h_reco_true_Ppi_sel = new TH2D("h_reco_true_Ppi_sel","Selected;True p_{#pi} (MeV/c); Reco/True -1", 120,0,1200,40,-1,1);
  h_res_Ppi_sel = new TH1D("h_res_Ppi_sel","Selected; Reco/True -1 (p_{#pi})", 40,-1,1);
  h_res_Ppi_sel->Sumw2();

  h_res_Ppi_michelscore = new TH2D("h_res_Ppi_michelscore",";Michel score;Reco/True -1 (p_{#pi})", 100,0,1,40,-1,1);

  h_true_Ppi_michel = new TH1D("h_true_Ppi_michel","With Michel score cut;True p_{#pi} (MeV/c);Events",120,0,1200);
  h_true_Ppi_michel->Sumw2();
  h_reco_true_Ppi_michel = new TH2D("h_reco_true_Ppi_michel","With Michel score cut;True p_{#pi} (MeV/c); Reco/True -1", 120,0,1200,40,-1,1);
  h_res_Ppi_michel = new TH1D("h_res_Ppi_michel","With Michel score cut; Reco/True -1 (p_{#pi})", 40,-1,1);
  h_res_Ppi_michel->Sumw2();
}

void EffEval::ProcessEvent(const anavar & evt){
  hadana.ProcessEvent(evt);
}

void EffEval::FillHistograms(const anavar & evt){
  for (size_t i = 0; i<evt.true_beam_daughter_ID->size(); ++i){
    if (std::abs((*evt.true_beam_daughter_PDG)[i]) == 211){
      double true_mom = (*evt.true_beam_daughter_startP)[i]*1000;
      h_true_Ppi_all->Fill(true_mom);
      bool foundreco = false;
      double besttrack = -1;
      double tracklength = 0;
      double michelscore = 0;
      for (size_t j = 0; j<evt.reco_daughter_PFP_ID->size(); ++j){
        if ((*evt.reco_daughter_PFP_true_byHits_ID)[j] == (*evt.true_beam_daughter_ID)[i]){
          foundreco = true;
          if ((*evt.reco_daughter_PFP_true_byHits_purity)[j]*(*evt.reco_daughter_PFP_true_byHits_completeness)[j]>besttrack){
            besttrack = (*evt.reco_daughter_PFP_true_byHits_purity)[j]*(*evt.reco_daughter_PFP_true_byHits_completeness)[j];
            tracklength = (*evt.reco_daughter_allTrack_alt_len)[j];
            if ((*evt.reco_daughter_allTrack_vertex_nHits)[j]){
              michelscore = (*evt.reco_daughter_allTrack_vertex_michel_score)[j]/(*evt.reco_daughter_allTrack_vertex_nHits)[j];
            }
          }
        }
      }
      if (foundreco){
        double reco_KE = GetPionKE(tracklength);
        double reco_mom = sqrt(pow(reco_KE+139.57,2)-pow(139.57,2));
        double res_mom = (reco_mom - true_mom)/true_mom;
        h_true_Ppi_sel->Fill(true_mom);
        h_reco_true_Ppi_sel->Fill(true_mom, res_mom);
        h_res_Ppi_sel->Fill(res_mom);
        if (michelscore<0) michelscore = 0;
        if (michelscore>=1) michelscore = 0.999999;
        h_res_Ppi_michelscore->Fill(michelscore, res_mom);
        if (michelscore>0.5){
          h_true_Ppi_michel->Fill(true_mom);
          h_reco_true_Ppi_michel->Fill(true_mom, res_mom);
          h_res_Ppi_michel->Fill(res_mom);
        }
      }
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

  if (TEfficiency::CheckConsistency(*h_true_Ppi_michel, *h_true_Ppi_all)){
    h_Eff_Ppi_michel = new TEfficiency(*h_true_Ppi_michel, *h_true_Ppi_all);
    h_Eff_Ppi_michel->Write("h_Eff_Ppi_michel");
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
