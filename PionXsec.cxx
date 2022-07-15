#include "PionXsec.h"
#include "HadAna.h"
#include "anavar.h"
#include "TH1D.h"
#include "util.h"
#include <iostream>

using namespace std;

PionXsec::PionXsec(){
  hadana.InitPi();
}

void PionXsec::BookHistograms(){
  //cout << "[PionXsec::BookHistograms] Start" << endl;
  outputFile = TFile::Open(fOutputFileName.c_str(), "recreate");

  for (int i = 0; i < pi::nIntTypes+1; ++i){
    htrack_BeamKE[i] = new TH1D(Form("htrack_BeamKE_%d",i),Form("%s;track_BeamKE", pi::intTypeName[i]), 2000, 0, 2000);
    htrack_BeamP[i] = new TH1D(Form("htrack_BeamP_%d",i),Form("%s;track_BeamP", pi::intTypeName[i]), 2000, 0, 2000);
    htrack_fittedP[i] = new TH1D(Form("htrack_fittedP_%d",i),Form("%s;track_fittedP", pi::intTypeName[i]), 2000, 0, 2000);
    htrack_fitted_dP[i] = new TH1D(Form("htrack_fitted_dP_%d",i),Form("%s;track_fitted_dP", pi::intTypeName[i]), 2000, -1000, 2000);
    htrack_fittedKE[i] = new TH1D(Form("htrack_fittedKE_%d",i),Form("%s;track_fittedKE", pi::intTypeName[i]), 500, 0, 2000);
    htrack_fitted_dKE[i] = new TH1D(Form("htrack_fitted_dKE_%d",i),Form("%s;track_fitted_dKE", pi::intTypeName[i]), 2000, -1000, 1000);
    htrack_KECalo[i] = new TH1D(Form("htrack_KECalo_%d",i),Form("%s;track_KECalo", pi::intTypeName[i]), 500, 0, 2000);
    htrack_dKE_fitted_vs_KECalo[i] = new TH1D(Form("htrack_dKE_fitted_vs_KECalo_%d",i),Form("%s;track_dKE_fitted_vs_KECalo", pi::intTypeName[i]), 2000, -1000, 1000);
    htrack_KETruth[i] = new TH1D(Form("htrack_KETruth_%d",i),Form("%s;track_KETruth", pi::intTypeName[i]), 500, 0, 2000);
    htrack_dKE_fitted_vs_Truth[i] = new TH1D(Form("htrack_dKE_fitted_vs_Truth_%d",i),Form("%s;track_dKE_fitted_vs_Truth", pi::intTypeName[i]), 2000, -1000, 1000);
    hend_energy[i] = new TH1D(Form("hend_energy_%d",i),Form("%s;End point energy", pi::intTypeName[i]), 100, -500, 500);

    htrack_dKE_fitted_vs_Truth_2D[i] = new TH2D(Form("htrack_dKE_fitted_vs_Truth_2D_%d",i),Form("%s;track_dKE_fitted_vs_Truth", pi::intTypeName[i]),600, 200., 800.,  2000, -1000, 1000);
  }

}

void PionXsec::ProcessEvent(const anavar & evt){
  hadana.ProcessEvent(evt);
}

void PionXsec::FillHistograms(const anavar & evt){
  //cout << "[PionXsec::FillHistograms] Start" <<endl;

  double beamP = evt.beam_inst_P*1000.;
  double beamKE = sqrt(pow(beamP, 2) + pow(139.57, 2)) - 139.57;
  double KE_calo = hadana.Integrate_dEdx(*evt.reco_beam_dEdX_SCE, *evt.reco_beam_TrkPitch_SCE);
  double end_energy = hadana.map_BB[211]->KEAtLength(beamKE, evt.reco_beam_alt_len);
  FillHistVec1D(htrack_BeamKE, beamKE, hadana.pitype);
  FillHistVec1D(htrack_BeamP, evt.beam_inst_P*1000., hadana.pitype);
  FillHistVec1D(htrack_KECalo, KE_calo, hadana.pitype);
  FillHistVec1D(hend_energy, end_energy, hadana.pitype);

  if((*evt.reco_beam_dEdX_SCE).size() > 20){
    //cout << "[PionXsec::FillHistograms] fitting" << endl;
    double fitted_length = hadana.Fit_dEdx_Residual_Length(evt, *evt.reco_beam_dEdX_SCE, *evt.reco_beam_resRange_SCE, 211, false);
    double fitted_KE = hadana.map_BB[211]->KEFromRangeSpline(fitted_length);
    double fitted_dKE = beamKE - fitted_KE;
    double dKE_calo_vs_fitted = KE_calo - fitted_KE;
    double fitted_P = hadana.map_BB[211]->KEtoMomentum(fitted_KE);
    double fitted_dP = evt.beam_inst_P*1000. - fitted_P;
    double true_KE_ff = hadana.true_ffKE;
    double dKE_true_KE_ff_vs_fitted = true_KE_ff - fitted_KE;
    FillHistVec1D(htrack_fittedP, fitted_P, hadana.pitype);
    FillHistVec1D(htrack_fitted_dP, fitted_dP, hadana.pitype);
    FillHistVec1D(htrack_fittedKE, fitted_KE, hadana.pitype);
    FillHistVec1D(htrack_fitted_dKE, fitted_dKE, hadana.pitype);
    FillHistVec1D(htrack_dKE_fitted_vs_KECalo, dKE_calo_vs_fitted, hadana.pitype);
    FillHistVec1D(htrack_KETruth, true_KE_ff, hadana.pitype);    
    FillHistVec1D(htrack_dKE_fitted_vs_Truth, dKE_true_KE_ff_vs_fitted, hadana.pitype);
    FillHistVec2D(htrack_dKE_fitted_vs_Truth_2D, true_KE_ff, dKE_true_KE_ff_vs_fitted, hadana.pitype);
  }
}

void PionXsec::SaveHistograms(){
  outputFile->cd();
  outputFile->Write();
}

void PionXsec::Run(anavar & evt, Long64_t nentries=-1){

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
  //Make_dEdx_Range_Profile();
}

