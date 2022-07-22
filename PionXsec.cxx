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

  // == Histograms for normalization
  h_cutflow = new TH1D("Cutflow", "Cutflow", 20, 0., 20.);

  
  // == Histograms for beam
  for (int i = 0; i < pi::nIntTypes+1; ++i){
    // == Histograms for beam before PiCuts
    htrack_BeamKE_precut[i] = new TH1D(Form("htrack_BeamKE_precut_%d",i),Form("%s;track_BeamKE", pi::intTypeName[i]), 5000, 0, 5000);
    htrack_BeamP_precut[i] = new TH1D(Form("htrack_BeamP_precut_%d",i),Form("%s;track_BeamP", pi::intTypeName[i]), 5000, 0, 5000);
    htrack_PandoraSlice_precut[i] = new TH1D(Form("htrack_PandoraSlice_precut_%d",i),Form("%s;PandoraSlice", pi::intTypeName[i]), 2., 0., 2.);
    htrack_CaloSize_precut[i] = new TH1D(Form("htrack_CaloSize_precut_%d",i),Form("%s;CaloSize", pi::intTypeName[i]), 2., 0., 2.);
    htrack_beam_dx_precut[i] = new TH1D(Form("htrack_beam_dx_precut_%d",i),Form("%s;beam_dx", pi::intTypeName[i]), 200., -10., 10.);
    htrack_beam_dy_precut[i] = new TH1D(Form("htrack_beam_dy_precut_%d",i),Form("%s;beam_dy", pi::intTypeName[i]), 200., -10., 10.);
    htrack_beam_dz_precut[i] = new TH1D(Form("htrack_beam_dz_precut_%d",i),Form("%s;beam_dz", pi::intTypeName[i]), 200., -10., 10.);
    htrack_beam_dxy_precut[i] = new TH1D(Form("htrack_beam_dxy_precut_%d",i),Form("%s;beam_dxy", pi::intTypeName[i]), 200., -10., 10.);
    htrack_beam_costh_precut[i] = new TH1D(Form("htrack_beam_costh_precut_%d",i),Form("%s;beam_dxy", pi::intTypeName[i]), 200., -1., 1.);
    htrack_daughter_michel_score_precut[i] = new TH1D(Form("htrack_daughter_michel_score_precut_%d",i),Form("%s;daughter_michel_score", pi::intTypeName[i]), 100., 0., 1.);
    htrack_chi2_proton_precut[i] = new TH1D(Form("htrack_chi2_proton_precut_%d",i),Form("%s;chi2_proton", pi::intTypeName[i]), 1000., 0., 100.);
    htrack_KEffTruth_precut[i] = new TH1D(Form("htrack_KEffTruth_precut_%d",i),Form("%s;track_KEffTruth_precut", pi::intTypeName[i]), 5000, 0, 5000);
    htrack_PffTruth_precut[i] = new TH1D(Form("htrack_PffTruth_precut_%d",i),Form("%s;track_PffTruth_precut", pi::intTypeName[i]), 5000, 0, 5000);
    htrack_beam_inst_XY_precut[i] = new TH2D(Form("htrack_beam_inst_XY_precut_%d",i),Form("%s;beam_inst_XY", pi::intTypeName[i]), 400., -50., -10., 400., 400., 440.);

    // == Histograms for beam after PiCuts
    htrack_BeamKE[i] = new TH1D(Form("htrack_BeamKE_%d",i),Form("%s;track_BeamKE", pi::intTypeName[i]), 5000, 0, 5000);
    htrack_BeamP[i] = new TH1D(Form("htrack_BeamP_%d",i),Form("%s;track_BeamP", pi::intTypeName[i]), 5000, 0, 5000);
    htrack_PandoraSlice[i] = new TH1D(Form("htrack_PandoraSlice_%d",i),Form("%s;PandoraSlice", pi::intTypeName[i]), 2., 0., 2.);
    htrack_CaloSize[i] = new TH1D(Form("htrack_CaloSize_%d",i),Form("%s;CaloSize", pi::intTypeName[i]), 2., 0., 2.);
    htrack_beam_dx[i] = new TH1D(Form("htrack_beam_dx_%d",i),Form("%s;beam_dx", pi::intTypeName[i]), 200., -10., 10.);
    htrack_beam_dy[i] = new TH1D(Form("htrack_beam_dy_%d",i),Form("%s;beam_dy", pi::intTypeName[i]), 200., -10., 10.);
    htrack_beam_dz[i] = new TH1D(Form("htrack_beam_dz_%d",i),Form("%s;beam_dz", pi::intTypeName[i]), 200., -10., 10.);
    htrack_beam_dxy[i] = new TH1D(Form("htrack_beam_dxy_%d",i),Form("%s;beam_dxy", pi::intTypeName[i]), 200., -10., 10.);
    htrack_beam_costh[i] = new TH1D(Form("htrack_beam_costh_%d",i),Form("%s;beam_dxy", pi::intTypeName[i]), 200., -1., 1.);
    htrack_daughter_michel_score[i] = new TH1D(Form("htrack_daughter_michel_score_%d",i),Form("%s;daughter_michel_score", pi::intTypeName[i]), 100., 0., 1.);
    htrack_chi2_proton[i] = new TH1D(Form("htrack_chi2_proton_%d",i),Form("%s;chi2_proton", pi::intTypeName[i]), 1000., 0., 100.);
    htrack_KEffTruth[i] = new TH1D(Form("htrack_KEffTruth_%d",i),Form("%s;track_KEffTruth", pi::intTypeName[i]), 5000, 0, 5000);
    htrack_PffTruth[i] = new TH1D(Form("htrack_PffTruth_%d",i),Form("%s;track_PffTruth", pi::intTypeName[i]), 5000, 0, 5000);
    htrack_beam_inst_XY[i] = new TH2D(Form("htrack_beam_inst_XY_%d",i),Form("%s;beam_inst_XY", pi::intTypeName[i]), 400., -50., -10., 400., 400., 440.);

    // == More plots
    htrack_BeamP_true[i] = new TH1D(Form("htrack_BeamP_true_%d",i),Form("%s;track_BeamP_true", pi::intTypeName[i]), 5000, 0, 5000);
    htrack_BeamKE_loss[i] = new TH1D(Form("htrack_BeamKE_loss_%d",i),Form("%s;track_BeamKE_loss", pi::intTypeName[i]), 1000, -500, 500);
    htrack_BeamKE_loss_true_mass[i] = new TH1D(Form("htrack_BeamKE_loss_true_mass_%d",i),Form("%s;track_BeamKE_loss_true_mass", pi::intTypeName[i]), 1000, -500, 500);
    htrack_Beam_alt_length[i] = new TH1D(Form("htrack_Beam_alt_length_%d",i),Form("%s;track_Beam_Length", pi::intTypeName[i]), 1000, 0, 1000);
    htrack_fittedP[i] = new TH1D(Form("htrack_fittedP_%d",i),Form("%s;track_fittedP", pi::intTypeName[i]), 2000, 0, 2000);
    htrack_fitted_dP[i] = new TH1D(Form("htrack_fitted_dP_%d",i),Form("%s;track_fitted_dP", pi::intTypeName[i]), 2000, -1000, 2000);
    htrack_fittedKE[i] = new TH1D(Form("htrack_fittedKE_%d",i),Form("%s;track_fittedKE", pi::intTypeName[i]), 500, 0, 2000);
    htrack_fitted_dKE[i] = new TH1D(Form("htrack_fitted_dKE_%d",i),Form("%s;track_fitted_dKE", pi::intTypeName[i]), 2000, -1000, 1000);
    htrack_KECalo[i] = new TH1D(Form("htrack_KECalo_%d",i),Form("%s;track_KECalo", pi::intTypeName[i]), 500, 0, 2000);
    htrack_dKE_fitted_vs_KECalo[i] = new TH1D(Form("htrack_dKE_fitted_vs_KECalo_%d",i),Form("%s;track_dKE_fitted_vs_KECalo", pi::intTypeName[i]), 2000, -1000, 1000);
    htrack_dKE_fitted_vs_Truth[i] = new TH1D(Form("htrack_dKE_fitted_vs_Truth_%d",i),Form("%s;track_dKE_fitted_vs_Truth", pi::intTypeName[i]), 2000, -1000, 1000);
  
    hend_energy[i] = new TH1D(Form("hend_energy_%d",i),Form("%s;End point energy", pi::intTypeName[i]), 400, -2000, 2000);

    htrack_dKE_fitted_vs_Truth_2D[i] = new TH2D(Form("htrack_dKE_fitted_vs_Truth_2D_%d",i),Form("%s;track_dKE_fitted_vs_Truth", pi::intTypeName[i]),600, 200., 800.,  2000, -1000, 1000);
  }

  // == Histograms for daughters
  hdaughter_proton_trackScore = new TH1D("hdaughter_proton_trackScore", "Daughter proton track score", 100, 0., 1.);
  hdaughter_proton_emScore = new TH1D("hdaughter_proton_emScore", "Daughter proton em score", 100, 0., 1.);
  hdaughter_proton_chi2_proton = new TH1D("hdaughter_proton_chi2_proton", "Daughter proton chi2 proton", 1000, 0., 100.);
  hdaughter_pion_trackScore = new TH1D("hdaughter_pion_trackScore", "Daughter pion track score", 100, 0., 1.);
  hdaughter_pion_emScore = new TH1D("hdaughter_pion_emScore", "Daughter pion em score", 100, 0., 1.);
  hdaughter_pion_chi2_proton = new TH1D("hdaughter_pion_chi2_proton", "Daughter pion chi2 proton", 1000, 0., 100.);
  hdaughter_other_trackScore = new TH1D("hdaughter_other_trackScore", "Daughter other track score", 100, 0., 1.);
  hdaughter_other_emScore = new TH1D("hdaughter_other_emScore", "Daughter other em score", 100, 0., 1.);
  hdaughter_other_chi2_proton = new TH1D("hdaughter_other_chi2_proton", "Daughter other chi2 proton", 1000, 0., 100.);
}

void PionXsec::ProcessEvent(const anavar & evt){
  hadana.ProcessEvent(evt);
}

void PionXsec::FillHistograms_precut(const anavar & evt){
  double beamP_scale = 1.0;
  if(evt.MC) beamP_scale = 0.5;
  double beamP = evt.beam_inst_P*1000. * beamP_scale;
  double beamKE = sqrt(pow(beamP, 2) + pow(139.57, 2)) - 139.57;
  FillHistVec1D(htrack_BeamKE_precut, beamKE, hadana.pitype);
  FillHistVec1D(htrack_BeamP_precut, beamP, hadana.pitype);
  FillHistVec1D(htrack_PandoraSlice_precut, hadana.PassPandoraSliceCut(evt), hadana.pitype);
  FillHistVec1D(htrack_CaloSize_precut, hadana.PassCaloSizeCut(evt), hadana.pitype);
  FillHistVec1D(htrack_beam_dx_precut, hadana.beam_dx, hadana.pitype);
  FillHistVec1D(htrack_beam_dy_precut, hadana.beam_dy, hadana.pitype);
  FillHistVec1D(htrack_beam_dz_precut, hadana.beam_dz, hadana.pitype);
  FillHistVec1D(htrack_beam_dxy_precut, hadana.beam_dxy, hadana.pitype);
  FillHistVec1D(htrack_beam_costh_precut, hadana.beam_costh, hadana.pitype);
  FillHistVec1D(htrack_daughter_michel_score_precut, hadana.daughter_michel_score, hadana.pitype);
  FillHistVec1D(htrack_chi2_proton_precut, hadana.chi2_proton, hadana.pitype);
  FillHistVec2D(htrack_beam_inst_XY_precut, evt.beam_inst_X, evt.beam_inst_Y, hadana.pitype);

  // == True Beam info
  double beamP_true = evt.true_beam_startP * 1000.;
  double true_KE_ff = hadana.true_ffKE;
  double true_P_ff = hadana.map_BB[abs(evt.true_beam_PDG)] -> KEtoMomentum(true_KE_ff);
  FillHistVec1D(htrack_BeamP_true, beamP_true, hadana.pitype);
  FillHistVec1D(htrack_KEffTruth_precut, true_KE_ff, hadana.pitype);
  FillHistVec1D(htrack_PffTruth_precut, true_P_ff, hadana.pitype);
}

void PionXsec::FillHistograms(const anavar & evt){
  //cout << "[PionXsec::FillHistograms] Start" <<endl;
  double beamP_scale = 1.0;
  if(evt.MC) beamP_scale = 0.5;
  double beamP = evt.beam_inst_P*1000. * beamP_scale;
  double beamKE = sqrt(pow(beamP, 2) + pow(139.57, 2)) - 139.57;
  double beamKE_true_mass = beamKE;
  if(evt.MC) beamKE_true_mass = hadana.map_BB[abs(evt.true_beam_PDG)] -> MomentumtoKE(beamP);
  double KE_loss = beamKE - hadana.true_ffKE;
  double KE_loss_true_mass = beamKE_true_mass - hadana.true_ffKE;
  double KE_calo = hadana.Integrate_dEdx(*evt.reco_beam_dEdX_SCE, *evt.reco_beam_TrkPitch_SCE);
  double end_energy = hadana.map_BB[211]->KEAtLength(beamKE, evt.reco_beam_alt_len);
  FillHistVec1D(htrack_BeamKE, beamKE, hadana.pitype);
  FillHistVec1D(htrack_BeamP, beamP, hadana.pitype);
  FillHistVec1D(htrack_PandoraSlice, hadana.PassPandoraSliceCut(evt), hadana.pitype);
  FillHistVec1D(htrack_CaloSize, hadana.PassCaloSizeCut(evt), hadana.pitype);
  FillHistVec1D(htrack_beam_dx, hadana.beam_dx, hadana.pitype);
  FillHistVec1D(htrack_beam_dy, hadana.beam_dy, hadana.pitype);
  FillHistVec1D(htrack_beam_dz, hadana.beam_dz, hadana.pitype);
  FillHistVec1D(htrack_beam_dxy, hadana.beam_dxy, hadana.pitype);
  FillHistVec1D(htrack_beam_costh, hadana.beam_costh, hadana.pitype);
  FillHistVec1D(htrack_daughter_michel_score, hadana.daughter_michel_score, hadana.pitype);
  FillHistVec1D(htrack_chi2_proton, hadana.chi2_proton, hadana.pitype);
  FillHistVec2D(htrack_beam_inst_XY, evt.beam_inst_X, evt.beam_inst_Y, hadana.pitype);

  FillHistVec1D(htrack_BeamKE_loss, KE_loss, hadana.pitype);
  FillHistVec1D(htrack_BeamKE_loss_true_mass, KE_loss_true_mass, hadana.pitype);
  FillHistVec1D(htrack_Beam_alt_length, evt.reco_beam_alt_len, hadana.pitype);
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
    FillHistVec1D(htrack_KEffTruth, true_KE_ff, hadana.pitype);    
    FillHistVec1D(htrack_dKE_fitted_vs_Truth, dKE_true_KE_ff_vs_fitted, hadana.pitype);
    FillHistVec2D(htrack_dKE_fitted_vs_Truth_2D, true_KE_ff, dKE_true_KE_ff_vs_fitted, hadana.pitype);
  }
  
  // == For loop for daugthers : pions and protons
  // == Reco
  vector<int> pion_index_vec;
  vector<int> p_index_vec;
  for (size_t i = 0; i < evt.reco_daughter_allTrack_ID->size(); i++){
    if((*evt.reco_daughter_allTrack_dQdX_SCE)[i].empty()) continue;
    //cout << "reco_daughter_allTrack_ID->at(" << i << ") : " << evt.reco_daughter_allTrack_ID->at(i) << endl;
    double michelscore = -1;
    unsigned int N_ranges = (*evt.reco_daughter_allTrack_resRange_SCE).at(i).size();
    unsigned int N_ranges_dEdX = (*evt.reco_daughter_allTrack_calibrated_dEdX_SCE).at(i).size();
    //cout << "[N_ranges:N_ranges_dEdX]\t[" << N_ranges << ":" << N_ranges_dEdX << "]" << endl;
    if(N_ranges < 20) continue; // == at leat 20 hits for a track
    

    // == Truth level distribution
    if(!evt.MC) continue;
    // = proton
    //cout << "evt.reco_daughter_PFP_true_byHits_PDG : " << (*evt.reco_daughter_PFP_true_byHits_PDG).at(i) << endl;
    if((*evt.reco_daughter_PFP_true_byHits_PDG).at(i) == 2212){
      hdaughter_proton_trackScore -> Fill((*evt.reco_daughter_PFP_trackScore).at(i));
      hdaughter_proton_emScore -> Fill((*evt.reco_daughter_PFP_emScore).at(i));
      double this_chi2 = (*evt.reco_daughter_allTrack_Chi2_proton).at(i) / (*evt.reco_daughter_allTrack_Chi2_ndof).at(i);
      //cout << "proton this_chi2 : " << this_chi2 << endl;
      hdaughter_proton_chi2_proton -> Fill(this_chi2);
    }
    else if(abs((*evt.reco_daughter_PFP_true_byHits_PDG).at(i)) == 211){
      hdaughter_pion_trackScore -> Fill((*evt.reco_daughter_PFP_trackScore).at(i));
      hdaughter_pion_emScore -> Fill((*evt.reco_daughter_PFP_emScore).at(i));
      double this_chi2 = (*evt.reco_daughter_allTrack_Chi2_proton).at(i) / (*evt.reco_daughter_allTrack_Chi2_ndof).at(i);
      hdaughter_pion_chi2_proton -> Fill(this_chi2);
      //cout << "pion this_chi2 : " << this_chi2 << endl;
    }
    else{
      hdaughter_other_trackScore -> Fill((*evt.reco_daughter_PFP_trackScore).at(i));
      hdaughter_other_emScore -> Fill((*evt.reco_daughter_PFP_emScore).at(i));
      double this_chi2 = (*evt.reco_daughter_allTrack_Chi2_proton).at(i) / (*evt.reco_daughter_allTrack_Chi2_ndof).at(i);
      hdaughter_other_chi2_proton -> Fill(this_chi2);
    }

    // == Truth matched
    for (size_t j = 0; j < evt.true_beam_daughter_ID->size(); j++){
      if (std::abs((*evt.true_beam_daughter_PDG)[j]) != 211 &&
	  (*evt.true_beam_daughter_PDG)[j] != 2212) continue;

      if ((*evt.reco_daughter_PFP_true_byHits_ID)[i] == (*evt.true_beam_daughter_ID)[j]){
	
	double true_mom = (*evt.true_beam_daughter_startP)[j]*1000;
	double true_theta = GetTheta((*evt.true_beam_daughter_startPx)[j],
				     (*evt.true_beam_daughter_startPy)[j],
				     (*evt.true_beam_daughter_startPz)[j]);
	double true_phi = GetPhi((*evt.true_beam_daughter_startPx)[j],
				 (*evt.true_beam_daughter_startPy)[j],
				 (*evt.true_beam_daughter_startPz)[j]);
	bool foundreco = false;
      }

    }
  }


  // == Truth
  for (size_t i = 0; i<evt.true_beam_daughter_ID->size(); ++i){
    if (std::abs((*evt.true_beam_daughter_PDG)[i]) != 211 &&
        (*evt.true_beam_daughter_PDG)[i] != 2212) continue;
    
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
    //if (jentry%100000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    if (jentry%1000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //std::cout<<evt.run<<" "<<evt.event<<" "<<evt.MC<<" "<<evt.reco_beam_true_byE_matched<<" "<<evt.true_beam_PDG<<" "<<(*evt.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    if (!hadana.isSelectedPart(evt)) continue;
    h_cutflow -> Fill(0.5);
    ProcessEvent(evt);
    FillHistograms_precut(evt);
    if (!hadana.PassPiCuts(evt)) continue;
    FillHistograms(evt);
  }
  SaveHistograms();
  //Make_dEdx_Range_Profile();
}

