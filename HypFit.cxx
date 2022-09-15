#include "HypFit.h"
#include "HadAna.h"
#include "anavar.h"
#include "TH1D.h"
#include "util.h"
#include <iostream>

using namespace std;

HypFit::HypFit(){
  hadana.InitPi();
  //hadana.InitP();
}

double HypFit::Gaussian_Reweight(double mu1, double sigma1, double mu2, double sigma2, double x_min, double x_max, double x){
  // == Estimate y of gaussian1 / gaussian2
  TF1 *gaus1 = new TF1("gaus1", "gaus", x_min, x_max);
  gaus1 -> SetParameter(0, 1. / (sigma1 * sqrt(2.0 * TMath::Pi())));
  gaus1 -> SetParameter(1, mu1);
  gaus1 -> SetParameter(2, sigma1);

  TF1 *gaus2 = new TF1("gaus2","gaus",x_min, x_max);
  gaus2 -> SetParameter(0, 1. / (sigma2 * sqrt(2.0 * TMath::Pi())));
  gaus2 -> SetParameter(1, mu2);
  gaus2 -> SetParameter(2, sigma2);

  double y1 = gaus1 -> Eval(x);
  double y2 = gaus2 -> Eval(x);
  double out = y1/y2;
  //cout << "[HypFit::Gaussian_Reweight] true beam P : " << x << ", y1 : " << y1 << ", y2 : " << y2 << ", y1/y2 : " << y1/y2 << endl;
  if(out > 10.) return 10.;
  return y1/y2;
}

vector<RecoDaughter> HypFit::GetAllRecoDaughters(const anavar & evt){

  // ==== Beam related variables
  TVector3 reco_beam_end(evt.reco_beam_endX, evt.reco_beam_endY, evt.reco_beam_endZ);
  TVector3 reco_unit_beam(evt.reco_beam_allTrack_trackDirX, evt.reco_beam_allTrack_trackDirY, evt.reco_beam_allTrack_trackDirZ);
  reco_unit_beam = (1. / reco_unit_beam.Mag() ) * reco_unit_beam;

  vector<RecoDaughter> out;

  //cout << "[HypFit::GetAllRecoDaughters] evt.reco_daughter_allTrack_ID->size() : " << evt.reco_daughter_allTrack_ID->size() << endl;
  for (size_t i = 0; i < evt.reco_daughter_allTrack_ID->size(); i++){
    if((*evt.reco_daughter_allTrack_dQdX_SCE)[i].empty()) continue;

    RecoDaughter this_RecoDaughter;
    //cout << "[HypFit::GetAllRecoDaughters] i : " << i << endl;
    this_RecoDaughter.SetIsEmpty(false);
    if(evt.MC){
      this_RecoDaughter.Set_PFP_true_byHits_PDG((*evt.reco_daughter_PFP_true_byHits_PDG).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_ID((*evt.reco_daughter_PFP_true_byHits_ID).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_origin((*evt.reco_daughter_PFP_true_byHits_origin).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_parID((*evt.reco_daughter_PFP_true_byHits_parID).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_parPDG((*evt.reco_daughter_PFP_true_byHits_parPDG).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_process((*evt.reco_daughter_PFP_true_byHits_process).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_sharedHits((*evt.reco_daughter_PFP_true_byHits_sharedHits).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_emHits((*evt.reco_daughter_PFP_true_byHits_emHits).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_len((*evt.reco_daughter_PFP_true_byHits_len).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_startX((*evt.reco_daughter_PFP_true_byHits_startX).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_startY((*evt.reco_daughter_PFP_true_byHits_startY).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_startZ((*evt.reco_daughter_PFP_true_byHits_startZ).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_endX((*evt.reco_daughter_PFP_true_byHits_endX).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_endY((*evt.reco_daughter_PFP_true_byHits_endY).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_endZ((*evt.reco_daughter_PFP_true_byHits_endZ).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_startPx((*evt.reco_daughter_PFP_true_byHits_startPx).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_startPy((*evt.reco_daughter_PFP_true_byHits_startPy).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_startPz((*evt.reco_daughter_PFP_true_byHits_startPz).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_startP((*evt.reco_daughter_PFP_true_byHits_startP).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_startE((*evt.reco_daughter_PFP_true_byHits_startE).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_endProcess((*evt.reco_daughter_PFP_true_byHits_endProcess).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_purity((*evt.reco_daughter_PFP_true_byHits_purity).at(i));
      this_RecoDaughter.Set_PFP_true_byHits_completeness((*evt.reco_daughter_PFP_true_byHits_completeness).at(i));
      this_RecoDaughter.Set_PFP_true_byE_PDG((*evt.reco_daughter_PFP_true_byE_PDG).at(i));
      this_RecoDaughter.Set_PFP_true_byE_len((*evt.reco_daughter_PFP_true_byE_len).at(i));
    }

    this_RecoDaughter.Set_PFP_ID((*evt.reco_daughter_PFP_ID).at(i));
    this_RecoDaughter.Set_PFP_nHits((*evt.reco_daughter_PFP_nHits).at(i));
    this_RecoDaughter.Set_PFP_nHits_collection((*evt.reco_daughter_PFP_nHits_collection).at(i));
    this_RecoDaughter.Set_PFP_trackScore((*evt.reco_daughter_PFP_trackScore).at(i));
    this_RecoDaughter.Set_PFP_emScore((*evt.reco_daughter_PFP_emScore).at(i));
    this_RecoDaughter.Set_PFP_michelScore((*evt.reco_daughter_PFP_michelScore).at(i));
    this_RecoDaughter.Set_PFP_trackScore_collection((*evt.reco_daughter_PFP_trackScore_collection).at(i));
    this_RecoDaughter.Set_PFP_emScore_collection((*evt.reco_daughter_PFP_emScore_collection).at(i));
    this_RecoDaughter.Set_PFP_michelScore_collection((*evt.reco_daughter_PFP_michelScore_collection).at(i));
    this_RecoDaughter.Set_allTrack_ID((*evt.reco_daughter_allTrack_ID).at(i));
    this_RecoDaughter.Set_allTrack_EField_SCE((*evt.reco_daughter_allTrack_EField_SCE).at(i));
    this_RecoDaughter.Set_allTrack_resRange_SCE((*evt.reco_daughter_allTrack_resRange_SCE).at(i));
    this_RecoDaughter.Set_allTrack_resRange_SCE_plane0((*evt.reco_daughter_allTrack_resRange_plane0).at(i)); // == FIXME, to SCE
    this_RecoDaughter.Set_allTrack_resRange_SCE_plane1((*evt.reco_daughter_allTrack_resRange_plane1).at(i)); // == FIXME, to SCE
    this_RecoDaughter.Set_allTrack_calibrated_dEdX_SCE((*evt.reco_daughter_allTrack_calibrated_dEdX_SCE).at(i));
    this_RecoDaughter.Set_allTrack_calibrated_dEdX_SCE_plane0((*evt.reco_daughter_allTrack_calibrated_dEdX_SCE_plane0).at(i));
    this_RecoDaughter.Set_allTrack_calibrated_dEdX_SCE_plane1((*evt.reco_daughter_allTrack_calibrated_dEdX_SCE_plane1).at(i));
    this_RecoDaughter.Set_allTrack_Chi2_proton((*evt.reco_daughter_allTrack_Chi2_proton).at(i));
    this_RecoDaughter.Set_allTrack_Chi2_pion((*evt.reco_daughter_allTrack_Chi2_pion).at(i));
    this_RecoDaughter.Set_allTrack_Chi2_muon((*evt.reco_daughter_allTrack_Chi2_muon).at(i));
    this_RecoDaughter.Set_allTrack_Chi2_ndof((*evt.reco_daughter_allTrack_Chi2_ndof).at(i));
    this_RecoDaughter.Set_allTrack_Chi2_ndof_pion((*evt.reco_daughter_allTrack_Chi2_ndof_pion).at(i));
    this_RecoDaughter.Set_allTrack_Chi2_ndof_muon((*evt.reco_daughter_allTrack_Chi2_ndof_muon).at(i));
    this_RecoDaughter.Set_allTrack_Theta((*evt.reco_daughter_allTrack_Theta).at(i));
    this_RecoDaughter.Set_allTrack_Phi((*evt.reco_daughter_allTrack_Phi).at(i));
    // == FIXME
    //this_RecoDaughter.Set_allTrack_startDirX((*evt.reco_daughter_allTrack_startDirX).at(i));
    //this_RecoDaughter.Set_allTrack_startDirY((*evt.reco_daughter_allTrack_startDirY).at(i));
    //this_RecoDaughter.Set_allTrack_startDirZ((*evt.reco_daughter_allTrack_startDirZ).at(i));
    this_RecoDaughter.Set_allTrack_alt_len((*evt.reco_daughter_allTrack_alt_len).at(i));
    this_RecoDaughter.Set_allTrack_startX((*evt.reco_daughter_allTrack_startX).at(i));
    this_RecoDaughter.Set_allTrack_startY((*evt.reco_daughter_allTrack_startY).at(i));
    this_RecoDaughter.Set_allTrack_startZ((*evt.reco_daughter_allTrack_startZ).at(i));
    this_RecoDaughter.Set_allTrack_endX((*evt.reco_daughter_allTrack_endX).at(i));
    this_RecoDaughter.Set_allTrack_endY((*evt.reco_daughter_allTrack_endY).at(i));
    this_RecoDaughter.Set_allTrack_endZ((*evt.reco_daughter_allTrack_endZ).at(i));
    this_RecoDaughter.Set_allTrack_vertex_michel_score((*evt.reco_daughter_allTrack_vertex_michel_score).at(i));
    this_RecoDaughter.Set_allTrack_vertex_nHits((*evt.reco_daughter_allTrack_vertex_nHits).at(i));
    this_RecoDaughter.Set_pandora_type((*evt.reco_daughter_pandora_type).at(i));

    TVector3 unit_daughter((*evt.reco_daughter_allTrack_endX).at(i) - (*evt.reco_daughter_allTrack_startX).at(i),
			   (*evt.reco_daughter_allTrack_endY).at(i) - (*evt.reco_daughter_allTrack_startY).at(i),
			   (*evt.reco_daughter_allTrack_endZ).at(i) - (*evt.reco_daughter_allTrack_startZ).at(i) );
    unit_daughter = (1./ unit_daughter.Mag() ) * unit_daughter;
    double cos_theta = cos(unit_daughter.Angle(reco_unit_beam));
    this_RecoDaughter.Set_Beam_Cos(cos_theta);

    TVector3 reco_daughter_start((*evt.reco_daughter_allTrack_startX).at(i), (*evt.reco_daughter_allTrack_startY).at(i), (*evt.reco_daughter_allTrack_startZ).at(i) );
    double dist_beam_end = (reco_daughter_start - reco_beam_end).Mag();
    this_RecoDaughter.Set_Beam_Dist(dist_beam_end);

    out.push_back(this_RecoDaughter);
  }

  return out;
}

vector<RecoDaughter> HypFit::GetPions(const vector<RecoDaughter> in){

  vector<RecoDaughter> out;

  double cut_cos_beam = 0.95;
  double cut_beam_dist = 10.;
  double cut_trackScore = 0.5;
  double cut_emScore = 0.5;
  double cut_chi2_proton = 60.;
  double cut_startZ = 220.;
  int cut_Nhit = 20;
  for(unsigned int i = 0; i < in.size(); i++){
    RecoDaughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() / this_in.allTrack_Chi2_ndof();
    if(this_in.PFP_trackScore() > cut_trackScore && this_in.PFP_emScore() < cut_emScore && this_chi2 > cut_chi2_proton
       && this_in.PFP_nHits() > cut_Nhit && this_in.Beam_Cos() < cut_cos_beam && this_in.Beam_Dist() < cut_beam_dist && this_in.allTrack_startZ() < cut_startZ){
      out.push_back(this_in);
    }
  }
  
  return out;
}

vector<RecoDaughter> HypFit::GetProtons(const vector<RecoDaughter> in){

  vector<RecoDaughter> out;

  double cut_cos_beam =0.95;
  double cut_beam_dist = 10.;
  double cut_trackScore = 0.5;
  double cut_emScore = 0.5;
  double cut_chi2_proton = 50.;
  int cut_Nhit = 20;
  for(unsigned int i = 0; i < in.size(); i++){
    RecoDaughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() /this_in.allTrack_Chi2_ndof();
    if(this_in.PFP_trackScore() > cut_trackScore && this_in.PFP_emScore() < cut_emScore && this_chi2 < cut_chi2_proton && this_in.PFP_nHits() > cut_Nhit && this_in.Beam_Cos() < cut_cos_beam && this_in.Beam_Dist() < cut_beam_dist){
      out.push_back(this_in);
    }
  }

  return out;
}

vector<RecoDaughter> HypFit::GetTruePions(const vector<RecoDaughter> in){

  vector<RecoDaughter> out;

  for(unsigned int i = 0; i < in.size(); i++){
    RecoDaughter this_in = in.at(i);
    int this_true_PID = this_in.PFP_true_byHits_PDG();
    if(abs(this_true_PID) == 211){
      out.push_back(this_in);
    }
  }

  return out;
}

vector<RecoDaughter> HypFit::GetTrueProtons(const vector<RecoDaughter> in){

  vector<RecoDaughter> out;

  for(unsigned int i = 0; i < in.size(); i++){
    RecoDaughter this_in = in.at(i);
    int this_true_PID = this_in.PFP_true_byHits_PDG();
    if(this_true_PID == 2212){
      out.push_back(this_in);
    }
  }

  return out;
}

void HypFit::BookHistograms(){
  //cout << "[HypFit::BookHistograms] Start" << endl;
  Hist.outfile = TFile::Open(fOutputFileName.c_str(), "recreate");
}

void HypFit::ProcessEvent(const anavar & evt){
  hadana.ProcessEvent(evt);
}

void HypFit::FillHistBeam(const anavar & evt, double weight, int PDGID, TString suffix){

  TString particle_str = "";
  if(PDGID == 2212) particle_str = "proton";
  else if(PDGID == 211) particle_str = "pion";
  else particle_str = "other";

  double beamP = evt.beam_inst_P*1000. * beamP_scale;
  double beamKE = hadana.map_BB[PDGID]->MomentumtoKE(beamP);

  // ==== Define strings
  int N_hits_X = (*evt.reco_beam_calibrated_dEdX_SCE).size();
  int N_hits_X_range = (*evt.reco_beam_resRange_SCE).size();
  //cout << "[HypFit::FillHistBeam] N_hits_X : " << N_hits_X << ", N_hits_X_range : " << N_hits_X_range << endl;
  TString beamKE_range_str = "";
  if(beamKE < 100.){
    int beamKE_int50 = beamKE / 50;
    beamKE_range_str = Form("KE%dto%d", beamKE_int50 * 50, (beamKE_int50 + 1) * 50);
  }
  else{
    int beamKE_int100 = beamKE / 100;
    beamKE_range_str = Form("KE%dto%d", beamKE_int100 * 100, (beamKE_int100 + 1) * 100 );
  }
  TString N_hits_X_str = "";
  int N_hits_int30 = N_hits_X / 30;
  if(N_hits_int30 == 0) N_hits_X_str = Form("NHits15to%d", (N_hits_int30 + 1) * 30);
  else N_hits_X_str = Form("NHits%dto%d", N_hits_int30 * 30, (N_hits_int30 + 1) * 30);

  // == Reordering dE/dx and range vectors
  vector<double> this_dEdx_vec;
  vector<double> this_range_vec;
  for(unsigned int i_hit = N_hits_X - 1; i_hit > 0; i_hit--){
    //cout << "[HypFit::FillHistBeam] i_hit : " << i_hit << endl;
    this_dEdx_vec.push_back((*evt.reco_beam_calibrated_dEdX_SCE).at(i_hit));
    this_range_vec.push_back((*evt.reco_beam_resRange_SCE).at(i_hit));
  }
  //for(unsigned int i_hit = 0; i_hit < this_dEdx_vec.size(); i_hit++){
  //cout << "this_dEdx_vec.at(" << i_hit << ") : " << this_dEdx_vec.at(i_hit) << ", this_range_vec.at(" << i_hit << ") : " << this_range_vec.at(i_hit) << endl;
  //}

  double KE_range = hadana.map_BB[PDGID]->KEFromRangeSpline(evt.reco_beam_alt_len);
  double KE_res_range = (KE_range - beamKE) / beamKE;

  double length_gaus = hadana.Fit_dEdx_Residual_Length(evt, this_dEdx_vec, this_range_vec, PDGID, false, false);
  double length_likelihood = hadana.Fit_Residual_Length_Likelihood(evt, this_dEdx_vec, this_range_vec, PDGID, false);
  double KE_fit_gaussian = -9999.;
  double KE_fit_likelihood = -9999.;
  if(length_gaus > 0.) KE_fit_gaussian = hadana.map_BB[PDGID] -> KEFromRangeSpline(length_gaus);
  if(length_likelihood > 0.) KE_fit_likelihood = hadana.map_BB[PDGID] -> KEFromRangeSpline(length_likelihood);
  double KE_res_gaussian = -9999.;
  double KE_res_likelihood = -9999.;
  if(KE_fit_gaussian > 0.) KE_res_gaussian = (KE_fit_gaussian - beamKE) / beamKE;
  if(KE_fit_likelihood > 0.) KE_res_likelihood = (KE_fit_likelihood - beamKE) / beamKE;
  
  if(KE_fit_gaussian > 0.) cout << "[HypFit::FillHistBeam] beamKE : " << beamKE << ", length_gaus : " << length_gaus << ", KE_fit_gaussian : " << KE_fit_gaussian << ", KE_res_gaussian : " << KE_res_gaussian << endl;
  if(KE_fit_likelihood > 0.) cout << "[HypFit::FillHistBeam] beamKE : " << beamKE << ", length_likelihood : " << length_likelihood << ", KE_fit_likelihood : " << KE_fit_likelihood << ", KE_res_likelihood : " << KE_res_likelihood << endl;

  Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_1D_" + particle_str + "_range_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_range, weight, 400., -2., 2.);
  Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_2D_" + particle_str + "_range_" + N_hits_X_str + "_" + suffix, beamKE, KE_res_range, weight, 2000., 0., 2000., 400., -2., 2.);
  if(KE_fit_gaussian > 0.){
    Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_1D_" + particle_str + "_gaussian_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_gaussian, weight, 400., -2., 2.);
    Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_2D_" + particle_str + "_gaussian_" + N_hits_X_str + "_" + suffix, beamKE, KE_res_gaussian, weight, 2000., 0., 2000., 400., -2., 2.);
  }
  if(KE_fit_likelihood > 0.){
    Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_1D_" + particle_str + "_likelihood_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_likelihood, weight, 400., -2., 2.);
    Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_2D_" + particle_str + "_likelihood_" + N_hits_X_str + "_" + suffix, beamKE, KE_res_likelihood, weight, 2000., 0., 2000., 400., -2., 2.);
  }

  if(hadana.daughter_michel_score > 0.5){
    Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_1D_" + particle_str + "_stopped_range_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_range, weight, 400., -2., 2.);
    Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_2D_" + particle_str + "_stopped_range_" + N_hits_X_str + "_" + suffix, beamKE, KE_res_range, weight, 2000., 0., 2000., 400., -2., 2.);
    if(KE_fit_gaussian > 0.){
      Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_1D_" + particle_str + "_stopped_gaussian_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_gaussian, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_2D_" + particle_str + "_stopped_gaussian_" + N_hits_X_str + "_" + suffix, beamKE, KE_res_gaussian, weight, 2000., 0., 2000., 400., -2., 2.);
    }
    if(KE_fit_likelihood > 0.){
      Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_1D_" + particle_str + "_stopped_likelihood_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_likelihood, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_2D_" + particle_str + "_stopped_likelihood_" + N_hits_X_str + "_" + suffix, beamKE, KE_res_likelihood, weight, 2000., 0., 2000., 400., -2., 2.);
    }
  }

  if(evt.MC){
    double true_KE_ff = hadana.true_ffKE;
    TString true_beam_PDG_str = Form("%d", abs(evt.true_beam_PDG));
    double KE_true_res_range = (KE_range - true_KE_ff) / true_KE_ff;
    double KE_true_res_gaussian = -9999.;
    double KE_true_res_likelihood = -9999.;
    if(KE_fit_gaussian > 0.) KE_true_res_gaussian = (KE_fit_likelihood - true_KE_ff) / true_KE_ff;
    if(KE_fit_likelihood > 0.) KE_true_res_likelihood = (KE_fit_likelihood - true_KE_ff) / true_KE_ff;
    Hist.JSFillHist(suffix, "Beam_true_KE_Res_1D_" + particle_str + "_range_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_range, weight, 400., -2., 2.);
    Hist.JSFillHist(suffix, "Beam_true_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_range_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_range, weight, 400., -2., 2.);
    Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_range_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_range, weight, 400., -2., 2.);

    if(KE_fit_gaussian > 0.){
      Hist.JSFillHist(suffix, "Beam_true_KE_Res_1D_" + particle_str + "_gaussian_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_gaussian, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Beam_true_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_gaussian_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_gaussian, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_gaussian_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_gaussian, weight, 400., -2., 2.);
    }
    if(KE_fit_likelihood > 0.){
      Hist.JSFillHist(suffix, "Beam_true_KE_Res_1D_" + particle_str + "_likelihood_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_likelihood, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Beam_true_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_likelihood_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_likelihood, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Beam_Reco_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_likelihood_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_likelihood, weight, 400., -2., 2.);
    }
  }
}

void HypFit::FillHistBeam_Micheless(const anavar & evt, double weight, int PDGID, TString suffix){

  if(hadana.daughter_michel_score < 0.85) return;

  TString particle_str = "";
  if(PDGID == 2212) particle_str = "proton";
  else if(PDGID == 211) particle_str = "pion";
  else if(PDGID == 13) particle_str = "muon";
  else particle_str = "other";

  double beamP = evt.beam_inst_P*1000. * beamP_scale;
  double beamKE = hadana.map_BB[PDGID]->MomentumtoKE(beamP);

  // ==== Define strings 
  int N_hits_X = (*evt.reco_beam_calibrated_dEdX_SCE).size();
  int N_hits_X_range = (*evt.reco_beam_resRange_SCE).size();
  TString beamKE_range_str = "";
  if(beamKE < 100.){
    int beamKE_int50 = beamKE / 50;
    beamKE_range_str = Form("KE%dto%d", beamKE_int50 * 50, (beamKE_int50 + 1) * 50);
  }
  else{
    int beamKE_int100 = beamKE / 100;
    beamKE_range_str = Form("KE%dto%d", beamKE_int100 * 100, (beamKE_int100 + 1) * 100 );
  }
  TString N_hits_X_str = "";
  int N_hits_int30 = N_hits_X / 30;
  if(N_hits_int30 == 0) N_hits_X_str = Form("NHits15to%d", (N_hits_int30 + 1) * 30);
  else N_hits_X_str = Form("NHits%dto%d", N_hits_int30 * 30, (N_hits_int30 + 1) * 30);

  // == Reordering dE/dx and range vectors
  vector<double> this_dEdx_vec; 
  vector<double> this_range_vec;
  for(unsigned int i_hit = N_hits_X - 1; i_hit > 0; i_hit--){
    this_dEdx_vec.push_back((*evt.reco_beam_calibrated_dEdX_SCE).at(i_hit));
    this_range_vec.push_back((*evt.reco_beam_resRange_SCE).at(i_hit));
  }
  double KE_range = hadana.map_BB[PDGID]->KEFromRangeSpline(evt.reco_beam_alt_len);
  double KE_res_range = (KE_range - beamKE) / beamKE;

  double length_gaus = hadana.Fit_dEdx_Residual_Length(evt, this_dEdx_vec, this_range_vec, PDGID, false, false);
  double length_likelihood = hadana.Fit_Residual_Length_Likelihood(evt, this_dEdx_vec, this_range_vec, PDGID, false);
  double KE_fit_gaussian = -9999.;
  double KE_fit_likelihood = -9999.;
  if(length_gaus > 0.) KE_fit_gaussian = hadana.map_BB[PDGID] -> KEFromRangeSpline(length_gaus);
  if(length_likelihood > 0.) KE_fit_likelihood = hadana.map_BB[PDGID] -> KEFromRangeSpline(length_likelihood);
  double KE_res_gaussian = -9999.;
  double KE_res_likelihood = -9999.;
  if(KE_fit_gaussian > 0.) KE_res_gaussian = (KE_fit_gaussian - beamKE) / beamKE;
  if(KE_fit_likelihood > 0.) KE_res_likelihood = (KE_fit_likelihood - beamKE) / beamKE;

  Hist.JSFillHist(suffix, "Micheless_Beam_Reco_KE_Res_1D_" + particle_str + "_range_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_range, weight, 400., -2., 2.);
  if(KE_fit_gaussian > 0.){
    Hist.JSFillHist(suffix, "Micheless_Beam_Reco_KE_Res_1D_" + particle_str + "_gaussian_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_gaussian, weight, 400., -2., 2.);
  }
  if(KE_fit_likelihood > 0.){
    Hist.JSFillHist(suffix, "Micheless_Beam_Reco_KE_Res_1D_" + particle_str + "_likelihood_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_likelihood, weight, 400., -2., 2.);
  }

  if(KE_fit_gaussian > 0.) cout << "[HypFit::FillHistBeam_Micheless] beamKE : " << beamKE << ", length_gaus : " << length_gaus << ", KE_fit_gaussian : " << KE_fit_gaussian << ", KE_res_gaussian : " << KE_res_gaussian << endl;
  if(KE_fit_likelihood > 0.) cout << "[HypFit::FillHistBeam_Micheless] beamKE : " << beamKE << ", length_likelihood : " << length_likelihood << ", KE_fit_likelihood : " << KE_fit_likelihood << ", KE_res_likelihood : " << KE_res_likelihood << endl;

  if(evt.MC){
    double true_KE_ff = hadana.true_ffKE;
    TString true_beam_PDG_str = Form("%d", abs(evt.true_beam_PDG));
    double KE_true_res_range = (KE_range - true_KE_ff) / true_KE_ff;
    double KE_true_res_gaussian = -9999.;
    double KE_true_res_likelihood = -9999.;
    if(KE_fit_gaussian > 0.) KE_true_res_gaussian = (KE_fit_likelihood - true_KE_ff) / true_KE_ff;
    if(KE_fit_likelihood > 0.) KE_true_res_likelihood = (KE_fit_likelihood - true_KE_ff) / true_KE_ff;
    Hist.JSFillHist(suffix, "Micheless_Beam_true_KE_Res_1D_" + particle_str + "_range_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_range, weight, 400., -2., 2.);
    Hist.JSFillHist(suffix, "Micheless_Beam_true_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_range_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_range, weight, 400., -2., 2.);
    Hist.JSFillHist(suffix, "Micheless_Beam_Reco_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_range_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_range, weight, 400., -2., 2.);

    if(KE_fit_gaussian > 0.){
      Hist.JSFillHist(suffix, "Micheless_Beam_true_KE_Res_1D_" + particle_str + "_gaussian_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_gaussian, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Micheless_Beam_true_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_gaussian_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_gaussian, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Micheless_Beam_Reco_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_gaussian_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_gaussian, weight, 400., -2., 2.);
    }
    if(KE_fit_likelihood > 0.){
      Hist.JSFillHist(suffix, "Micheless_Beam_true_KE_Res_1D_" + particle_str + "_likelihood_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_likelihood, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Micheless_Beam_true_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_likelihood_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_true_res_likelihood, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Micheless_Beam_Reco_KE_Res_1D_" + particle_str + "_PID" + true_beam_PDG_str + "_likelihood_" + beamKE_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_likelihood, weight, 400., -2., 2.);
    }
  }
}  

void HypFit::FillHistDaughterTrue(const vector<RecoDaughter> daughters, const anavar & evt, int PDGID, double weight, TString suffix){

  // ==== Basic variables for daughter selection
  for(unsigned int i = 0; i < daughters.size(); i++){
    RecoDaughter this_daughter = daughters.at(i);
    TString particle_str = "";
    int this_PdgID = abs(this_daughter.PFP_true_byHits_PDG());
    if(this_PdgID == 2212) particle_str = "proton";
    else if(abs(this_PdgID) == 211) particle_str = "pion";
    else if(abs(this_PdgID) == 13) particle_str = "muon";
    else particle_str = "other";

    // == Daughter hits : this_daughter.allTrack_calibrated_dEdX_SCE() and this_daughter.allTrack_resRange_SCE()
    // == Truth level info
    // ==== Momentum : this_daughter.PFP_true_byHits_startP() * 1000.
    // ==== KE : hadana.map_BB[this_PdgID]->MomentumtoKE(Momentum)
    //if(this_PdgID != 2212 && this_PdgID != 211) continue; // == Run only for protons and pions
    if(this_PdgID != PDGID) continue;
    double P_true = this_daughter.PFP_true_byHits_startP() * 1000.;
    double KE_true = hadana.map_BB[this_PdgID]->MomentumtoKE(P_true);
    int N_hits_X = this_daughter.allTrack_calibrated_dEdX_SCE().size();
    double KE_range = hadana.map_BB[this_PdgID]->KEFromRangeSpline(this_daughter.allTrack_alt_len());
    double KE_res_range = (KE_range - KE_true) / KE_true;

    double Expected_range = hadana.map_BB[this_PdgID]->RangeFromKESpline(KE_true);
    double range_ratio = this_daughter.allTrack_alt_len() / Expected_range;
    
    // == Define strings for histogram names
    TString KE_true_range_str = "";
    if(KE_true < 100.){
      int KE_true_int50 = KE_true / 50;
      KE_true_range_str = Form("KE%dto%d", KE_true_int50 * 50, (KE_true_int50 + 1) * 50);
    }
    else{
      int KE_true_int100 = KE_true / 100;
      KE_true_range_str = Form("KE%dto%d", KE_true_int100 * 100, (KE_true_int100 + 1) * 100 );
    }
    TString N_hits_X_str = "";
    int N_hits_int30 = N_hits_X / 30;
    if(N_hits_int30 == 0) N_hits_X_str = Form("NHits15to%d", (N_hits_int30 + 1) * 30);
    else N_hits_X_str = Form("NHits%dto%d", N_hits_int30 * 30, (N_hits_int30 + 1) * 30);
    TString range_ratio_str = "";
    if(range_ratio < 0.1) range_ratio_str = "LessThan0p1";
    else if(range_ratio > 0.1 && range_ratio < 0.3) range_ratio_str = "0p1to0p3";
    else if(range_ratio > 0.3 && range_ratio < 0.5) range_ratio_str = "0p3to0p5";
    else if(range_ratio > 0.5 && range_ratio < 0.7) range_ratio_str = "0p5to0p7";
    else if(range_ratio > 0.7 && range_ratio < 0.9) range_ratio_str = "0p7to0p9";
    else range_ratio_str = "BiggerThan0p9";

    if(N_hits_X <= 15) continue; // == Remove tracks with too small number of hits 
    if(evt.true_beam_ID == this_daughter.PFP_true_byHits_ID()) continue; // == Remove broken tracks
    //cout << "[HypFit::FillHistDaughterTrue] P_true : " << P_true << ", KE_true : " << KE_true << ", KE_true_range_str : " << KE_true_range_str << ", N_hits_X : " << N_hits_X << ", N_hists_X_str : " << N_hits_X_str << endl;

    double daughter_length_gaus = hadana.Fit_dEdx_Residual_Length(evt, this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), this_PdgID, false, false);
    double daughter_length_likelihood = hadana.Fit_Residual_Length_Likelihood(evt, this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), this_PdgID, false);
    double KE_fit_gaussian = -9999.;
    double KE_fit_likelihood = -9999.;
    if(daughter_length_gaus > 0.) KE_fit_gaussian = hadana.map_BB[this_PdgID] -> KEFromRangeSpline(daughter_length_gaus);
    if(daughter_length_likelihood > 0.) KE_fit_likelihood = hadana.map_BB[this_PdgID] -> KEFromRangeSpline(daughter_length_likelihood);
    double KE_res_gaussian = -9999.;
    double KE_res_likelihood = -9999.;
    if(KE_fit_gaussian > 0.) KE_res_gaussian = (KE_fit_gaussian - KE_true) / KE_true;
    if(KE_fit_likelihood > 0.) KE_res_likelihood = (KE_fit_likelihood - KE_true) / KE_true;
    
    // == Fill histograms for efficiency
    Hist.JSFillHist(suffix, "KE_eff_" + particle_str + "_range_" + range_ratio_str + "_" + suffix, KE_true, weight, 100., 0., 1000.);
    if(KE_fit_gaussian > 0.){
      Hist.JSFillHist(suffix, "KE_eff_" + particle_str + "_gaussian_" + range_ratio_str + "_" + suffix, KE_true, weight, 100., 0., 1000.);
      //cout << "[HypFit::FillHistDaughterTrue] Fitting Efficienty fill for Gaussian with range ratio " << range_ratio_str << endl;
    }
    if(KE_fit_likelihood > 0.){
      Hist.JSFillHist(suffix, "KE_eff_" + particle_str + "_likelihood_" + range_ratio_str + "_" + suffix, KE_true, weight, 100., 0., 1000.);
      //cout << "[HypFit::FillHistDaughterTrue] Fitting Efficienty fill for Likelihood with range ratio " << range_ratio_str << endl;
    }

    // == Fill histograms for all
    Hist.JSFillHist(suffix, "KE_Res_1D_" + particle_str + "_range_" + KE_true_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_range, weight, 400., -2., 2.);
    Hist.JSFillHist(suffix, "KE_Res_2D_" + particle_str + "_range_" + N_hits_X_str + "_" + suffix, KE_true, KE_res_range, weight, 1000., 0., 1000., 400., -2., 2.);
    if(KE_fit_gaussian > 0.){
      Hist.JSFillHist(suffix, "KE_Res_1D_" + particle_str + "_gaussian_" + KE_true_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_gaussian, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "KE_Res_2D_" + particle_str + "_gaussian_" + N_hits_X_str + "_" + suffix, KE_true, KE_res_gaussian, weight, 1000., 0., 1000., 400., -2., 2.);
    }
    if(KE_fit_likelihood > 0.){
      Hist.JSFillHist(suffix, "KE_Res_1D_" + particle_str + "_likelihood_" + KE_true_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_likelihood, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "KE_Res_2D_" + particle_str + "_likelihood_" + N_hits_X_str + "_" + suffix, KE_true, KE_res_likelihood, weight, 1000., 0., 1000., 400., -2., 2.);
    }

    // == Fill histograms for stopped and interacted particles
    double this_michelscore = this_daughter.allTrack_vertex_michel_score() / this_daughter.allTrack_vertex_nHits(); 
    if(this_michelscore > 0.5){
      Hist.JSFillHist(suffix, "KE_Res_1D_" + particle_str + "_stopped_range_" + KE_true_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_range, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "KE_Res_2D_" + particle_str + "_stopped_range_" + N_hits_X_str + "_" + suffix, KE_true, KE_res_range, weight, 1000., 0., 1000., 400., -2., 2.);
      if(KE_fit_gaussian > 0.){
	Hist.JSFillHist(suffix, "KE_Res_1D_" + particle_str + "_stopped_gaussian_" + KE_true_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_gaussian, weight, 400., -2., 2.);
	Hist.JSFillHist(suffix, "KE_Res_2D_" + particle_str + "_stopped_gaussian_" + N_hits_X_str + "_" + suffix, KE_true, KE_res_gaussian, weight, 1000., 0., 1000., 400., -2., 2.);
      }
      if(KE_fit_likelihood > 0.){
	Hist.JSFillHist(suffix, "KE_Res_1D_" + particle_str + "_stopped_likelihood_" + KE_true_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_likelihood, weight, 400., -2., 2.);
	Hist.JSFillHist(suffix, "KE_Res_2D_" + particle_str + "_stopped_likelihood_" + N_hits_X_str + "_" + suffix, KE_true, KE_res_likelihood, weight, 1000., 0., 1000., 400., -2., 2.);
      }
    }
    else{
      Hist.JSFillHist(suffix, "KE_Res_1D_" + particle_str + "_interacted_range_" + KE_true_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_range, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "KE_Res_2D_" + particle_str + "_interacted_range_" + N_hits_X_str + "_" + suffix, KE_true, KE_res_range, weight, 1000., 0., 1000., 400., -2., 2.);
      if(KE_fit_gaussian > 0.){
	Hist.JSFillHist(suffix, "KE_Res_1D_" + particle_str + "_interacted_gaussian_" + KE_true_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_gaussian, weight, 400., -2., 2.);
	Hist.JSFillHist(suffix, "KE_Res_2D_" + particle_str + "_interacted_gaussian_" + N_hits_X_str + "_" + suffix, KE_true, KE_res_gaussian, weight, 1000., 0., 1000., 400., -2., 2.);

      }
      if(KE_fit_likelihood > 0.){
	Hist.JSFillHist(suffix, "KE_Res_1D_" + particle_str + "_interacted_likelihood_" + KE_true_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_likelihood, weight, 400., -2., 2.);
	Hist.JSFillHist(suffix, "KE_Res_2D_" + particle_str + "_interacted_likelihood_" + N_hits_X_str + "_" + suffix, KE_true, KE_res_likelihood, weight, 1000., 0., 1000., 400., -2., 2.);
      }
    }
  }
}

void HypFit::FillHistDaughterReco(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix){

  TString pitype_str = Form("%d", hadana.pitype);
  vector<RecoDaughter> pions = GetPions(daughters);

  // == Loop for pions
  for(unsigned int i = 0; i < pions.size(); i++){
    RecoDaughter this_pion = pions.at(i);
    double this_michelscore = this_pion.allTrack_vertex_michel_score() / this_pion.allTrack_vertex_nHits();
    if(this_michelscore < 0.5) continue;
    double KE_daughter = hadana.map_BB[211]->KEFromRangeSpline(this_pion.allTrack_alt_len());

    int N_hits_X = this_pion.allTrack_calibrated_dEdX_SCE().size();

    vector<double> this_dEdx_half;
    vector<double> this_range_half;
    double last_range = this_pion.allTrack_resRange_SCE().at(N_hits_X / 2);
    for(unsigned int i_hit = N_hits_X / 2; i_hit < N_hits_X; i_hit++){
      this_dEdx_half.push_back(this_pion.allTrack_calibrated_dEdX_SCE().at(i_hit));
      this_range_half.push_back(this_pion.allTrack_resRange_SCE().at(i_hit));
    }
    for(unsigned int i_hit = 0; i_hit < this_range_half.size(); i_hit++){
      this_range_half[i_hit] = this_range_half[i_hit] - last_range;
    }
    
    double daughter_length_gaus = hadana.Fit_dEdx_Residual_Length(evt, this_dEdx_half, this_range_half, 211, false, false);
    double daughter_length_likelihood = hadana.Fit_Residual_Length_Likelihood(evt, this_dEdx_half, this_range_half, 211, false);
    double KE_fit_gaussian = -9999.;
    double KE_fit_likelihood = -9999.;
    if(daughter_length_gaus > 0.) KE_fit_gaussian = hadana.map_BB[211] -> KEFromRangeSpline(daughter_length_gaus);
    if(daughter_length_likelihood > 0.) KE_fit_likelihood = hadana.map_BB[211] -> KEFromRangeSpline(daughter_length_likelihood);
    double KE_res_gaussian = -9999.;
    double KE_res_likelihood = -9999.;
    if(KE_fit_gaussian > 0.) KE_res_gaussian = (KE_fit_gaussian - KE_daughter) / KE_daughter;
    if(KE_fit_likelihood > 0.) KE_res_likelihood = (KE_fit_likelihood - KE_daughter) / KE_daughter;

    TString KE_daughter_range_str = "";
    if(KE_daughter < 100.){
      int KE_daughter_int50 = KE_daughter / 50;
      KE_daughter_range_str = Form("KE%dto%d", KE_daughter_int50 * 50, (KE_daughter_int50 + 1) * 50);
    }
    else{
      int KE_daughter_int100 = KE_daughter / 100;
      KE_daughter_range_str = Form("KE%dto%d", KE_daughter_int100 * 100, (KE_daughter_int100 + 1) * 100 );
    }
    TString N_hits_X_str = "";
    int N_hits_int30 = N_hits_X / 30;
    if(N_hits_int30 == 0) N_hits_X_str = Form("NHits15to%d", (N_hits_int30 + 1) * 30);
    else N_hits_X_str = Form("NHits%dto%d", N_hits_int30 * 30, (N_hits_int30 + 1) * 30);

    if(KE_fit_gaussian > 0.){
      //cout << KE_daughter_range_str <<", KE_res_gaussian : " << KE_res_gaussian << endl;
      Hist.JSFillHist(suffix, "Reco_KE_Res_1D_pion_gaussian_" + KE_daughter_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_gaussian, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Reco_KE_Res_2D_pion_gaussian_" + N_hits_X_str + "_" + suffix, KE_daughter, KE_res_gaussian, weight, 1000., 0., 1000., 400., -2., 2.);
    }
    if(KE_fit_likelihood > 0.){
      //cout << KE_daughter_range_str << ", KE_res_likelihood : " << KE_res_likelihood << endl;
      Hist.JSFillHist(suffix, "Reco_KE_Res_1D_pion_likelihood_" + KE_daughter_range_str + "_" + N_hits_X_str + "_" + suffix, KE_res_likelihood, weight, 400., -2., 2.);
      Hist.JSFillHist(suffix, "Reco_KE_Res_2D_pion_likelihood_" + N_hits_X_str + "_" + suffix, KE_daughter, KE_res_likelihood, weight, 1000., 0., 1000., 400., -2., 2.);
    }
  }
}




void HypFit::SaveHistograms(){
  
  //outputFile->cd();
  //outputFile->Write();
  Hist.WriteHist();
}

void HypFit::Run(anavar & evt, Long64_t nentries=-1){

  BookHistograms();
  TString beam_P_str = "2.0";
  double beam_P_inst_cut_upper = 9999.;
  double beam_P_inst_cut_lower = 0.;
  if(beam_P_str == "2.0"){
    beam_P_inst_cut_upper = 2400.;
    beam_P_inst_cut_lower = 1600.;
  }
  else if(beam_P_str == "1.0"){
    beam_P_inst_cut_upper = 1200.;
    beam_P_inst_cut_lower = 800.;
  }
  else if(beam_P_str == "0.5"){
    beam_P_inst_cut_upper = 600.;
    beam_P_inst_cut_lower = 400.;
  }
  else{
    cout << "[HypFit::Run] Not Valid Beam Momentum" << endl;
    return;
  }
  cout << "[HypFit::Run] Beam Momentum is " << beam_P_str << " GeV" << endl;


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
    Hist.FillHist("Cutflow", 0.5, 1., 10., 0., 10.);
    ProcessEvent(evt);

    double weight = 1.;
    double weight_piOnly = 1.;
    if(evt.MC){
      // == Beam P scale
      if(beam_P_str == "2.0") beamP_scale = 2.0;
      if(beam_P_str == "1.0") beamP_scale = 1.0;
      if(beam_P_str == "0.5") beamP_scale = 0.5;

      double beamP_true = evt.true_beam_startP * 1000.;

      double data_mu = 1.;
      double data_sigma = 1.;
      double true_mu = 1.;
      double true_sigma = 1.;
      // == 0.5 GeV
      if(beam_P_str == "0.5"){
	data_mu = 508.9;
	data_sigma = 36.66;
	true_mu = 504.7;
	true_sigma = 28.3;
      }
      // == 1.0 GeV
      if(beam_P_str == "1.0"){
	data_mu = 1007.;
	data_sigma = 68.17;
	true_mu = 1007.;
	true_sigma = 57.7;
      }
      // == 2.0 GeV
      if(beam_P_str == "2.0"){
	data_mu = 2013.;
	data_sigma = 141.02;
	true_mu = 2013.;
	true_sigma = 113.4;
      }
      
      double P_reweight = Gaussian_Reweight(data_mu, data_sigma, true_mu, true_sigma, 0., 2000., beamP_true);
      weight = P_reweight;
      if(evt.true_beam_PDG == 211){
	weight_piOnly = P_reweight;
      }
    }
    //FillHistBeam(evt, 1., "precut_noweight");
    //FillHistBeam(evt, weight, "precut_Preweight");
    //FillHistBeam(evt, weight_piOnly, "precut_Preweight_piOnly");

    bool pass_picuts = false;
    bool pass_pcuts = false;
    if (hadana.PassPiCuts(evt)) pass_picuts = true;
    if (hadana.PassPCuts(evt)) pass_pcuts = true;
    if (!pass_picuts && !pass_pcuts) continue;
    //FillHistBeam(evt, 1., "noweight");
    //FillHistBeam(evt, weight, "Preweight");
    //FillHistBeam(evt, weight_piOnly, "Preweight_piOnly");

    double beamP = evt.beam_inst_P*1000. * beamP_scale;
    if(beamP < beam_P_inst_cut_lower || beamP > beam_P_inst_cut_upper) continue;
    /*
    if(pass_picuts){
      //FillHistBeam(evt, 1., 211, "beam_window_noweight");
      FillHistBeam_Micheless(evt, 1., 211, "beam_window_noweight");
      FillHistBeam_Micheless(evt, 1., 13, "beam_window_noweight");
    }
    */

    //if(pass_pcuts) FillHistBeam(evt, 1., 2212, "beam_window_noweight");

    //FillHistBeam(evt, weight, "beam_window_Preweight");
    //FillHistBeam(evt, weight_piOnly, "beam_window_Preweight_piOnly");

    //cout << evt.run << ":" << evt.event << endl;
    vector<RecoDaughter> RecoDaughters_all = GetAllRecoDaughters(evt);
    //cout << "Before continue" << endl;
    if(RecoDaughters_all.size() < 1) continue;

    vector<RecoDaughter> pions = GetPions(RecoDaughters_all);
    vector<RecoDaughter> protons = GetProtons(RecoDaughters_all);

    if(evt.MC){
      //vector<RecoDaughter> true_pions = GetTruePions(RecoDaughters_all);
      //vector<RecoDaughter> true_protons = GetTrueProtons(RecoDaughters_all);
      //cout << "(run:event)\t(" << evt.run << ":" << evt.event << ")" << endl;
      FillHistDaughterTrue(RecoDaughters_all, evt, 211, 1., "beam_window_noweight");
      FillHistDaughterTrue(RecoDaughters_all, evt, 2212, 1., "beam_window_noweight");
     
      //FillHistDaughterTrue(RecoDaughters_all, evt, weight_piOnly, "beam_window_Preweight_piOnly");
    }
    //FillHistDaughterReco(RecoDaughters_all, evt, 1., "beam_window_noweight");
    
    RecoDaughters_all.clear();
    pions.clear();
    protons.clear();
    //cout << "(run:event)\t(" << evt.run << ":" << evt.event << ")\tdaughter_michel_score\t" << hadana.daughter_michel_score << endl;
  }
  SaveHistograms();
  //Make_dEdx_Range_Profile();
}

