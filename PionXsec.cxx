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

double PionXsec::Get_EQE(double P_pion, double cos_theta){
  double m_proton = 938.272;
  double m_neutron = 939.565;
  double m_pion = 139.57;
  double E_binding = 40.;

  double E_pion = sqrt( pow(P_pion, 2.0) + pow(m_pion, 2.0) );

  double numer = pow(m_proton, 2.0)  - pow(m_neutron - E_binding, 2.0) - pow(m_pion, 2.0) + 2.0 * (m_neutron - E_binding) * E_pion;
  double denom = 2.0 * ( m_neutron - E_binding - E_pion + P_pion * cos_theta );

  double EQE = numer / denom;

  //cout << "[PionXsec::Get_EQE]\tEQE\t" << EQE << "\tP_pion\t" << P_pion << "\tcos_theta\t" << cos_theta << endl;

  return EQE;
}

double PionXsec::Get_EQE_NC_Pion(double P_pion, double cos_theta, double E_binding, int which_sol){
  double m_proton = 938.272;
  double m_neutron = 939.565;
  double m_pion = 139.57;

  double E_pion = sqrt( pow(P_pion, 2.0) + pow(m_pion, 2.0) );

  double A = m_proton - E_binding - E_pion;
  double B = pow(m_pion, 2.) - pow(P_pion, 2.) - pow(m_proton, 2.);

  // == ax^2 + bx + c = 0
  double a = 4. * (A*A - P_pion * P_pion * cos_theta * cos_theta);
  double b = 4. * A * (A*A + B);
  double c = pow(A*A + B, 2.) + 4. * m_pion * m_pion * P_pion * P_pion * cos_theta * cos_theta;


  double numer1 = (-1.) * b;
  double numer_sqrt = sqrt(b*b - 4. * a * c);
  double denom = 2. * a;

  double EQE = (numer1 + (which_sol + 0.) * numer_sqrt ) / denom;

  return EQE;
}

double PionXsec::Gaussian_Reweight(double mu1, double sigma1, double mu2, double sigma2, double x_min, double x_max, double x){
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
  //cout << "[PionXsec::Gaussian_Reweight] true beam P : " << x << ", y1 : " << y1 << ", y2 : " << y2 << ", y1/y2 : " << y1/y2 << endl;
  if(out > 10.) return 10.;
  return y1/y2;
}

vector<RecoDaughter> PionXsec::GetAllRecoDaughters(const anavar & evt){

  // ==== Beam related variables
  TVector3 reco_beam_end(evt.reco_beam_endX, evt.reco_beam_endY, evt.reco_beam_endZ);
  TVector3 reco_unit_beam(evt.reco_beam_allTrack_trackDirX, evt.reco_beam_allTrack_trackDirY, evt.reco_beam_allTrack_trackDirZ);
  reco_unit_beam = (1. / reco_unit_beam.Mag() ) * reco_unit_beam;

  vector<RecoDaughter> out;

  //cout << "[PionXsec::GetAllRecoDaughters] evt.reco_daughter_allTrack_ID->size() : " << evt.reco_daughter_allTrack_ID->size() << endl;
  for (size_t i = 0; i < evt.reco_daughter_allTrack_ID->size(); i++){
    if((*evt.reco_daughter_allTrack_dQdX_SCE)[i].empty()) continue;

    RecoDaughter this_RecoDaughter;
    //cout << "[PionXsec::GetAllRecoDaughters] i : " << i << endl;
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

    TVector3 unit_daughter((*evt.reco_daughter_allTrack_endX).at(i) - (*evt.reco_daughter_allTrack_startX).at(i), (*evt.reco_daughter_allTrack_endY).at(i) - (*evt.reco_daughter_allTrack_startY).at(i), (*evt.reco_daughter_allTrack_endZ).at(i) - (*evt.reco_daughter_allTrack_startZ).at(i) );
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

vector<RecoDaughter> PionXsec::GetPions(const vector<RecoDaughter> in){

  vector<RecoDaughter> out;

  double cut_cos_beam = 0.95;
  double cut_beam_dist = 10.;
  double cut_trackScore = 0.5;
  double cut_emScore = 0.5;
  double cut_chi2_proton = 60.;
  int cut_Nhit = 20;
  for(unsigned int i = 0; i < in.size(); i++){
    RecoDaughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() / this_in.allTrack_Chi2_ndof();
    if(this_in.PFP_trackScore() > cut_trackScore && this_in.PFP_emScore() < cut_emScore && this_chi2 > cut_chi2_proton && this_in.PFP_nHits() > cut_Nhit && this_in.Beam_Cos() < cut_cos_beam && this_in.Beam_Dist() < cut_beam_dist){
      out.push_back(this_in);
    }
  }
  
  return out;
}

vector<RecoDaughter> PionXsec::GetProtons(const vector<RecoDaughter> in){

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

vector<RecoDaughter> PionXsec::GetTruePions(const vector<RecoDaughter> in){

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

vector<RecoDaughter> PionXsec::GetTrueProtons(const vector<RecoDaughter> in){

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

void PionXsec::BookHistograms(){
  //cout << "[PionXsec::BookHistograms] Start" << endl;
  Hist.outfile = TFile::Open(fOutputFileName.c_str(), "recreate");
}

void PionXsec::ProcessEvent(const anavar & evt){
  hadana.ProcessEvent(evt);
}

void PionXsec::FillHistBeam(const anavar & evt, double weight, TString suffix){

  double beamP = evt.beam_inst_P*1000. * beamP_scale;
  double beamKE = sqrt(pow(beamP, 2) + pow(139.57, 2)) - 139.57;
  TString pitype_str = Form("%d", hadana.pitype);

  // ==== Beam
  // == Reco
  Hist.JSFillHist(suffix, "htrack_BeamP_" + suffix + "_" + pitype_str, beamP, weight, 5000., 0., 5000.);
  Hist.JSFillHist(suffix, "htrack_PandoraSlice_" + suffix + "_" + pitype_str, hadana.PassPandoraSliceCut(evt), weight, 2., 0., 2.);
  Hist.JSFillHist(suffix, "htrack_CaloSize_" + suffix + "_" + pitype_str, hadana.PassCaloSizeCut(evt),weight, 2., 0.,2.);
  Hist.JSFillHist(suffix, "htrack_beam_dx_" + suffix + "_" + pitype_str, hadana.beam_dx, weight, 200., -10., 10.);
  Hist.JSFillHist(suffix, "htrack_beam_dy_" + suffix + "_" + pitype_str, hadana.beam_dy, weight, 200., -10.,10.);
  Hist.JSFillHist(suffix, "htrack_beam_dz_" + suffix + "_" + pitype_str, hadana.beam_dz, weight, 200., -10.,10.);
  Hist.JSFillHist(suffix, "htrack_beam_dxy_" + suffix + "_" + pitype_str, hadana.beam_dxy, weight, 200., -10.,10.);
  Hist.JSFillHist(suffix, "htrack_beam_costh_" + suffix + "_" + pitype_str, hadana.beam_costh, weight, 200., -1.,1.);
  Hist.JSFillHist(suffix, "htrack_daughter_michel_score_" + suffix + "_" + pitype_str, hadana.daughter_michel_score, weight, 100., 0.,1.);
  Hist.JSFillHist(suffix, "htrack_chi2_proton_" + suffix + "_" + pitype_str, hadana.chi2_proton, weight, 1000., 0., 100.);
  Hist.JSFillHist(suffix, "htrack_beam_inst_XY_" + suffix + "_" + pitype_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 400., -50., -10., 400., 400., 440.); 

  // == Truth
  if(evt.MC){
    double beamP_true = evt.true_beam_startP * 1000.;
    double true_KE_ff = hadana.true_ffKE;
    double true_P_ff = hadana.map_BB[abs(evt.true_beam_PDG)] -> KEtoMomentum(true_KE_ff);
    double beamP_res = (beamP_true - beamP) / beamP_true;
    Hist.JSFillHist(suffix, "htrack_BeamP_Res_" + suffix + "_" + pitype_str, beamP_res, weight, 4000, -2, 2);
    Hist.JSFillHist(suffix, "htrack_BeamP_true_" + suffix + "_" + pitype_str, beamP_true, weight, 5000, 0, 5000);
    Hist.JSFillHist(suffix, "htrack_KEffTruth_" + suffix + "_" + pitype_str, true_KE_ff, weight, 5000, 0, 5000);
    Hist.JSFillHist(suffix, "htrack_PffTruth_" + suffix + "_" + pitype_str, true_P_ff, weight, 5000, 0, 5000);
  }
}

void PionXsec::FillHistDaughterTrue(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix){

  // ==== Basic variables for daughter selection
  for(unsigned int i = 0; i < daughters.size(); i++){
    RecoDaughter this_daughter = daughters.at(i);
    TString particle_str = "";
    int this_PdgID = this_daughter.PFP_true_byHits_PDG();
    if(this_PdgID == 2212) particle_str = "proton";
    else if(abs(this_PdgID) == 211) particle_str = "pion";
    else if(abs(this_PdgID) == 13) particle_str = "muon";
    else particle_str = "other";

    // == Daughter selection variables
    double this_beam_cos = this_daughter.Beam_Cos();
    double this_beam_dist = this_daughter.Beam_Dist();
    double this_Nhits = this_daughter.PFP_nHits();
    double this_trackScore = this_daughter.PFP_trackScore();
    double this_emScore = this_daughter.PFP_emScore();
    double this_chi2 = this_daughter.allTrack_Chi2_proton() / this_daughter.allTrack_Chi2_ndof();
    double this_mean_dEdx = hadana.Truncatd_Mean_dEdx(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE());
    double this_proton_chi2 = hadana.Particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), false, 2212);
    double this_pion_chi2 = hadana.Particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), false, 211);
    double this_muon_chi2 = hadana.Particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), false, 13);
    double this_proton_chi2_offset = hadana.Particle_chi2_with_offset(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), false, 2212);
    double this_pion_chi2_offset = hadana.Particle_chi2_with_offset(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), false, 211);
    double this_muon_chi2_offset = hadana.Particle_chi2_with_offset(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), false, 13);

    //cout << "this_PdgID\t" << this_PdgID << "\tchi2_p\t" << this_proton_chi2 << "\tchi2_pi\t" << this_pion_chi2 << "\tchi2_mu\t" << this_muon_chi2 << endl;
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_beam_cos_" + suffix, this_beam_cos, weight, 200., -1., 1.);
    // == Not broken beam track
    if(evt.true_beam_ID != this_daughter.PFP_true_byHits_ID()){
      Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_beam_cos_not_broken_track_" + suffix, this_beam_cos, weight, 200., -1., 1.);
    }
    else{ // == Broken beam track
      Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_beam_cos_broken_track_" + suffix, this_beam_cos, weight, 200., -1., 1.);
    }

    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_beam_dist_" + suffix, this_beam_dist, weight, 100., 0., 100.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_trackScore_" + suffix, this_trackScore, weight, 100., 0., 1.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_emScore_" + suffix, this_emScore, weight, 100., 0., 1.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_chi2_proton_" + suffix, this_chi2, weight, 1000., 0., 100.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_PFP_nHits_" + suffix, this_Nhits, weight, 1000., 0., 1000.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_mean_dEdx_" + suffix, this_mean_dEdx, weight, 200., 0., 20.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_this_proton_chi2_" + suffix, this_proton_chi2, weight, 1000., 0., 200.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_this_pion_chi2_" + suffix, this_pion_chi2, weight, 1000., 0., 200.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_this_muon_chi2_" + suffix, this_muon_chi2, weight, 1000., 0., 200.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_this_proton_chi2_offset_" + suffix, this_proton_chi2_offset, weight, 1000., 0., 200.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_this_pion_chi2_offset_" + suffix, this_pion_chi2_offset, weight, 1000., 0., 200.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_this_muon_chi2_offset_" + suffix, this_muon_chi2_offset, weight, 1000., 0., 200.);

    /*
    double fitted_proton_chi2 = hadana.Fit_particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), false, 2212);
    double fitted_pion_chi2 = hadana.Fit_particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), false, 211);
    double fitted_muon_chi2 = hadana.Fit_particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), false, 13);
    double fitted_proton_chi2_plane1 = hadana.Fit_particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE_plane1(), this_daughter.allTrack_resRange_SCE_plane1(), false, 2212);
    double fitted_pion_chi2_plane1 = hadana.Fit_particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE_plane1(), this_daughter.allTrack_resRange_SCE_plane1(), false, 211);
    double fitted_muon_chi2_plane1 = hadana.Fit_particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE_plane1(), this_daughter.allTrack_resRange_SCE_plane1(), false, 13);
    double min_fitted_proton_chi2 = min(fitted_proton_chi2, fitted_proton_chi2_plane1);
    double max_fitted_proton_chi2 = max(fitted_proton_chi2, fitted_proton_chi2_plane1);
    double min_fitted_pion_chi2 = min(this_pion_chi2, fitted_pion_chi2_plane1);
    double max_fitted_pion_chi2 = max(this_pion_chi2, fitted_pion_chi2_plane1);

    //cout << "this_PdgID\t" << this_PdgID << "\tchi2_p\t" << min_fitted_proton_chi2 << "\tchi2_pi\t" << max_fitted_pion_chi2 << endl;
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_min_fitted_proton_chi2_" + suffix, min_fitted_proton_chi2, weight, 1000., 0., 100.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_max_fitted_proton_chi2_" + suffix, max_fitted_proton_chi2, weight, 1000., 0., 100.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_min_fitted_pion_chi2_" + suffix, min_fitted_pion_chi2, weight, 1000., 0., 100.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_max_fitted_pion_chi2_" + suffix, max_fitted_pion_chi2, weight, 1000., 0., 100.);
    */
    // == Cutflow for daughter proton & pion selection
    //if(!(this_PdgID == 2212 || abs(this_PdgID) == 211)) continue;
    double true_start_P = this_daughter.PFP_true_byHits_startP() * 1000.;
    //double true_start_KE = hadana.map_BB[abs(this_PdgID)] -> MomentumtoKE(true_start_P);
    //cout << "[PionXsec::FillHistDaughterTrue]\t" << i << "\tthis_PdgID\t" << this_PdgID << "\ttrue_start_P : " << true_start_P << "\tTrue_beam_ID\t" << evt.true_beam_ID << "\tPFP_true_byHits_ID\t" << this_daughter.PFP_true_byHits_ID() << "\treco_beam_true_byHits_ID\t" << evt.reco_beam_true_byHits_ID << endl;

    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_TrueP_nocut_" + suffix, true_start_P, weight, 5000., 0., 5000.);
    
    if(this_beam_cos < 0.95){    
      Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_TrueP_cos_" + suffix, true_start_P, weight, 5000., 0., 5000.);
      if(this_beam_dist < 10.){
	Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_TrueP_dist_" + suffix, true_start_P, weight, 5000., 0., 5000.);
	if(this_Nhits > 20){
	  Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_TrueP_Nhits_" + suffix, true_start_P, weight, 5000., 0., 5000.);
	  if(this_emScore < 0.5){
	    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_TrueP_emScore_" + suffix, true_start_P, weight, 5000., 0., 5000.);
	    if(this_trackScore > 0.5){
	      Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_TrueP_TrackScore_" + suffix, true_start_P, weight, 5000., 0., 5000.);
	      Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_TrueP_vs_proton_chi2_" + suffix, true_start_P, this_proton_chi2, weight, 5000., 0.,5000., 1000., 0, 1000.);
	      Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_TrueP_vs_proton_chi2_offset_" + suffix, true_start_P, this_proton_chi2_offset, weight, 5000., 0.,5000., 1000., 0, 1000.);

	      // == Proton ID
	      if(this_chi2 < 50.){
		Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_TrueP_chi2_proton_protonID_" + suffix, true_start_P, weight, 5000., 0., 5000.);
	      }

	      // == Pion ID
	      if(this_chi2 > 60.){
		Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_TrueP_chi2_proton_pionID_" + suffix, true_start_P, weight, 5000., 0., 5000.);
		if(this_mean_dEdx < 6.){
		  Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_TrueP_mean_dEdx_pionID_" + suffix, true_start_P, weight, 5000., 0., 5000.);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  // ==== Matching between Reco daughters and true daughters
  //if ((*evt.reco_daughter_PFP_true_byHits_ID)[i] == (*evt.true_beam_daughter_ID)[j]){
}

void PionXsec::FillHistDaughterPurity(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix){

  vector<RecoDaughter> pions = GetPions(daughters);
  vector<RecoDaughter> protons = GetProtons(daughters);

  //cout << "[PionXsec::FillHistDaughterPurity] Start" << endl;
  // ==== Loop for pions
  for(unsigned int i = 0; i < pions.size(); i++){
    RecoDaughter this_daughter = pions.at(i);
    int this_PdgID = this_daughter.PFP_true_byHits_PDG();
    TString particle_str = "";
    if(this_PdgID == 2212) particle_str = "proton";
    else if(abs(this_PdgID) == 211) particle_str = "pion";
    else if(abs(this_PdgID) == 13) particle_str = "muon";
    else particle_str = "other";
    double true_start_P = this_daughter.PFP_true_byHits_startP() * 1000.;
    //double true_start_KE = hadana.map_BB[abs(this_PdgID)] -> MomentumtoKE(true_start_P);

    Hist.JSFillHist(suffix, "hdaughter_pion_purity_TrueP_all_" + suffix, true_start_P, weight, 5000., 0., 5000.);
    Hist.JSFillHist(suffix, "hdaughter_pion_purity_TrueP_" + particle_str + "_" + suffix, true_start_P, weight, 5000., 0., 5000.);
  }

  // ==== Loop for protons
  for(unsigned int i = 0; i < protons.size(); i++){
    RecoDaughter this_daughter = protons.at(i);
    int this_PdgID = this_daughter.PFP_true_byHits_PDG();
    TString particle_str = "";
    if(this_PdgID == 2212) particle_str = "proton";
    else if(abs(this_PdgID) == 211) particle_str = "pion";
    else if(abs(this_PdgID) == 13) particle_str = "muon";
    else particle_str = "other";
    double true_start_P = this_daughter.PFP_true_byHits_startP() * 1000.;
    //double true_start_KE = hadana.map_BB[abs(this_PdgID)] -> MomentumtoKE(true_start_P);

    Hist.JSFillHist(suffix, "hdaughter_proton_purity_TrueP_all_"+ suffix, true_start_P, weight, 5000., 0., 5000.);
    Hist.JSFillHist(suffix, "hdaughter_proton_purity_TrueP_" + particle_str + "_" + suffix, true_start_P, weight, 5000., 0., 5000.);
  }

}

void PionXsec::FillHistQE_MCstudy(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix){

  // ==== Beam related variables
  TVector3 reco_unit_beam(evt.reco_beam_allTrack_trackDirX, evt.reco_beam_allTrack_trackDirY, evt.reco_beam_allTrack_trackDirZ);
  reco_unit_beam = (1. / reco_unit_beam.Mag() ) * reco_unit_beam;
  double true_beam_End_E = evt.reco_beam_true_byHits_endE;

  // ==== Draw truth level distributions for QE using true_beam_daughter
  for(unsigned int i = 0; i < (*evt.true_beam_daughter_ID).size(); i++){
    if(evt.true_beam_ID == (*evt.true_beam_daughter_ID).at(i)) continue;
    //if((*evt.true_beam_daughter_PDG).at(i) == 211){ // == only for pion+
    int this_PdgID = (*evt.true_beam_daughter_PDG).at(i);
    TString particle_str = "";
    if(this_PdgID == 2212) particle_str = "proton";
    else if(abs(this_PdgID) == 211) particle_str = "pion";
    else if(abs(this_PdgID) == 13) particle_str = "muon";
    else particle_str = "other";
    TVector3 unit_daughter((*evt.true_beam_daughter_startPx).at(i), (*evt.true_beam_daughter_startPy).at(i), (*evt.true_beam_daughter_startPz).at(i) );
    unit_daughter = (1./ unit_daughter.Mag() ) * unit_daughter;

    TVector3 unit_beam(evt.true_beam_endPx, evt.true_beam_endPy, evt.true_beam_endPz);
    unit_beam = (1. / unit_beam.Mag() ) * unit_beam;

    double cos_theta = cos(unit_beam.Angle(unit_daughter));
    double P_beam = evt.true_beam_endP * 1000.;
    double m_beam = evt.true_beam_mass;
    double E_beam = sqrt(P_beam * P_beam + m_beam * m_beam);
    double P_daughter = (*evt.true_beam_daughter_startP).at(i) * 1000.;
    double EQE = Get_EQE(P_daughter, cos_theta);
    double EQE_NC_40 = Get_EQE_NC_Pion(P_daughter, cos_theta, 40., -1.);
    double EQE_NC_4 = Get_EQE_NC_Pion(P_daughter, cos_theta, 4., -1.);
    /*
      cout << "[PionXsec::FillHistQE_MCstudy] cos_theta : " << cos_theta << ", P_daughter : " << P_daughter << ", EQE_NC_4 : "
      << EQE_NC_4 << ", E_beam : " << E_beam << ", mEQE_NC_4 : " << EQE_NC_4 - E_beam << endl;
    */
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_EQE_" + suffix, EQE, weight, 6000., -3000., 3000.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_EQE_NC_40_" + suffix, EQE_NC_40, weight, 6000., -3000., 3000.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_EQE_NC_4_" + suffix, EQE_NC_4, weight, 6000., -3000., 3000.);
   
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_EQEmE_" + suffix, EQE - E_beam, weight, 6000., -3000., 3000.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_EQEmE_NC_40_" + suffix, EQE_NC_40 - E_beam, weight, 6000., -3000., 3000.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_EQEmE_NC_4_" + suffix, EQE_NC_4 - E_beam, weight, 6000., -3000., 3000.);

    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_E_beam_vs_EQE_" + suffix, E_beam, EQE, 1., 240., 100., 1300., 240., 100., 1300.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_E_beam_vs_EQE_NC_40_" + suffix, E_beam, EQE_NC_40, 1., 240., 100., 1300., 240., 100., 1300.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_E_beam_vs_EQE_NC_4_" + suffix, E_beam, EQE_NC_4, 1., 240., 100., 1300., 240., 100., 1300.);

  }
}

void PionXsec::FillHistQE_Reco(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix){
 
  TString pitype_str = Form("%d", hadana.pitype);
  // ==== Beam related variables 
  TVector3 reco_unit_beam(evt.reco_beam_allTrack_trackDirX, evt.reco_beam_allTrack_trackDirY, evt.reco_beam_allTrack_trackDirZ);
  reco_unit_beam = (1. / reco_unit_beam.Mag() ) * reco_unit_beam;
  // == Beam kinematics under pion+ assumption
  double beamP = evt.beam_inst_P*1000. * beamP_scale;
  double beamKE = sqrt(pow(beamP, 2) + pow(139.57, 2)) - 139.57;
  double beamKE_ff = beamKE - 10.;
  double beamKE_end = hadana.map_BB[211]->KEAtLength(beamKE_ff, evt.reco_beam_alt_len);
  double beamE_end = beamKE_end + 139.57;

  vector<RecoDaughter> pions = GetPions(daughters);

  //if(pions.size() < 1) return;

  int N_pion = 0;
  int N_pion_michel = 0;
  int N_pion_fitted = 0;
  for(unsigned int i = 0; i < pions.size(); i++){
    RecoDaughter this_pion = pions.at(i);
    double this_michelscore = this_pion.allTrack_vertex_michel_score() / this_pion.allTrack_vertex_nHits();
    TVector3 unit_daughter( this_pion.allTrack_endX() - this_pion.allTrack_startX(),
                            this_pion.allTrack_endY() - this_pion.allTrack_startY(),
                            this_pion.allTrack_endZ() - this_pion.allTrack_startZ());
    unit_daughter = (1. / unit_daughter.Mag() ) * unit_daughter;
    double cos_theta = cos(reco_unit_beam.Angle(unit_daughter));

    if(this_michelscore > 0.5){
      double KE_daughter = hadana.map_BB[211]->KEFromRangeSpline(this_pion.allTrack_alt_len());
      double P_daughter = hadana.map_BB[211]->KEtoMomentum(KE_daughter);
      double EQE = Get_EQE(P_daughter, cos_theta);
      double EQE_NC_40 = Get_EQE_NC_Pion(P_daughter, cos_theta, 40., -1.);
      double EQE_NC_4 = Get_EQE_NC_Pion(P_daughter, cos_theta, 4., -1.);
      double mEQE = EQE - beamE_end;
      double mEQE_NC_40 = EQE_NC_40 - beamE_end;
      double mEQE_NC_4 = EQE_NC_4 - beamE_end;
     
      // == Basic kinematics for reco pion
      Hist.JSFillHist(suffix, "hdaughter_reco_pion_KE_" + suffix + "_" + pitype_str, KE_daughter, weight, 3000., 0., 3000.);
      Hist.JSFillHist(suffix, "hdaughter_reco_pion_cos_" + suffix + "_" + pitype_str, cos_theta, weight, 200., -1., 1.);
      Hist.JSFillHist(suffix, "hdaughter_reco_pion_length_" + suffix + "_" + pitype_str, this_pion.allTrack_alt_len(), weight, 700., 0., 700.);

      // == QE variables
      Hist.JSFillHist(suffix, "hdaughter_reco_EQEmE_" + suffix + "_" + pitype_str, mEQE, weight, 6000., -3000., 3000.);
      Hist.JSFillHist(suffix, "hdaughter_reco_EQEmE_NC_40_" + suffix + "_" + pitype_str, mEQE_NC_40, weight, 6000., -3000., 3000.);
      Hist.JSFillHist(suffix, "hdaughter_reco_EQEmE_NC_4_" + suffix + "_" + pitype_str, mEQE_NC_4, weight, 6000., -3000., 3000.);
      Hist.JSFillHist(suffix, "hdaughter_reco_EQEmE_NC_4_vs_KE_daughter_" + suffix + "_" + pitype_str, mEQE_NC_4, KE_daughter, weight, 1500., -1000., 500., 1000., 0., 1000.);

      if(KE_daughter > 400.){
	Hist.JSFillHist(suffix, "hdaughter_reco_pion_KE_KE400_" + suffix + "_" + pitype_str, KE_daughter, weight, 3000., 0., 3000.);
	Hist.JSFillHist(suffix, "hdaughter_reco_pion_cos_KE400_" + suffix + "_" + pitype_str, cos_theta, weight, 200., -1., 1.);
	Hist.JSFillHist(suffix, "hdaughter_reco_pion_length_KE400_" + suffix + "_" + pitype_str, this_pion.allTrack_alt_len(), weight, 700., 0., 700.);

	Hist.JSFillHist(suffix, "hdaughter_reco_EQEmE_KE400_" + suffix + "_" + pitype_str, mEQE, weight, 6000., -3000., 3000.);
	Hist.JSFillHist(suffix, "hdaughter_reco_EQEmE_NC_40_KE400_" + suffix + "_" + pitype_str, mEQE_NC_40, weight, 6000., -3000., 3000.);
	Hist.JSFillHist(suffix, "hdaughter_reco_EQEmE_NC_4_KE400_" + suffix + "_" + pitype_str, mEQE_NC_4, weight, 6000., -3000., 3000.);
      }

      N_pion_michel++;
    }
  }

  Hist.JSFillHist(suffix, "hdaughter_reco_QE_N_pion_michel_" + suffix + "_" + pitype_str, N_pion_michel, weight, 10., -0.5, 9.5);
  
}

void PionXsec::SaveHistograms(){
  
  //outputFile->cd();
  //outputFile->Write();
  Hist.WriteHist();
}

void PionXsec::Run(anavar & evt, Long64_t nentries=-1){

  BookHistograms();
  TString beam_P_str = "1.0";
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
    cout << "[PionXsec::Run] Not Valid Beam Momentum" << endl;
    return;
  }
  cout << "[PionXsec::Run] Beam Momentum is " << beam_P_str << " GeV" << endl;


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
    FillHistBeam(evt, 1., "precut_noweight");
    FillHistBeam(evt, weight, "precut_Preweight");
    FillHistBeam(evt, weight_piOnly, "precut_Preweight_piOnly");

    if (!hadana.PassPiCuts(evt)) continue;
    FillHistBeam(evt, 1., "noweight");
    FillHistBeam(evt, weight, "Preweight");
    FillHistBeam(evt, weight_piOnly, "Preweight_piOnly");

    double beamP = evt.beam_inst_P*1000. * beamP_scale;
    if(beamP < beam_P_inst_cut_lower || beamP > beam_P_inst_cut_upper) continue;
    FillHistBeam(evt, 1., "beam_window_noweight");
    FillHistBeam(evt, weight, "beam_window_Preweight");
    FillHistBeam(evt, weight_piOnly, "beam_window_Preweight_piOnly");

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
      //FillHistDaughterTrue(RecoDaughters_all, evt, 1., "noweight");
      FillHistDaughterTrue(RecoDaughters_all, evt, weight_piOnly, "Preweight_piOnly");

      //FillHistDaughterPurity(RecoDaughters_all, evt, 1., "noweight");
      FillHistDaughterPurity(RecoDaughters_all, evt, weight_piOnly, "Preweight_piOnly");

      //FillHistQE_MCstudy(RecoDaughters_all, evt, 1., "noweight");
      FillHistQE_MCstudy(RecoDaughters_all, evt, weight_piOnly, "Preweight_piOnly");
    }
    FillHistQE_Reco(RecoDaughters_all, evt, weight_piOnly, "Preweight_piOnly");

    RecoDaughters_all.clear();
    pions.clear();
    protons.clear();
    //cout << "(run:event)\t(" << evt.run << ":" << evt.event << ")\tdaughter_michel_score\t" << hadana.daughter_michel_score << endl;
  }
  SaveHistograms();
  //Make_dEdx_Range_Profile();
}

