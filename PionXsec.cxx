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

bool PionXsec::Is_EQE(const anavar & evt, double window){
  bool is_QE = false;

  TVector3 unit_beam(evt.true_beam_endPx, evt.true_beam_endPy, evt.true_beam_endPz);
  double P_beam = evt.true_beam_endP * 1000.;
  double m_beam = evt.true_beam_mass;
  double E_beam = sqrt(P_beam * P_beam + m_beam * m_beam);

  int N_daugther_piplus = 0;
  for(unsigned int i = 0; i < (*evt.true_beam_daughter_ID).size(); i++){
    int this_daughter_PID = (*evt.true_beam_daughter_PDG).at(i);
    if(this_daughter_PID == 211){
      N_daugther_piplus++;

      TVector3 unit_daughter((*evt.true_beam_daughter_startPx).at(i), (*evt.true_beam_daughter_startPy).at(i), (*evt.true_beam_daughter_startPz).at(i) );
      unit_daughter = (1./ unit_daughter.Mag() ) * unit_daughter;
      double cos_theta = cos(unit_beam.Angle(unit_daughter));
      double P_daughter = (*evt.true_beam_daughter_startP).at(i) * 1000.;
      double EQE_NC_4 = Get_EQE_NC_Pion(P_daughter, cos_theta, 4., -1.);
      double this_EQEmE = EQE_NC_4 - E_beam;
      if(this_EQEmE > (-1.) * window && this_EQEmE < window){
        is_QE = true;
      }
    }
  }

  
  return is_QE && (N_daugther_piplus == 1) && (hadana.pitype == 1);
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

bool PionXsec::Pass_Beam_PID(const anavar & evt, int PID){

  // == Study only proton and pion/muon beams
  if(PID != 2212 && PID != 211 && abs(PID) != 13) return false;

  vector<int> PDGID_candiate_vec;
  PDGID_candiate_vec.clear();
  if(PID == 2212){
    PDGID_candiate_vec.push_back(2212);
  }

  if(PID == 211 || PID == 13){
    PDGID_candiate_vec.push_back(211);
    PDGID_candiate_vec.push_back(13);
    PDGID_candiate_vec.push_back(-13);
  }

  if (evt.reco_reconstructable_beam_event == 0) return false; // remove empty events first
  if (evt.MC){
    for (size_t i = 0; i < PDGID_candiate_vec.size(); ++i){
      if (evt.true_beam_PDG == PDGID_candiate_vec[i]) return true; // truth matched
    }
    return false;
  }
  else{ // real data
    if (evt.beam_inst_trigger == 8) return false; // is cosmics
    if (evt.beam_inst_nMomenta != 1 || evt.beam_inst_nTracks != 1) return false;
    for (size_t i = 0; i<PDGID_candiate_vec.size(); ++i){
      for (size_t j = 0; j<evt.beam_inst_PDG_candidates->size(); ++j){
        if ((*evt.beam_inst_PDG_candidates)[j] == PDGID_candiate_vec[i]) return true; // Matching data PID result and PDGID candidates
      }
    }
    return false;
  }
}

double PionXsec::Convert_P_Spectrometer_to_P_ff(const anavar & evt, double P_beam_inst, TString suffix, TString key, int syst){

  double out = P_beam_inst;
  double delta_P = 0.;
  double p0 = 0., p1 = 0., p2 = 0.;

  if(suffix.Contains("pion")){
    if(key == "AllTrue"){
      if(syst == 0){
        p0 = 182.7;
        p1 = -0.4893;
        p2 = 0.0003186;
      }
      else if(syst == -1){
        p0 = 140.2;
        p1 = -0.3979;
        p2 = 0.0002678;
      }
      else if(syst == 1){
        p0 = 235.2;
        p1 = -0.6022;
        p2 = 0.0003807;
      }
      else{
        return out;
      }
    }
  }
  else{
    return out;
  }

  delta_P = p0 + p1 * out + p2 * out * out;
  
  //cout << "[PionXsec::Convert_P_Spectrometer_to_P_ff] P_beam_inst : " << P_beam_inst << ", suffix : " << suffix << ", key : " << key << endl;

  return out - delta_P;
}

void PionXsec::Fill_Eslice_Histograms(TString histname, TString suffix, double KE_init, double KE_end, double KE_int, double weight, const vector<double> binning, int N_bin, bool fill_int){

  // == Fill KE_init, KE_inc, and KE_end distributions
  int i_init = -1;
  int i_end = -1;
  int i_int = -1;
  for(unsigned int i_bin = 0; i_bin < binning.size() - 1; i_bin++){
    if(KE_init > binning.at(i_bin) && KE_init < binning.at(i_bin + 1)) i_init = i_bin;
    if(KE_end > binning.at(i_bin) && KE_end < binning.at(i_bin + 1)) i_end = i_bin;
    if(KE_int > binning.at(i_bin) && KE_int < binning.at(i_bin + 1)) i_int = i_bin;
  }

  if(i_init == -1 || i_end == -1 || i_init == i_end) return;

  if(fill_int && i_int > 0){
    // == Fill KE_int
    Hist.JSFillHist(suffix, "Xsec_KE_int_" + histname + "_" + suffix, binning.at(i_int) + 0.1, weight, N_bin, &binning[0]);
  }
  else{
    // == Fill KE_init
    Hist.JSFillHist(suffix, "Xsec_KE_init_" + histname + "_" + suffix, binning.at(i_init - 1) + 0.1, weight, N_bin, &binning[0]);

    // == Fill KE_end
    Hist.JSFillHist(suffix, "Xsec_KE_end_" + histname + "_" + suffix, binning.at(i_end) + 0.1, weight, N_bin, &binning[0]);

    // == Fill KE_inc
    for(unsigned int i_bin = i_end; i_bin < i_init; i_bin++){
      Hist.JSFillHist(suffix, "Xsec_KE_inc_" + histname + "_" + suffix, binning.at(i_bin) + 0.1, weight, N_bin, &binning[0]);
    }
  }

  return;
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

vector<RecoDaughter> PionXsec::GetPions(const vector<RecoDaughter> in){

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

  //cout << "[PionXsec::FillHistBeam] " << suffix << endl;
  double P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  double mass_beam = 139.57;
  double KE_beam_inst = sqrt(pow(P_beam_inst, 2) + pow(mass_beam, 2)) - mass_beam;

  double P_ff_reco = Convert_P_Spectrometer_to_P_ff(evt, P_beam_inst, suffix, "AllTrue", 0);
  //cout << "[PionXsec::FillHistBeam] P_ff_reco : " << P_ff_reco << endl;
  double KE_ff_reco = sqrt(pow(P_ff_reco, 2) + pow(mass_beam, 2)) - mass_beam;
  //cout << "[PionXsec::FillHistBeam] KE_ff_reco : " << KE_ff_reco << endl;
  double KE_end_reco = hadana.map_BB[211]->KEAtLength(KE_ff_reco, evt.reco_beam_alt_len);
  //cout << "[PionXsec::FillHistBeam] KE_end_reco : " << KE_end_reco << endl;
  double E_end_reco = KE_end_reco + mass_beam;
  //cout << "[PionXsec::FillHistBeam] P_beam_inst : " << P_beam_inst << ", P_ff_reco : " << P_ff_reco << ", KE_ff_reco : " << KE_ff_reco << endl;

  // == Angles
  TVector3 unit_TPC_beam(evt.reco_beam_allTrack_trackDirX, evt.reco_beam_allTrack_trackDirY, evt.reco_beam_allTrack_trackDirZ);
  unit_TPC_beam = (1. / unit_TPC_beam.Mag() ) * unit_TPC_beam;
  TVector3 unit_spec_beam(evt.beam_inst_dirX, evt.beam_inst_dirY, evt.beam_inst_dirZ);
  unit_spec_beam = (1. / unit_spec_beam.Mag() ) * unit_spec_beam; 
  double cos_theta_spec_TPC = unit_TPC_beam.Dot(unit_spec_beam);
  
  TVector3 unit_TPC_yz(0., evt.reco_beam_allTrack_trackDirY, evt.reco_beam_allTrack_trackDirZ);
  unit_TPC_yz = (1. / unit_TPC_yz.Mag() ) * unit_TPC_yz;
  TVector3 unit_spec_yz(0., evt.beam_inst_dirY, evt.beam_inst_dirZ);
  unit_spec_yz = (1. / unit_spec_yz.Mag() ) * unit_spec_yz;
  double cos_phi_spec_TPC = unit_TPC_yz.Dot(unit_spec_yz);
  double phi_spec_TPC = unit_TPC_yz.Angle(unit_spec_yz) * (unit_TPC_yz.Cross(unit_spec_yz)).X() / fabs((unit_TPC_yz.Cross(unit_spec_yz)).X());

  TVector3 unit_TPC_xz(evt.reco_beam_allTrack_trackDirX, 0., evt.reco_beam_allTrack_trackDirZ);
  unit_TPC_xz = (1. / unit_TPC_xz.Mag() ) * unit_TPC_xz;
  TVector3 unit_spec_xz(evt.beam_inst_dirX, 0., evt.beam_inst_dirZ);
  unit_spec_xz = (1. / unit_spec_xz.Mag() ) * unit_spec_xz;
  double cos_azim_spec_TPC = unit_TPC_xz.Dot(unit_spec_xz);
  double azim_spec_TPC = unit_TPC_xz.Angle(unit_spec_xz) * (unit_TPC_xz.Cross(unit_spec_xz)).Y() / fabs((unit_TPC_xz.Cross(unit_spec_xz)).Y());

  //if(fabs(azim_spec_TPC > 2.)) cout << "[PionXsec::FillHistBeam] cos_theta_spec_TPC : " << cos_theta_spec_TPC << ", phi_spec_TPC : " << phi_spec_TPC << ", azim_spec_TPC : " << azim_spec_TPC << endl;

  TString pitype_str = Form("%d", hadana.pitype);

  TString Is_QE_str = "";
  if(Is_QE_true) Is_QE_str = "_IsQE";

  //cout << "hadana.chi2_proton : " << hadana.chi2_proton << ", SB chi2 : " << hadana.Particle_chi2( (*evt.reco_beam_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 2212) << endl;
  //cout << "hadana.chi2_proton : " << hadana.chi2_proton << ", SB chi2 : " << hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 2212) << endl;
  double chi2_0p90_res = hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 2212, 0.9);
  double chi2_0p95_res = hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 2212, 0.95);
  double chi2_1p05_res = hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 2212, 1.05);
  double chi2_1p10_res = hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 2212, 1.10);

  double chi2_proton = hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 2212);
  double chi2_pion = hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 211);
  double chi2_muon = hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 13);
  double chi2_diff_pion_proton = chi2_pion - chi2_proton;

  // ==== Beam
  // == Reco
  Hist.JSFillHist(suffix, "htrack_P_beam_inst_" + suffix + "_" + pitype_str + Is_QE_str, P_beam_inst, weight, 5000., 0., 5000.);
  Hist.JSFillHist(suffix, "htrack_PandoraSlice_" + suffix + "_" + pitype_str + Is_QE_str, hadana.PassPandoraSliceCut(evt), weight, 2., 0., 2.);
  Hist.JSFillHist(suffix, "htrack_CaloSize_" + suffix + "_" + pitype_str + Is_QE_str, hadana.PassCaloSizeCut(evt),weight, 2., 0.,2.);
  Hist.JSFillHist(suffix, "htrack_beam_dx_" + suffix + "_" + pitype_str + Is_QE_str, hadana.beam_dx, weight, 200., -10., 10.);
  Hist.JSFillHist(suffix, "htrack_beam_dy_" + suffix + "_" + pitype_str + Is_QE_str, hadana.beam_dy, weight, 200., -10., 10.);
  Hist.JSFillHist(suffix, "htrack_beam_dz_" + suffix + "_" + pitype_str + Is_QE_str, hadana.beam_dz, weight, 200., -10., 10.);
  Hist.JSFillHist(suffix, "htrack_beam_dxy_" + suffix + "_" + pitype_str + Is_QE_str, hadana.beam_dxy, weight, 200., -10., 10.);
  Hist.JSFillHist(suffix, "htrack_beam_costh_" + suffix + "_" + pitype_str + Is_QE_str, hadana.beam_costh, weight, 2000., -1., 1.);
  Hist.JSFillHist(suffix, "htrack_beam_phi_" + suffix + "_" + pitype_str + Is_QE_str, phi_spec_TPC, weight, 800., -4., 4.);
  Hist.JSFillHist(suffix, "htrack_beam_azim_" + suffix + "_" + pitype_str + Is_QE_str, azim_spec_TPC, weight, 800., -4., 4.);
  
  Hist.JSFillHist(suffix, "htrack_daughter_michel_score_" + suffix + "_" + pitype_str + Is_QE_str, hadana.daughter_michel_score, weight, 100., 0.,1.);

  Hist.JSFillHist(suffix, "htrack_chi2_proton_" + suffix + "_" + pitype_str + Is_QE_str, hadana.chi2_proton, weight, 1000., 0., 1000.);
  Hist.JSFillHist(suffix, "htrack_chi2_pion_" + suffix + "_" + pitype_str + Is_QE_str, chi2_pion, weight, 1000., 0., 1000.);
  Hist.JSFillHist(suffix, "htrack_chi2_muon_" + suffix + "_" + pitype_str + Is_QE_str, chi2_muon, weight, 1000., 0., 1000.);
  Hist.JSFillHist(suffix, "htrack_chi2_diff_pion_proton_" + suffix + "_" + pitype_str + Is_QE_str, chi2_diff_pion_proton, weight, 2000., -1000., 1000.);

  Hist.JSFillHist(suffix, "htrack_chi2_proton_0p90_res_" + suffix + "_" + pitype_str + Is_QE_str, chi2_0p90_res, weight, 1000., 0., 1000.);
  Hist.JSFillHist(suffix, "htrack_chi2_proton_0p95_res_" + suffix + "_" + pitype_str + Is_QE_str, chi2_0p95_res, weight, 1000., 0., 1000.);
  Hist.JSFillHist(suffix, "htrack_chi2_proton_1p05_res_" + suffix + "_" + pitype_str + Is_QE_str, chi2_1p05_res, weight, 1000., 0., 1000.);
  Hist.JSFillHist(suffix, "htrack_chi2_proton_1p10_res_" + suffix + "_" + pitype_str + Is_QE_str, chi2_1p10_res, weight, 1000., 0., 1000.);

  Hist.JSFillHist(suffix, "htrack_beam_alt_len_" + suffix + "_" + pitype_str + Is_QE_str, evt.reco_beam_alt_len, weight, 600., 0., 600.);
  Hist.JSFillHist(suffix, "htrack_beam_inst_XY_" + suffix + "_" + pitype_str + Is_QE_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 400., -50., -10., 400., 400., 440.); 
  Hist.JSFillHist(suffix, "htrack_beam_start_XY_" + suffix + "_" + pitype_str + Is_QE_str, evt.reco_beam_calo_startX, evt.reco_beam_calo_startY, weight, 400., -50., -10., 400., 400., 440.);
  Hist.JSFillHist(suffix, "htrack_beam_start_Z_" + suffix + "_" + pitype_str + Is_QE_str, evt.reco_beam_calo_startZ, weight, 1000., -100., 900.);
  Hist.JSFillHist(suffix, "htrack_beam_end_Z_" + suffix + "_" + pitype_str + Is_QE_str, evt.reco_beam_calo_endZ, weight, 1000., -100., 900.);

  Hist.JSFillHist(suffix, "htrack_P_ff_reco_" + suffix + "_" + pitype_str + Is_QE_str, P_ff_reco, weight, 5000., 0., 5000.);
  Hist.JSFillHist(suffix, "htrack_KE_ff_reco_" + suffix + "_" + pitype_str + Is_QE_str, KE_ff_reco, weight, 5000., 0., 5000.);
  Hist.JSFillHist(suffix, "htrack_KE_end_reco_" + suffix + "_" + pitype_str + Is_QE_str, KE_end_reco, weight, 5000., 0., 5000.);


  // == Comaprison between Beam Inst. and TPC track
  TVector3 Beam_Inst_r(evt.beam_inst_X, evt.beam_inst_Y, evt.beam_inst_Z);
  TVector3 Beam_Inst_n(evt.beam_inst_dirX, evt.beam_inst_dirY, evt.beam_inst_dirZ);
  TVector3 TPC_r(evt.reco_beam_startX, evt.reco_beam_startY, evt.reco_beam_startZ);
  TVector3 TPC_n(evt.reco_beam_trackDirX, evt.reco_beam_trackDirY, evt.reco_beam_trackDirZ);
  double distance = fabs( (Beam_Inst_n.Cross(TPC_n)).Dot(Beam_Inst_r - TPC_r) ) / ( (Beam_Inst_n.Cross(TPC_n)).Mag() );
  //cout << "[PionXsec::FillHistBeam] " << suffix << ", Beam_Inst_n.Mag() : " << Beam_Inst_n.Mag() << ", TPC_n.Mag() : " << TPC_n.Mag() << ", Dist. : " << distance << endl;

  int i_KE_50 =KE_beam_inst / 50.;
  int i_KE_100 = KE_beam_inst / 100.;
  int i_P_50 = P_beam_inst / 50.;
  int i_P_100 = P_beam_inst / 100.;
  TString KE_range_str_50MeV = Form("KE_beam_inst%dto%dMeV_", 50 * i_KE_50, 50 * (i_KE_50  + 1));
  TString KE_range_str_100MeV = Form("KE_beam_inst%dto%dMeV_", 100 * i_KE_100, 100 * (i_KE_100 + 1));
  TString P_range_str_50MeV = Form("P_beam_inst%dto%dMeV_", 50 * i_P_50, 50 * (i_P_50 + 1));
  TString P_range_str_100MeV = Form("P_beam_inst%dto%dMeV_", 100 * i_P_100, 100 * (i_P_100 + 1));


  // == Truth
  if(evt.MC){
    double P_beam_inst_true = evt.true_beam_startP * 1000.;
    double KE_beam_inst_true = sqrt(pow(P_beam_inst_true, 2) + pow(mass_beam, 2)) - mass_beam;

    double true_KE_ff = hadana.true_ffKE;
    double true_P_ff = hadana.map_BB[abs(evt.true_beam_PDG)] -> KEtoMomentum(true_KE_ff);

    double KE_diff = KE_beam_inst - true_KE_ff;
    double KE_diff_true = KE_beam_inst_true - true_KE_ff;
    double KE_diff_beam_true = KE_beam_inst - KE_beam_inst_true;

    double diff_KE_ff_reco_KE_true = KE_ff_reco - true_KE_ff;
    double diff_P_ff_reco_P_true = P_ff_reco - true_P_ff;

    TString true_beam_PDG_str = Form("%d", abs(evt.true_beam_PDG));
    double P_beam_inst_res = (P_beam_inst_true - P_beam_inst) / P_beam_inst_true;
    Hist.JSFillHist(suffix, "htrack_P_beam_inst_Res_" + suffix + "_" + pitype_str, P_beam_inst_res, weight, 4000, -2, 2);
    Hist.JSFillHist(suffix, "htrack_P_beam_inst_true_" + suffix + "_" + pitype_str, P_beam_inst_true, weight, 5000, 0, 5000);
    Hist.JSFillHist(suffix, "htrack_P_beam_inst_true_" + suffix + "_PID" + true_beam_PDG_str, P_beam_inst_true, weight, 5000, 0, 5000);
    
    Hist.JSFillHist(suffix, "htrack_KEffTruth_" + suffix + "_" + pitype_str, true_KE_ff, weight, 5000, 0, 5000);
    Hist.JSFillHist(suffix, "htrack_PffTruth_" + suffix + "_" + pitype_str, true_P_ff, weight, 5000, 0, 5000);
    
    Hist.JSFillHist(suffix, "htrack_diff_KE_ff_reco_KE_true_" + suffix + "_" + pitype_str, diff_KE_ff_reco_KE_true, weight, 200., -100., 100.);
    Hist.JSFillHist(suffix, "htrack_diff_P_ff_reco_P_true_" + suffix + "_" + pitype_str, diff_P_ff_reco_P_true, weight, 200., -100., 100.);
    Hist.JSFillHist(suffix, "htrack_diff_KE_ff_reco_KE_true_" + KE_range_str_50MeV + suffix + "_" + pitype_str, diff_KE_ff_reco_KE_true, weight, 200., -100., 100.);
    Hist.JSFillHist(suffix, "htrack_diff_KE_ff_reco_KE_true_" + KE_range_str_100MeV + suffix + "_" + pitype_str, diff_KE_ff_reco_KE_true, weight, 200., -100., 100.);
    Hist.JSFillHist(suffix, "htrack_diff_P_ff_reco_P_true_" + P_range_str_50MeV + suffix + "_" + pitype_str, diff_P_ff_reco_P_true, weight, 200., -100., 100.);
    Hist.JSFillHist(suffix, "htrack_diff_P_ff_reco_P_true_" + P_range_str_100MeV + suffix + "_" + pitype_str, diff_P_ff_reco_P_true, weight, 200., -100., 100.);

    Hist.JSFillHist(suffix, "htrack_diff_P_ff_reco_P_true_" + suffix + "_" + pitype_str, diff_P_ff_reco_P_true, weight, 200., -100., 100.);
    Hist.JSFillHist(suffix, "htrack_diff_KE_ff_reco_KE_true_" +suffix + "_PID" + true_beam_PDG_str, diff_KE_ff_reco_KE_true, weight, 200., -100., 100.);
    Hist.JSFillHist(suffix, "htrack_diff_P_ff_reco_P_true_" + suffix + "_PID" + true_beam_PDG_str, diff_P_ff_reco_P_true, weight, 200., -100., 100.);

    Hist.JSFillHist(suffix, "htrack_P_beam_inst_true_vs_D_" + suffix + "_PID" + true_beam_PDG_str, P_beam_inst_true, distance, weight, 300., 0., 1500., 100., 0., 100.);
    Hist.JSFillHist(suffix, "htrack_BeamKE_vs_KEffTruth_" + suffix + "_" + pitype_str, KE_beam_inst, true_KE_ff, weight, 300., 0., 1500., 300., 0., 1500.);
    Hist.JSFillHist(suffix, "htrack_BeamKE_vs_KEffTruth_" + suffix + "_PID" + true_beam_PDG_str, KE_beam_inst, true_KE_ff, weight, 300., 0., 1500., 300., 0., 1500.);
    Hist.JSFillHist(suffix, "htrack_BeamKE_vs_KE_diff_" + suffix + "_" + pitype_str, KE_beam_inst, KE_diff, weight, 300., 0., 1500., 500., -500., 500.);
    Hist.JSFillHist(suffix, "htrack_BeamKE_vs_KE_diff_" + suffix + "_PID" + true_beam_PDG_str, KE_beam_inst, KE_diff, weight, 300., 0., 1500., 500., -500., 500.);
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

void PionXsec::Fill_Eslice_Study(const anavar & evt, double weight, TString suffix){

  TString true_beam_PDG_str = Form("%d", abs(evt.true_beam_PDG));
  TString pitype_str = Form("%d", hadana.pitype);
  double true_KE_ff = hadana.true_ffKE;
  double KE_end = -9999.;
  double P_beam = evt.true_beam_endP * 1000.;
  double m_beam = evt.true_beam_mass;
  double E_beam = sqrt(P_beam * P_beam + m_beam * m_beam);

  for(unsigned int i = 0; i < (evt.true_beam_traj_KE)->size(); i++){
    if((evt.true_beam_traj_KE)->at(i) > 0.1){
      KE_end = (evt.true_beam_traj_KE)->at(i);
    }
    //cout << "[true_beam_traj_KE.at(" << i << ") : " << (evt.true_beam_traj_KE)->at(i) << endl;
  }

  int bin_init = true_KE_ff / 50;
  int bin_end = KE_end / 50;
  if(true_KE_ff > 10000.) return; // == Remove cases that KE_ff could be defined
  if(bin_init == bin_end){
    //cout << "true_KE_ff : " << true_KE_ff << ", KE_end : " << KE_end << endl;
    return;
  }

  // == E slice with 50 MeV binning
  double init_KE_center = 50. * (bin_init + 0.) + 25.;
  double end_KE_center = 50. * (bin_end + 0.) + 25.;
  Hist.JSFillHist(suffix, "hdaughter_KE_pion_init_" + suffix + "_" + pitype_str, init_KE_center - 50., weight, 30., 0., 1500.);
  if(pitype_str != "0") Hist.JSFillHist(suffix, "hdaughter_KE_pion_init_" + suffix + "_" + true_beam_PDG_str, init_KE_center - 50., weight, 30., 0., 1500.);
  if(bin_init - 1 == bin_end){
    //cout << Form("(true_KE_ff, KE_end) : (%f, %f), (init_KE_center - 50, end_KE_center) : (%f, %f)", true_KE_ff, KE_end, init_KE_center -50., end_KE_center) << endl;
    Hist.JSFillHist(suffix, "hdaughter_KE_pion_end_at_next_slice_" + suffix + "_" + pitype_str, end_KE_center, weight, 30., 0., 1500.);
    if(pitype_str != "0"){
      Hist.JSFillHist(suffix, "hdaughter_KE_pion_end_at_next_slice_" + suffix + "_" + true_beam_PDG_str, end_KE_center, weight, 30., 0., 1500.);
    }
  }
  Hist.JSFillHist(suffix, "hdaughter_KE_pion_end_" + suffix + "_" + pitype_str, end_KE_center, weight, 30., 0., 1500.);
  if(pitype_str != "0"){
    Hist.JSFillHist(suffix, "hdaughter_KE_pion_end_" + suffix + "_" + true_beam_PDG_str, end_KE_center, weight, 30., 0., 1500.);
    Hist.JSFillHist(suffix, "hdaughter_KE_pion_int_InElas_" + suffix + "_" + true_beam_PDG_str, end_KE_center, weight, 30., 0., 1500.);
  }
  for(int i = bin_end; i < bin_init; i++){
    double this_KE_center = 50. * (i + 0.) + 25.;
    Hist.JSFillHist(suffix, "hdaughter_KE_pion_inc_" + suffix + "_" + pitype_str, this_KE_center, weight, 30., 0., 1500.);
    if(pitype_str != "0") Hist.JSFillHist(suffix, "hdaughter_KE_pion_inc_" + suffix + "_" + true_beam_PDG_str, this_KE_center, weight, 30., 0., 1500.);
  }

  // ==== Study detailed categories of inelastic scatterings
  TVector3 unit_beam(evt.true_beam_endPx, evt.true_beam_endPy, evt.true_beam_endPz);

  int N_daughter_piplus = 0;
  int N_daughter_pizero = 0;
  int N_daughter_piminus = 0;
  int N_daughter_proton = 0;
  bool is_QE = false;
  double EQEmE = -9999.;
  for(unsigned int i = 0; i < (*evt.true_beam_daughter_ID).size(); i++){
    int this_daughter_PID = (*evt.true_beam_daughter_PDG).at(i);
    if(this_daughter_PID == 211){
      N_daughter_piplus++;
      TVector3 unit_daughter((*evt.true_beam_daughter_startPx).at(i), (*evt.true_beam_daughter_startPy).at(i), (*evt.true_beam_daughter_startPz).at(i) );
      unit_daughter = (1./ unit_daughter.Mag() ) * unit_daughter;
      double cos_theta = cos(unit_beam.Angle(unit_daughter));
      double P_daughter = (*evt.true_beam_daughter_startP).at(i) * 1000.;
      double EQE_NC_4 = Get_EQE_NC_Pion(P_daughter, cos_theta, 4., -1.);
      double this_EQEmE = EQE_NC_4 - E_beam;
      if(this_EQEmE > -50. && this_EQEmE < 50.){
        is_QE = true;
        EQEmE = this_EQEmE;
      }
    }
    else if(this_daughter_PID == -211) N_daughter_piminus++;
    else if(this_daughter_PID == 111) N_daughter_pizero++;
    else if(this_daughter_PID == 2212) N_daughter_proton++;
    else continue;
  }
  // == charge exchange
  if(N_daughter_piplus == 0 && N_daughter_pizero == 1 && N_daughter_piminus == 0){
    if(pitype_str != "0") Hist.JSFillHist(suffix, "hdaughter_KE_pion_int_ChargeEx_" + suffix + "_" + true_beam_PDG_str, end_KE_center, weight, 30., 0., 1500.);
  }
  else if(N_daughter_piplus == 0 && N_daughter_pizero == 0 && N_daughter_piminus == 1){
    if(pitype_str != "0") Hist.JSFillHist(suffix, "hdaughter_KE_pion_int_DoubleChargeEx_" + suffix + "_" + true_beam_PDG_str, end_KE_center, weight, 30., 0., 1500.);
  }
  else if(N_daughter_piplus == 1 && N_daughter_pizero == 0 && N_daughter_piminus == 0){
    if(pitype_str != "0") Hist.JSFillHist(suffix, "hdaughter_KE_pion_int_QuasiElas_" + suffix + "_" + true_beam_PDG_str, end_KE_center, weight, 30., 0., 1500.);
  }
  else if(N_daughter_piplus == 0 && N_daughter_pizero == 0 && N_daughter_piminus == 0){
    if(pitype_str != "0") Hist.JSFillHist(suffix, "hdaughter_KE_pion_int_Absorption_" + suffix + "_" + true_beam_PDG_str, end_KE_center, weight, 30., 0., 1500.);
  }
  else if(N_daughter_piplus + N_daughter_pizero + N_daughter_piminus > 1){
    if(pitype_str != "0") Hist.JSFillHist(suffix, "hdaughter_KE_pion_int_PiProd_" + suffix + "_" + true_beam_PDG_str, end_KE_center, weight, 30., 0., 1500.);
  }
  else{
    cout << "[PionXsec::Fill_Eslice_Study] No pion interaction Category" << endl;
  }

  if(is_QE){
    if(pitype_str != "0") Hist.JSFillHist(suffix, "hdaughter_KE_pion_int_EQE_pass_" + suffix + "_" + true_beam_PDG_str, end_KE_center, weight, 30., 0., 1500.);
    cout << "[PionXsec::Fill_Eslice_Study] EQE_pass" << endl;
  }
}

void PionXsec::FillHistQE_MCstudy(const anavar & evt, double weight, TString suffix){

  TString pitype_str = Form("%d", hadana.pitype);
  // ==== Beam related variables
  TVector3 reco_unit_beam(evt.reco_beam_allTrack_trackDirX, evt.reco_beam_allTrack_trackDirY, evt.reco_beam_allTrack_trackDirZ);
  reco_unit_beam = (1. / reco_unit_beam.Mag() ) * reco_unit_beam;
  double true_beam_End_E = evt.reco_beam_true_byHits_endE;


  int N_daughter_piplus = 0;
  int N_daughter_pizero = 0;
  int N_daughter_piminus = 0;
  int N_daughter_proton = 0;
  int N_daughter_neutron = 0;
  bool is_QE = false;
  double EQEmE = -9999.;
  for(unsigned int i = 0; i < (*evt.true_beam_daughter_ID).size(); i++){
    int this_daughter_PID = (*evt.true_beam_daughter_PDG).at(i);
    if(this_daughter_PID == 211) N_daughter_piplus++;
    else if(this_daughter_PID == -211) N_daughter_piminus++;
    else if(this_daughter_PID == 111) N_daughter_pizero++;
    else if(this_daughter_PID == 2212) N_daughter_proton++;
    else if(this_daughter_PID == 2112) N_daughter_neutron++;
    else continue;
  }
  TString N_pions_str = Form("Npiplus%d_Npizero%d_Npiminus%d_Nproton%d_Nneutron%d", N_daughter_piplus, N_daughter_pizero, N_daughter_piminus, N_daughter_proton, N_daughter_neutron);

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
    //double KE_daughter = hadana.map_BB[abs(this_PdgID)] -> MomentumtoKE(P_daughter);
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

    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_EQEmE_" + suffix + "_" + pitype_str, EQE - E_beam, weight, 6000., -3000., 3000.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_EQEmE_NC_40_" + suffix + "_" + pitype_str, EQE_NC_40 - E_beam, weight, 6000., -3000., 3000.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_EQEmE_NC_4_" + suffix + "_" + pitype_str, EQE_NC_4 - E_beam, weight, 6000., -3000., 3000.);

    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_E_beam_vs_EQE_" + suffix + "_" + pitype_str, E_beam, EQE, weight, 240., 100., 1300., 240., 100., 1300.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_E_beam_vs_EQE_NC_40_" + suffix + "_" + pitype_str, E_beam, EQE_NC_40, weight, 240., 100., 1300., 240., 100., 1300.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_E_beam_vs_EQE_NC_4_" + suffix + "_" + pitype_str, E_beam, EQE_NC_4, weight, 240., 100., 1300., 240., 100., 1300.);

    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_P_vs_EQE_" + suffix, P_daughter, EQE, weight, 1000., 0., 1000., 600., 0., 3000.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_P_vs_EQE_NC_40_" + suffix, P_daughter, EQE_NC_40, weight, 1000., 0., 1000., 600., 0., 3000.);
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_P_vs_EQE_NC_4_" + suffix, P_daughter, EQE_NC_4, weight, 1000., 0., 1000., 600., 0., 3000.);

    // == Fill events with exatly one + charged pion at final state
    if(this_PdgID != 211) continue;
    Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_EQEmE_NC_4_" + N_pions_str + "_" + suffix + "_" + pitype_str, EQE_NC_4 - E_beam, weight, 6000., -3000., 3000.);
    if(N_daughter_piplus == 1 && N_daughter_pizero == 0 && N_daughter_piminus == 0){
      Hist.JSFillHist(suffix, "hdaughter_" + particle_str + "_EQEmE_NC_4_NpiplusSel_" + suffix + "_" + pitype_str, EQE_NC_4 - E_beam, weight, 6000., -3000., 3000.);
    }
  }
}

void PionXsec::FillHistQE_Reco(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix){

  TString pitype_str = Form("%d", hadana.pitype);
  // ==== Beam related variables
  double P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  double mass_beam = 139.57;

  double P_ff_reco = P_beam_inst;
  TString Is_QE_str = "";
  if(evt.MC){
    P_ff_reco = Convert_P_Spectrometer_to_P_ff(evt, P_beam_inst, suffix, "AllTrue", 0);
  }
  double KE_ff_reco = sqrt(pow(P_ff_reco, 2) + pow(mass_beam, 2)) - mass_beam;
  double KE_end_reco = hadana.map_BB[211]->KEAtLength(KE_ff_reco, evt.reco_beam_alt_len);
  double E_end_reco = KE_end_reco + mass_beam;
  int E_end_reco_50bin = E_end_reco / 50;
  TString E_end_reco_50bin_str = Form("E_end_reco_%dto%dMeV_", 50 * E_end_reco_50bin, 50 * (E_end_reco_50bin + 1));

  Fill_Eslice_Histograms("QE_test", suffix, KE_ff_reco, KE_end_reco, -1., weight, binning_100MeV, N_binning_100MeV, false);

  TVector3 reco_unit_beam(evt.reco_beam_allTrack_trackDirX, evt.reco_beam_allTrack_trackDirY, evt.reco_beam_allTrack_trackDirZ);
  reco_unit_beam = (1. / reco_unit_beam.Mag() ) * reco_unit_beam;

  vector<RecoDaughter> pions = GetPions(daughters);

  //if(pions.size() < 1) return;

  int N_pion = 0;
  int N_pion_michel = 0;
  int N_pion_gaus_fitted = 0;
  int N_pion_likelihood_fitted = 0;
  for(unsigned int i = 0; i < pions.size(); i++){
    RecoDaughter this_pion = pions.at(i);
    double this_michelscore = this_pion.allTrack_vertex_michel_score() / this_pion.allTrack_vertex_nHits();
    TVector3 unit_daughter( this_pion.allTrack_endX() - this_pion.allTrack_startX(),
			    this_pion.allTrack_endY() - this_pion.allTrack_startY(),
			    this_pion.allTrack_endZ() - this_pion.allTrack_startZ());
    unit_daughter = (1. / unit_daughter.Mag() ) * unit_daughter;
    double cos_theta = cos(reco_unit_beam.Angle(unit_daughter));

    int pion_true_PID = -1;

    double P_daughter_true = -1.;
    double KE_daughter_true = -1.;
    int KE_daughter_50bin = -1;
    TString KE_daughter_50bin_str = "";
    if(evt.MC){
      int this_PID = this_pion.PFP_true_byHits_PDG();
      // 1 : pion, 2 : proton, 3 : muon, 4 : other
      if(this_PID == 211){
	pion_true_PID = 1;
	if(Is_QE_true) Is_QE_str = "_IsQE";
      }
      else if(this_PID == 2212) pion_true_PID = 2;
      else if(abs(this_PID) == 13) pion_true_PID = 3;
      else pion_true_PID = 4;

      P_daughter_true = this_pion.PFP_true_byHits_startP() * 1000.;
      KE_daughter_true = hadana.map_BB[211] -> MomentumtoKE(P_daughter_true);
      KE_daughter_50bin =KE_daughter_true / 50;
      KE_daughter_50bin_str = Form("KE_true_%dto%dMeV_", 50 * KE_daughter_50bin, 50 * (KE_daughter_50bin + 1));
    }

    if(this_michelscore > 0.5){
      double KE_daughter = hadana.map_BB[211]->KEFromRangeSpline(this_pion.allTrack_alt_len());
      double P_daughter = hadana.map_BB[211]->KEtoMomentum(KE_daughter);
      if(evt.MC){
	double P_res = (P_daughter - P_daughter_true) / P_daughter_true;
	Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_P_res_" + suffix + "_" + pitype_str + Is_QE_str, P_daughter_true, P_res,weight,1200., 0., 1200., 1000., -5., 5.);
        Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_PID_" + suffix + "_" + pitype_str + Is_QE_str, pion_true_PID, weight, 10, -2.5, 7.5);
        //Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_PID_" + KE_daughter_50bin_str + suffix + "_" + pitype_str + Is_QE_str, pion_true_PID, weight, 10, -2.5, 7.5);

	FillHistQE_MCstudy(evt, weight, "true_stop_pion_" + suffix);
      }

      double EQE = Get_EQE(P_daughter, cos_theta);
      double EQE_NC_40 = Get_EQE_NC_Pion(P_daughter, cos_theta, 40., -1.);
      double EQE_NC_4 = Get_EQE_NC_Pion(P_daughter, cos_theta, 4., -1.);
      double mEQE = EQE - E_end_reco;
      double mEQE_NC_40 = EQE_NC_40 - E_end_reco;
      double mEQE_NC_4 = EQE_NC_4 - E_end_reco;
      
      Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_EQEmE_NC_4_vs_daughter_StartZ_" + suffix + "_" + pitype_str, mEQE_NC_4, this_pion.allTrack_startZ(), weight, 1500., -1000., 500., 700., 0., 700);
      /*
      if(mEQE_NC_4 > 0. && mEQE_NC_4 < 200.){
	cout << "[PionXsec::FillHistQE_Reco] (run:event)\t(" << evt.run << ":" << evt.event << "), mEQE_NC_4 : " << mEQE_NC_4 << endl;
      }
      else continue;
      */
      // == Basic kinematics for reco pion
      Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_KE_" + suffix + "_" + pitype_str + Is_QE_str, KE_daughter, weight, 3000., 0., 3000.);
      Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_cos_" + suffix + "_" + pitype_str + Is_QE_str, cos_theta, weight, 200., -1., 1.);
      Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_length_" + suffix + "_" + pitype_str + Is_QE_str, this_pion.allTrack_alt_len(), weight, 700., 0., 700.);

      // == QE variables
      Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_EQEmE_NC_4_" + suffix + "_" + pitype_str + Is_QE_str, mEQE_NC_4, weight, 6000., -3000., 3000.);
      Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_EQEmE_NC_4_vs_KE_daughter_" + suffix + "_" + pitype_str + Is_QE_str, mEQE_NC_4, KE_daughter, weight, 1500., -1000., 500., 1000., 0., 1000.);
      
      Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_EQEmE_NC_4_" + suffix + "_" + E_end_reco_50bin_str + pitype_str + Is_QE_str, mEQE_NC_4, weight, 6000., -3000., 3000.);

      if(KE_daughter > 400.){
	Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_KE_KE400_" + suffix + "_" + pitype_str + Is_QE_str, KE_daughter, weight, 3000., 0., 3000.);
	Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_cos_KE400_" + suffix + "_" + pitype_str + Is_QE_str, cos_theta, weight, 200., -1., 1.);
	Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_length_KE400_" + suffix + "_" + pitype_str + Is_QE_str, this_pion.allTrack_alt_len(), weight, 700., 0., 700.);

	Hist.JSFillHist(suffix, "hdaughter_reco_stop_pion_EQEmE_NC_4_KE400_" + suffix + "_" + pitype_str + Is_QE_str, mEQE_NC_4, weight, 6000., -3000., 3000.);
	if(evt.MC){
	  FillHistQE_MCstudy(evt, weight, "true_stop_pion_KE400_" + suffix);
	  if(mEQE_NC_4 > -100. && mEQE_NC_4 < 100.) FillHistQE_MCstudy(evt, weight, "true_stop_pion_mEQEcut_" + suffix);
	}
      }

      N_pion_michel++;
    }
    else{ // == Pions not passing Michel score cut

      //if(!(evt.run == 18825510 && evt.event == 1056))   continue;
      double daughter_length_gaus = hadana.Fit_dEdx_Residual_Length(evt, this_pion.allTrack_calibrated_dEdX_SCE(), this_pion.allTrack_resRange_SCE(), 211, false, false);
      //cout << evt.run << ":" << evt.event << endl;
      double daughter_length_likelihood = hadana.Fit_Residual_Length_Likelihood(evt, this_pion.allTrack_calibrated_dEdX_SCE(), this_pion.allTrack_resRange_SCE(), 211, false);
      //double daughter_length_likelihood = -1.;

      if(daughter_length_gaus > 0.){
	double KE_daughter = hadana.map_BB[211]->KEFromRangeSpline(daughter_length_gaus);
	double P_daughter = hadana.map_BB[211]->KEtoMomentum(KE_daughter);

	if(evt.MC){
	  int true_PDG = this_pion.PFP_true_byHits_PDG();
	  TString true_PDG_str = Form("%d", true_PDG);
	  double P_res = (P_daughter - P_daughter_true) / P_daughter_true;
	  Hist.JSFillHist(suffix, "hdaughter_reco_gaus_fit_pion_P_res_" + suffix + "_" + true_PDG_str + Is_QE_str, P_daughter_true, P_res, weight, 1200., 0., 1200., 1000., -5., 5.);
	  Hist.JSFillHist(suffix, "hdaughter_reco_gaus_fit_pion_PID_" + suffix + "_" + pitype_str + Is_QE_str, pion_true_PID, weight, 10, -2.5, 7.5);
	  //Hist.JSFillHist(suffix, "hdaughter_reco_gaus_fit_pion_PID_" + KE_daughter_50bin_str + suffix + "_" + pitype_str + Is_QE_str, pion_true_PID, weight, 10, -2.5, 7.5);
	  
          FillHistQE_MCstudy(evt, weight, "true_gaus_fit_pion_" + suffix);
	}

	double EQE = Get_EQE(P_daughter, cos_theta);
	double EQE_NC_40 = Get_EQE_NC_Pion(P_daughter, cos_theta, 40., -1.);
	double EQE_NC_4 = Get_EQE_NC_Pion(P_daughter, cos_theta, 4., -1.);
	double mEQE = EQE - E_end_reco;
	double mEQE_NC_40 = EQE_NC_40 - E_end_reco;
	double mEQE_NC_4 = EQE_NC_4 - E_end_reco;
      
	// == Basic kinematics for reco pion
	Hist.JSFillHist(suffix, "hdaughter_reco_gaus_fit_pion_KE_" + suffix + "_" + pitype_str + Is_QE_str, KE_daughter, weight, 3000., 0., 3000.);
	Hist.JSFillHist(suffix, "hdaughter_reco_gaus_fit_pion_cos_" + suffix + "_" + pitype_str + Is_QE_str, cos_theta, weight, 200., -1., 1.);
	Hist.JSFillHist(suffix, "hdaughter_reco_gaus_fit_pion_length_" + suffix + "_" + pitype_str + Is_QE_str, this_pion.allTrack_alt_len(), weight, 700., 0., 700.);
	Hist.JSFillHist(suffix, "hdaughter_reco_gaus_fit_pion_EQEmE_NC_4_vs_daughter_StartZ_" + suffix + "_" + pitype_str + Is_QE_str, mEQE_NC_4, this_pion.allTrack_startZ(), weight, 1500., -1000., 500., 700., 0., 700);

	// == QE variables
	Hist.JSFillHist(suffix, "hdaughter_reco_gaus_fit_pion_EQEmE_NC_4_" + suffix + "_" + pitype_str + Is_QE_str, mEQE_NC_4, weight, 6000., -3000., 3000.);
	Hist.JSFillHist(suffix, "hdaughter_reco_gaus_fit_pion_EQEmE_NC_4_vs_KE_daughter_" + suffix + "_" + pitype_str + Is_QE_str, mEQE_NC_4, KE_daughter, weight, 1500., -1000., 500., 1000., 0., 1000.);
	Hist.JSFillHist(suffix, "hdaughter_reco_gaus_fit_pion_EQEmE_NC_4_" + suffix + "_" + E_end_reco_50bin_str + pitype_str + Is_QE_str, mEQE_NC_4, weight, 6000., -3000., 3000.);

	if(evt.MC){
          if(mEQE_NC_4 > -100. && mEQE_NC_4 < 100.) FillHistQE_MCstudy(evt, weight, "true_gaus_fit_pion_mEQEcut_" + suffix);
        }
	N_pion_gaus_fitted++;
      }
      if(daughter_length_likelihood > 0.){
	double KE_daughter = hadana.map_BB[211]->KEFromRangeSpline(daughter_length_likelihood);
        double P_daughter = hadana.map_BB[211]->KEtoMomentum(KE_daughter);

        if(evt.MC){
	  int true_PDG = this_pion.PFP_true_byHits_PDG();
          TString true_PDG_str = Form("%d", true_PDG);
	  double P_res = (P_daughter - P_daughter_true) / P_daughter_true;
          Hist.JSFillHist(suffix, "hdaughter_reco_likelihood_fit_pion_P_res_" + suffix + "_" + true_PDG_str + Is_QE_str, P_daughter_true, P_res, weight, 1200., 0., 1200., 1000., -5., 5.);
	  Hist.JSFillHist(suffix, "hdaughter_reco_likelihood_fit_pion_PID_" + suffix + "_" + pitype_str + Is_QE_str, pion_true_PID, weight, 10, -2.5, 7.5);
          //Hist.JSFillHist(suffix, "hdaughter_reco_likelihood_fit_pion_PID_" + KE_daughter_50bin_str + suffix + "_" + pitype_str + Is_QE_str, pion_true_PID, weight, 10, -2.5, 7.5);

          FillHistQE_MCstudy(evt, weight, "true_likelihood_fit_pion_" + suffix);
	}

        double EQE = Get_EQE(P_daughter, cos_theta);
        double EQE_NC_40 = Get_EQE_NC_Pion(P_daughter, cos_theta, 40., -1.);
        double EQE_NC_4 = Get_EQE_NC_Pion(P_daughter, cos_theta, 4., -1.);
        double mEQE = EQE - E_end_reco;
        double mEQE_NC_40 = EQE_NC_40 - E_end_reco;
        double mEQE_NC_4 = EQE_NC_4 - E_end_reco;

	// == Basic kinematics for reco pion
	Hist.JSFillHist(suffix, "hdaughter_reco_likelihood_fit_pion_KE_" + suffix + "_" + pitype_str + Is_QE_str, KE_daughter, weight, 3000., 0., 3000.);
	Hist.JSFillHist(suffix, "hdaughter_reco_likelihood_fit_pion_cos_" + suffix + "_" + pitype_str + Is_QE_str, cos_theta, weight, 200., -1., 1.);
        Hist.JSFillHist(suffix, "hdaughter_reco_likelihood_fit_pion_length_" + suffix + "_" + pitype_str + Is_QE_str, this_pion.allTrack_alt_len(), weight, 700., 0., 700.);
	Hist.JSFillHist(suffix, "hdaughter_reco_likelihood_fit_pion_EQEmE_NC_4_vs_daughter_StartZ_" + suffix + "_" + pitype_str + Is_QE_str, mEQE_NC_4, this_pion.allTrack_startZ(), weight, 1500., -1000., 500., 700., 0., 700);

        // == QE variables
	Hist.JSFillHist(suffix, "hdaughter_reco_likelihood_fit_pion_EQEmE_NC_4_" + suffix + "_" + pitype_str + Is_QE_str, mEQE_NC_4, weight, 6000., -3000., 3000.);
        Hist.JSFillHist(suffix, "hdaughter_reco_likelihood_fit_pion_EQEmE_NC_4_vs_KE_daughter_" + suffix + "_" + pitype_str + Is_QE_str, mEQE_NC_4, KE_daughter, weight, 1500., -1000., 500., 1000., 0., 1000.);
	Hist.JSFillHist(suffix, "hdaughter_reco_likelihood_fit_pion_EQEmE_NC_4_" + suffix + "_" + E_end_reco_50bin_str + pitype_str + Is_QE_str, mEQE_NC_4, weight, 6000., -3000., 3000.);

	if(evt.MC){
          if(mEQE_NC_4 > -100. && mEQE_NC_4 < 100.) FillHistQE_MCstudy(evt, weight, "true_likelihood_fit_pion_mEQEcut_" + suffix);
        }
        N_pion_likelihood_fitted++;
      }
  
    }
  }

  Hist.JSFillHist(suffix, "hdaughter_reco_QE_N_pion_michel_" + suffix + "_" + pitype_str + Is_QE_str, N_pion_michel, weight, 10., -0.5, 9.5);
  Hist.JSFillHist(suffix, "hdaughter_reco_QE_N_pion_gaus_fitted_" + suffix + "_" + pitype_str + Is_QE_str, N_pion_gaus_fitted, weight, 10., -0.5, 9.5);
  Hist.JSFillHist(suffix, "hdaughter_reco_QE_N_pion_likelihood_fitted_" + suffix + "_" + pitype_str + Is_QE_str, N_pion_likelihood_fitted, weight, 10., -0.5, 9.5);
}

void PionXsec::Run_Beam(const anavar & evt, double weight, TString suffix){

  if(!Pass_Beam_PID(evt, 211)) return;

  TString beam_particle_str = "pion";

  //// == Cuts for beam instrument variables
  // == Beam scraper cut
  if(!hadana.PassBeamScraperCut(evt)) return;
  //FillHistBeam(evt, weight, beam_particle_str + "_BeamScraper_" + suffix);

  // == Beam window cut
  double P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  if(P_beam_inst < P_beam_inst_cut_lower || P_beam_inst > P_beam_inst_cut_upper) return;
  FillHistBeam(evt, weight, beam_particle_str + "_BeamWindow_" + suffix);

  //// == Cuts for TPC variables
  if(!hadana.PassPandoraSliceCut(evt)) return;
  if(!hadana.PassCaloSizeCut(evt)) return;

  // == APA3 cut                             
  if(!hadana.PassAPA3Cut(evt)) return;
  FillHistBeam(evt, weight, beam_particle_str + "_APA3_" + suffix);

  // == Michel score cut                 
  if(!hadana.PassMichelScoreCut()) return;
  FillHistBeam(evt, weight, beam_particle_str + "_MichelScore_" + suffix);

  // == Track length cut
  if(evt.reco_beam_alt_len < 10.) return;
  FillHistBeam(evt, weight, beam_particle_str + "_TrkLength_" + suffix);

  // == Track angle cut
  TVector3 unit_TPC_beam(evt.reco_beam_allTrack_trackDirX, evt.reco_beam_allTrack_trackDirY, evt.reco_beam_allTrack_trackDirZ);
  unit_TPC_beam = (1. / unit_TPC_beam.Mag() ) * unit_TPC_beam;
  TVector3 unit_spec_beam(evt.beam_inst_dirX, evt.beam_inst_dirY, evt.beam_inst_dirZ);
  unit_spec_beam = (1. / unit_spec_beam.Mag() ) * unit_spec_beam;
  double cos_theta_spec_TPC = unit_TPC_beam.Dot(unit_spec_beam);
  //if(cos_theta_spec_TPC < 0.95) return;
  if(hadana.beam_costh < 0.95) return;
  FillHistBeam(evt, weight, beam_particle_str + "_TrkAngle_" + suffix);

  // == Track Start Z cut
  if(evt.reco_beam_calo_startZ < -5.) return;
  FillHistBeam(evt, weight, beam_particle_str + "_StartZ_" + suffix);

  // == Proton veto cut
  if(!hadana.PassProtonCut()) return;
  FillHistBeam(evt, weight, beam_particle_str + "_ProtonVeto_" + suffix);

  // == Muon veto cut
  double chi2_muon = hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 13);
  if(chi2_muon < 8) return;
  FillHistBeam(evt, weight, beam_particle_str + "_MuonVeto_" + suffix);

  // == Elastic pion veto cut
  double chi2_pion = hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 211);
  if(chi2_pion < 13) return;
  FillHistBeam(evt, weight, beam_particle_str + "_ElasPiVeto_" + suffix);

  // == End of selection
  return;
}

void PionXsec::Run_Daughter(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix){

  if(!Pass_Beam_PID(evt, 211)) return;

  TString beam_particle_str = "pion";
  if(evt.MC){
    FillHistQE_MCstudy(evt, weight, "nocut_" + suffix);
  }

  //// == Cuts for beam instrument variables
  // == Beam scraper cut
  if(!hadana.PassBeamScraperCut(evt)) return;

  // == Beam window cut
  double P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  if(P_beam_inst < P_beam_inst_cut_lower || P_beam_inst > P_beam_inst_cut_upper) return;

  //// == Cuts for TPC variables
  if(!hadana.PassPandoraSliceCut(evt)) return;
  if(!hadana.PassCaloSizeCut(evt)) return;

  // == APA3 cut
  if(!hadana.PassAPA3Cut(evt)) return;

  // == Michel score cut
  if(!hadana.PassMichelScoreCut()) return;

  // == Track length cut
  if(evt.reco_beam_alt_len < 10.) return;

  // == Track angle cut
  if(hadana.beam_costh < 0.95) return;

  // == Track Start Z cut
  if(evt.reco_beam_calo_startZ < -5.) return;

  // == Proton veto cut
  if(!hadana.PassProtonCut()) return;

  // == Muon veto cut
  double chi2_muon = hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 13);
  if(chi2_muon < 8) return;

  // == Elastic pion veto cut
  double chi2_pion = hadana.Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 211);
  if(chi2_pion < 13) return;

  FillHistQE_Reco(daughters, evt, weight, suffix);

  if(evt.MC){
    FillHistDaughterTrue(daughters, evt, weight, suffix);
    FillHistDaughterPurity(daughters, evt, weight, suffix);
    FillHistQE_MCstudy(evt, weight, suffix);
  }

  // == End of selection
  return;
}

void PionXsec::SaveHistograms(){
  //outputFile->cd();
  //outputFile->Write();
  Hist.WriteHist();
}



void PionXsec::Run(anavar & evt, Long64_t nentries=-1){

  BookHistograms();

  // == Beam Momentum Setting
  TString beam_P_str = "1.0";
  if(beam_P_str == "2.0"){
    P_beam_inst_cut_upper = 2400.;
    P_beam_inst_cut_lower = 1600.;
  }
  else if(beam_P_str == "1.0"){
    P_beam_inst_cut_upper = 1200.;
    P_beam_inst_cut_lower = 800.;
  }
  else if(beam_P_str == "0.5"){
    P_beam_inst_cut_upper = 600.;
    P_beam_inst_cut_lower = 400.;
  }
  else{
    cout << "[PionXsec::Run] Not Valid Beam Momentum" << endl;
    return;
  }
  cout << "[PionXsec::Run] Beam Momentum is " << beam_P_str << " GeV" << endl;

  
  // == Run
  if (nentries == -1) nentries = evt.fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%1000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    
    //if(jentry < 646324) continue;
    //cout << jentry<<"/"<<nentries<<std::endl;

    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;

    // if (Cut(ientry) < 0) continue;
    //std::cout<<evt.run<<" "<<evt.event<<" "<<evt.MC<<" "<<evt.reco_beam_true_byE_matched<<" "<<evt.true_beam_PDG<<" "<<(*evt.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    Hist.FillHist("Cutflow", 0.5, 1., 10., 0., 10.);
    ProcessEvent(evt);
    //Fill_Eslice_Study(evt, 1., "nocut_noweight");

    if(evt.MC){
      Is_QE_true = Is_EQE(evt, 200.);
    }

    double weight = 1.;
    double weight_piOnly = 1.;
    if(evt.MC){
      // == Beam P scale
      if(P_beam_str == "2.0") P_beam_inst_scale = 2.0;
      if(P_beam_str == "1.0") P_beam_inst_scale = 1.0;
      if(P_beam_str == "0.5") P_beam_inst_scale = 0.5;
    }

    //Run_Beam(evt, 1., "noweight");
    double P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
    double momentum_reweight = 1.;
    if(evt.MC) momentum_reweight = hadana.Beam_Mom_Reweight("TrkLength", P_beam_inst);
    //cout << "[PionXsec::Run] P_beam_inst : " << P_beam_inst << ", momentum_reweight : " << momentum_reweight << endl;
    //Run_Beam(evt, momentum_reweight, "P_reweight_TrkLength");

    vector<RecoDaughter> RecoDaughters_all = GetAllRecoDaughters(evt);

    if(RecoDaughters_all.size() == 0) continue;
    Run_Daughter(RecoDaughters_all, evt, momentum_reweight, "P_reweight_TrkLength");

    RecoDaughters_all.clear();
  }
  SaveHistograms();
  //Make_dEdx_Range_Profile();
}

