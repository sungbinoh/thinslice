#include "anavar.h"
#include "HadAna.h"
#include "TMath.h"
#include "TH1D.h"
#include "TF1.h"
#include "TVector3.h"
#include <iostream>

HadAna::HadAna(){
  if (fProtonCSDACheck) {
    TFile *file_mom2csda = TFile::Open("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/archive/sam_managed_users/tjyang/data/3/a/7/a/51075a8b-cc10-4d2d-b4c6-8e4c4d2817d7-proton_mom_csda_converter.root");
    csda_range_vs_mom_sm = (TGraph *)file_mom2csda->Get("csda_range_vs_mom_sm");
  }
}

void HadAna::InitPi(){

  AddTruePDG(-13);
  AddTruePDG(13);
  AddTruePDG(211);

  SetPandoraSlicePDG(13);

  SetBeamQualityCuts();
}

void HadAna::InitP(){

  AddTruePDG(2212);

  SetPandoraSlicePDG(13);

  SetBeamQualityCuts();
}


void HadAna::AddTruePDG(int pdg){
  truepdglist.push_back(pdg);
};

bool HadAna::isSelectedPart(const anavar& evt) const{
  if (evt.reco_reconstructable_beam_event == 0) return false; // remove empty events first
  if (evt.MC){
    for (size_t i = 0; i<truepdglist.size(); ++i){
      if (evt.true_beam_PDG == truepdglist[i]) return true; // truth matched
    }
    return false;
  }
  else{ // real data
    if (evt.beam_inst_trigger == 8) return false; // is cosmics
    if (evt.beam_inst_nMomenta != 1 || evt.beam_inst_nTracks != 1) return false;
    for (size_t i = 0; i<truepdglist.size(); ++i){
      for (size_t j = 0; j<evt.beam_inst_PDG_candidates->size(); ++j){
        if ((*evt.beam_inst_PDG_candidates)[j] == truepdglist[i]) return true; // what's this used for?
      }
    }
    return false;
  }
}

bool HadAna::isCosmics(const anavar& evt) const{
  if (evt.MC) return false;
  else{
    if (evt.beam_inst_trigger == 8) return true;
    else return false;
  }
  return false;
}

void HadAna::SetPandoraSlicePDG(int pdg){
  pandora_slice_pdg = pdg;
};

void HadAna::SetBeamQualityCuts(double dx_min, double dx_max,
                                double dy_min, double dy_max,
                                double dz_min, double dz_max,
                                double dxy_min, double dxy_max,
                                double costh_min, double costh_max){
  beamcut_dx_min = dx_min; beamcut_dx_max = dx_max;
  beamcut_dy_min = dy_min; beamcut_dy_max = dy_max;
  beamcut_dz_min = dz_min; beamcut_dz_max = dz_max;
  beamcut_dxy_min = dxy_min; beamcut_dxy_max = dxy_max;
  beamcut_costh_min = costh_min; beamcut_costh_max = costh_max;
}

int HadAna::GetPiParType(const anavar& evt){

  if (!evt.MC){
    return pi::kData;
  }
  else if (evt.event%2){ // divide half of MC as fake data
    return pi::kData;
  }
  else if (!evt.reco_beam_true_byE_matched){ // the true beam track is not selected
    if (evt.reco_beam_true_byE_origin == 2) {
      return pi::kMIDcosmic;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 211){ // the selected track is a pion (but not true beam pion, so it is a secondary pion)
      return pi::kMIDpi;
    }
    else if (evt.reco_beam_true_byE_PDG == 2212){
      return pi::kMIDp;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 13){
      return pi::kMIDmu;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 11 ||
             evt.reco_beam_true_byE_PDG == 22){
      return pi::kMIDeg;
    }
    else {
      //cout<<evt.reco_beam_true_byE_PDG<<endl;
      return pi::kMIDother;
    }
  }
  else if (evt.true_beam_PDG == -13){
    return pi::kMuon;
  }
  else if (evt.true_beam_PDG == 211){
    if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
      return pi::kPiInel;
    }
    else return pi::kPiElas;
  }
  
  return pi::kMIDother;
}

int HadAna::GetPParType(const anavar& evt){

  if (!evt.MC){
    return p::kData;
  }
  else if (evt.event%2){
    return p::kData;
  }
  else if (!evt.reco_beam_true_byE_matched){
    if (evt.reco_beam_true_byE_origin == 2) {
      return p::kMIDcosmic;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 211){
      return p::kMIDpi;
    }
    else if (evt.reco_beam_true_byE_PDG == 2212){
      return p::kMIDp;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 13){
      return p::kMIDmu;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 11 ||
             evt.reco_beam_true_byE_PDG == 22){
      return p::kMIDeg;
    }
    else {
      //std::cout<<reco_beam_true_byE_PDG<<std::endl;
      return p::kMIDother;
    }
  }
  else if (evt.true_beam_PDG == 2212){
    if ((*evt.true_beam_endProcess) == "protonInelastic"){
      return p::kPInel;
    }
    else return p::kPElas;
  }
  
  return p::kMIDother;
}

bool HadAna::PassPandoraSliceCut(const anavar& evt) const{ // whether recognized by Pandora correctly

  if (fAllTrackCheck) return evt.reco_beam_type > 0;
  else return (evt.reco_beam_type == pandora_slice_pdg);
}

bool HadAna::PassBeamQualityCut(bool has_angle_cut) const{ // cut on beam entrance location and beam angle

  if (beamcut_dx_min<beamcut_dx_max){
    if (beam_dx<beamcut_dx_min)
      return false;
    if (beam_dx>beamcut_dx_max)
      return false;
  }

  if (beamcut_dy_min<beamcut_dy_max){
    if (beam_dy<beamcut_dy_min)
      return false;
    if (beam_dy>beamcut_dy_max)
      return false;
  }

  if (beamcut_dz_min<beamcut_dz_max){
    if (beam_dz<beamcut_dz_min)
      return false;
    if (beam_dz>beamcut_dz_max)
      return false;
  }

  if (beamcut_dxy_min<beamcut_dxy_max){
    if (beam_dxy<beamcut_dxy_min)
      return false;
    if (beam_dxy>beamcut_dxy_max)
      return false;
  }

  if (has_angle_cut && beamcut_costh_min<beamcut_costh_max){
    if (beam_costh<beamcut_costh_min)
      return false;
    if (beam_costh>beamcut_costh_max)
      return false;
  }

  return true;
}

bool HadAna::PassAPA3Cut(const anavar& evt) const{ // only use track in the first TPC
  return true;
  double cutAPA3_Z = 220.;
  
  if (fAllTrackCheck) return evt.reco_beam_calo_endZ_allTrack < cutAPA3_Z;
  else return evt.reco_beam_calo_endZ < cutAPA3_Z;
}

bool HadAna::PassCaloSizeCut(const anavar& evt) const{ // Require hits information in collection plane
  
  if (fAllTrackCheck) return !(evt.reco_beam_calo_wire_allTrack->empty());
  else return !(evt.reco_beam_calo_wire->empty());
}

bool HadAna::PassMichelScoreCut() const{ // further veto muon tracks according to Michel score
  
  return daughter_michel_score < 0.55;
}

bool HadAna::PassProtonCut() const{ // to remove proton background

  return chi2_proton > 80; //median_dEdx < 2.4;
}

bool HadAna::PassPiCuts(const anavar& evt) const{
  return PassPandoraSliceCut(evt)&&
    PassCaloSizeCut(evt)&&
    PassBeamQualityCut()&&
    PassAPA3Cut(evt)&&
    PassMichelScoreCut()&&
    PassProtonCut();
}

bool HadAna::PassPCuts(const anavar& evt) const{
  return PassPandoraSliceCut(evt)&&
    PassCaloSizeCut(evt)&&
    PassBeamQualityCut();
}

// == Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf
const double rho = 1.39; // [g/cm3], density of LAr
const double K = 0.307075; // [MeV cm2 / mol]
const double A = 39.948; // [g / mol], atomic mass of Ar
const double I = 188.0e-6; // [MeV], mean excitation energy
const double Me = 0.511; // [Mev], mass of electron 
const double density_C = 5.2146;
const double density_y0 = 0.2;
const double density_y1 = 3.0;
const double density_a = 0.19559;
const double density_k = 3.0;
double HadAna::Density_Correction(double beta, double gamma){
  // == Estimate the density correction
  double density_y = TMath::Log10(beta * gamma);
  double ln10 = TMath::Log(10);
  double this_delta = 0.;
  if(density_y > density_y1){
    this_delta = 2.0 * ln10 * density_y - density_C;
  }
  else if (density_y < density_y0){
    this_delta = 0.;
  }
  else{
    this_delta = 2.0 * ln10 * density_y - density_C + density_a * pow(density_y1 - density_y, density_k);
  }

  return this_delta;
}
double HadAna::dEdx_Bethe_Bloch(double KE, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));
  double delta = Density_Correction(beta, gamma);

  // == dE/dx with the density correction
  double f = rho * K * (18.0 / A) * pow(1. / beta, 2);
  double a0 = 0.5 * TMath::Log(2.0 * Me * pow(beta * gamma, 2) * Wmax / (I * I));
  double this_dEdx = f * ( a0 - pow(beta, 2) - delta / 2.0); // [MeV/cm]

  return this_dEdx;
}
double HadAna::Get_Landau_xi(double KE, double dx, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double xi = rho * dx * 0.5 * K * (18.0 / A) * pow(1. / beta, 2);
  return xi;
}
double HadAna::dpdx_Bethe_Bloch(double KE, double dx, double mass){
  // == Energy loss by a single particle should be described by the Landau-Vavilov distribution
  // == It provides MPV of dE/dx for a pitch
  // == https://pdg.lbl.gov/2021/reviews/rpp2020-rev-passage-particles-matter.pdf
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double delta = Density_Correction(beta, gamma);
  double xi = rho * dx * 0.5 * K * (18.0 / A) * pow(1. / beta, 2);
  double a0 = 2.0 * Me * pow(beta * gamma, 2) / I;
  double this_dpdx = (xi / dx) * (TMath::Log(a0) + + TMath::Log(xi / I) + 0.2 - pow(beta, 2) - Density_Correction(beta, gamma)) ;

  return this_dpdx;
}

void HadAna::Draw_dpdx_vs_KE(double dx, double mass, TString suffix){

  double KE_min = 100; 
  double KE_max = 10000; // 10 GeV
  double KE_step = 10.0;
  int N_Step = (KE_max - KE_min) / KE_step + 1;
  vector<double> dpdx_collection;
  vector<double> KE_collection;
  for(int i = 0; i < N_Step; i++){
    double this_KE = KE_min + KE_step * (i + 0.);
    double this_dpdx = dpdx_Bethe_Bloch(this_KE, dx, mass);
    KE_collection.push_back(this_KE);
    dpdx_collection.push_back(this_dpdx);
  }
  TGraph *dpdx_vs_KE_gr = new TGraph(N_Step, &KE_collection[0], &dpdx_collection[0]);
  dpdx_vs_KE_gr -> SetName("dpdx_vs_KE_gr_" + suffix);
  dpdx_vs_KE_gr -> Write();
  dpdx_collection.clear();
  KE_collection.clear();

}

double HadAna::Get_Landau_P(double MPV, double FWHM, double x){

  TF1 *this_landau = new TF1("landauDistr", "TMath::Landau(x,[0],[1],1)", 0, 10);
  this_landau -> SetParameter(0, MPV);
  this_landau -> SetParameter(1, FWHM);
  int N_step = 1000;
  double x_step = 10. / N_step;
  double area = 0.;
  for(int i = 0; i < N_step; i++){
    double this_x = 0. + (i + 0.) * x_step;
    double this_y = this_landau -> Eval(this_x);
    area = area + this_y * x_step;
    if(this_x > x) break;
  }
  
  delete this_landau;

  return area;

}

double HadAna::Fit_Beam_Hit_dEdx_Bethe_Bloch(const anavar& evt, int PID){

  double this_KE = 0.; // == default KE value : 0 MeV
  int N_max = 100; // == Maximum number of hits used for the Bethe-Bloch fitting
  
  // == PID input : mass hypothesis
  double this_mass = 105.658; // == default : muon mass, https://pdg.lbl.gov/2022/listings/rpp2022-list-muon.pdf

  if(PID == 2212) this_mass = 938.272; // == https://pdg.lbl.gov/2022/listings/rpp2022-list-p.pdf
  else if(abs(PID) == 211) this_mass = 139.57; // == https://pdg.lbl.gov/2022/listings/rpp2022-list-pi-plus-minus.pdf
  else best_fit_KE = this_KE;

  // == Scan KE and fit
  double best_KE = 100.0; // [MeV]
  double best_likelihood = 99999.;
  double initial_KE = 100.0; // [MeV]
  double final_KE = 10000.0; // [MeV], 10 GeV
  double KE_step = 10.0; // [MeV]
  int N_KE_trial = (final_KE - initial_KE) / KE_step;
  int this_N_calo = evt.reco_beam_calo_Z->size();
  int this_N_hits = TMath::Min(this_N_calo, N_max);
  vector<double> likelihood_vector;
  vector<double> KE_vector;

  for(int i = 0; i < N_KE_trial; i++){
    label_i_KE : 
    double this_KE = initial_KE + KE_step * (i + 0.);
    double dEdx_measured = 0.;
    double dE_measured = 0.;
    double this_KE_likelihood = 0.;
    for(int j = 0; j < this_N_hits; j++){
      this_KE = this_KE - dE_measured;
      if(this_KE < 0.){ // == leave this int j for loop
	i = i + 1;
	j = 0;
	goto label_i_KE;
      }
      double this_pitch = (*evt.reco_beam_TrkPitch_SCE)[j];
      dEdx_measured =  (*evt.reco_beam_calibrated_dEdX_SCE)[j];
      double this_dEdx_theory = dpdx_Bethe_Bloch(this_KE, this_pitch, this_mass);
      double FWHM = 4 * Get_Landau_xi(this_KE, this_pitch, this_mass);
      double P_MPV = Get_Landau_P(this_dEdx_theory, FWHM, this_dEdx_theory);
      double P_dEdx_measured = Get_Landau_P(this_dEdx_theory, FWHM, dEdx_measured);
      if(dEdx_measured > this_dEdx_theory){
	P_dEdx_measured = 1 - P_dEdx_measured;
	P_MPV = 1 - P_MPV;
      }
      double this_likelihood = TMath::Log(P_dEdx_measured / P_MPV);
      this_KE_likelihood = this_KE_likelihood + this_likelihood;
      //cout << "SB debug, dEdx_measured : " << dEdx_measured << ", this_dEdx_theory : " << this_dEdx_theory << ", this_KE : " << this_KE << " / " << initial_KE + KE_step * (i + 0.) << ", this_pitch : " << this_pitch << endl;
      //cout << "SB debug, this_likelihood : " << this_likelihood << endl;
      dE_measured = dEdx_measured * this_pitch;
    }

    this_KE_likelihood = -2.0 * this_KE_likelihood;
    //cout << "SB debug, " << initial_KE + KE_step * (i + 0.) << ", this_KE_likelihood : " << this_KE_likelihood << endl;
    if(this_KE_likelihood < best_likelihood){
      best_likelihood = this_KE_likelihood;
      best_KE = initial_KE + KE_step * (i + 0.);
    }
    likelihood_vector.push_back(this_KE_likelihood);
    KE_vector.push_back(initial_KE + KE_step * (i + 0.));
  }

  for(unsigned int i = 0; i < likelihood_vector.size(); i++){
    likelihood_vector[i] = likelihood_vector[i] - best_likelihood;
  }

  double best_mom = sqrt(pow(best_KE, 2) + 2.0 * best_KE * this_mass);
  cout << "SB debug, [" << evt.run << ":" << evt.event << "] best_KE : " << best_KE << ", best_likelihood : " << best_likelihood << ", best_mom : " << best_mom << ", this_mass : " << this_mass << endl;

  TGraph *likelihood_gr = new TGraph(KE_vector.size(), &KE_vector[0], &likelihood_vector[0]);
  likelihood_gr -> SetName(Form("Run%d_Evt%d_PID%d", evt.run, evt.event, PID));
  //likelihood_gr -> Write();
  
  likelihood_vector.clear();
  KE_vector.clear();

  delete likelihood_gr;

  return best_mom;
}

void HadAna::Draw_Landau(double KE, double mass, TString suffix){

  double MPV = dpdx_Bethe_Bloch(KE, 0.65, mass);
  double FWHM = 4 * Get_Landau_xi(KE, 0.65, mass);

  TF1 *landauDistr = new TF1("landauDistr", "TMath::Landau(x,[0],[1],1)", 0, 10);
  landauDistr -> SetName("Landau_" + suffix);
  landauDistr -> SetParameter(0, MPV);
  landauDistr -> SetParameter(1, FWHM);
  landauDistr -> Write();

  int N_step = 1000;
  double x_step = 10. / N_step;
  double area = 0.;
  for(int i = 0; i < N_step; i++){
    double this_x = 0. + (i + 0.) * x_step;
    double this_y = landauDistr -> Eval(this_x);
    area = area + this_y * x_step;
    //cout << "SB debug, KE : " << KE << ", mass : " << mass << ", (x / MPV, y) : (" << this_x << " / " << MPV << ", " << this_y << "), area = " << area << endl;
  }

  delete landauDistr;
  
}

void HadAna::ProcessEvent(const anavar& evt){
  
  pitype = GetPiParType(evt);
  ptype = GetPParType(evt);

  median_dEdx = -1;
  chi2_proton = -1;
  daughter_michel_score = -999;

  if (fAllTrackCheck) {} // removed for not in use
  else {
    if (!evt.reco_beam_calo_wire->empty()){
      chi2_proton = evt.reco_beam_Chi2_proton/evt.reco_beam_Chi2_ndof;
      median_dEdx = TMath::Median(evt.reco_beam_calibrated_dEdX_SCE->size(), &(*evt.reco_beam_calibrated_dEdX_SCE)[0]);//TMath::Median(evt.reco_beam_calibrated_dEdX_SCE->size(), &(*evt.reco_beam_calibrated_dEdX_SCE)[0]);
  //    daughter_michel_score = 0;
  //    int nhits = 0;
  //    for (size_t i = 0; i<reco_daughter_PFP_michelScore_collection->size(); ++i){
  //      nhits += (*reco_daughter_PFP_nHits_collection)[i];
  //      daughter_michel_score += (*reco_daughter_PFP_michelScore_collection)[0] * (*reco_daughter_PFP_nHits_collection)[i];
  //    }
  //    if (nhits) daughter_michel_score/=nhits;
  //    else daughter_michel_score = -999;
      if (evt.reco_beam_vertex_nHits) daughter_michel_score = evt.reco_beam_vertex_michel_score_weight_by_charge;//evt.reco_beam_vertex_michel_score/evt.reco_beam_vertex_nHits;
    }

    beam_dx = -999;
    beam_dy = -999;
    beam_dz = -999;
    beam_dxy = -999;
    beam_costh = -999;

    if (!evt.reco_beam_calo_wire->empty()){

      //if(evt.reco_beam_true_byE_PDG == 211) Fit_Hit_dEdx_Bethe_Bloch(evt, 211);
      //if(evt.reco_beam_true_byE_PDG == 2212) Fit_Hit_dEdx_Bethe_Bloch(evt, 2212);

      /*TVector3 pt0(evt.reco_beam_calo_startX,
                   evt.reco_beam_calo_startY,
                   evt.reco_beam_calo_startZ);
      TVector3 pt1(evt.reco_beam_calo_endX,
                   evt.reco_beam_calo_endY,
                   evt.reco_beam_calo_endZ);
      TVector3 dir = pt1 - pt0;
      dir = dir.Unit();*/
      //forced track info
      TVector3 pt0(evt.reco_beam_calo_startX,
                   evt.reco_beam_calo_startY,
                   evt.reco_beam_calo_startZ);
      TVector3 pt1(evt.reco_beam_calo_endX,
                   evt.reco_beam_calo_endY,
                   evt.reco_beam_calo_endZ);
      TVector3 dir = pt1 - pt0;
      dir = dir.Unit();

      if (evt.MC){
        TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180),
                         cos(beam_angleY_mc*TMath::Pi()/180),
                         cos(beam_angleZ_mc*TMath::Pi()/180));
        beamdir = beamdir.Unit();
        /*beam_dx = (evt.reco_beam_calo_startX - beam_startX_mc)/beam_startX_rms_mc;
        beam_dy = (evt.reco_beam_calo_startY - beam_startY_mc)/beam_startY_rms_mc;
        beam_dz = (evt.reco_beam_calo_startZ - beam_startZ_mc)/beam_startZ_rms_mc;
        beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
        beam_costh = dir.Dot(beamdir);*/
        //forced track info
        beam_dx = (evt.reco_beam_calo_startX - beam_startX_mc)/beam_startX_rms_mc;
        beam_dy = (evt.reco_beam_calo_startY - beam_startY_mc)/beam_startY_rms_mc;
        beam_dz = (evt.reco_beam_calo_startZ - beam_startZ_mc)/beam_startZ_rms_mc;
        beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
        beam_costh = dir.Dot(beamdir);
      }
      else{
        TVector3 beamdir(cos(beam_angleX_data*TMath::Pi()/180),
                         cos(beam_angleY_data*TMath::Pi()/180),
                         cos(beam_angleZ_data*TMath::Pi()/180));
        beamdir = beamdir.Unit();
        /*beam_dx = (evt.reco_beam_calo_startX - beam_startX_data)/beam_startX_rms_data;
        beam_dy = (evt.reco_beam_calo_startY - beam_startY_data)/beam_startY_rms_data;
        beam_dz = (evt.reco_beam_calo_startZ - beam_startZ_data)/beam_startZ_rms_data;
        beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
        beam_costh = dir.Dot(beamdir);*/
        //forced track info
        beam_dx = (evt.reco_beam_calo_startX - beam_startX_data)/beam_startX_rms_data;
        beam_dy = (evt.reco_beam_calo_startY - beam_startY_data)/beam_startY_rms_data;
        beam_dz = (evt.reco_beam_calo_startZ - beam_startZ_data)/beam_startZ_rms_data;
        beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
        beam_costh = dir.Dot(beamdir);
      }
    }

    dEdx_5cm = -1;
    /*
    if (!reco_beam_calibrated_dEdX_SCE->empty()){
      dEdx_5cm = 0;
      int nhits = 0;
      for (int i = 0; i<reco_beam_calibrated_dEdX_SCE->size(); ++i){
        std::cout<<i<<" "<<reco_beam_resRange->back()-(*reco_beam_resRange)[i]<<" "<<(*reco_beam_calibrated_dEdX_SCE)[i]<<endl;
        if (std::abs(reco_beam_resRange->back()-(*reco_beam_resRange)[i])<5){
          dEdx_5cm += (*reco_beam_calibrated_dEdX_SCE)[i];
          ++nhits;
        }
      }
      if (nhits) dEdx_5cm/=nhits;
      else dEdx_5cm = -1;
    }
    */

    //if (event == 78467) cout<<reco_beam_calibrated_dEdX_SCE->size()<<endl;
    //cout<<reco_beam_calibrated_dEdX_SCE->size()<<endl;
    if (!evt.reco_beam_calibrated_dEdX_SCE->empty()){ // what's this used for?
      //dEdx_5cm = 0;
      //int nhits = 0;
      std::vector<double> vdEdx;
      for (int i = 0; i<evt.reco_beam_calibrated_dEdX_SCE->size(); ++i){
        //std::cout<<i<<" "<<reco_beam_resRange->back()-(*reco_beam_resRange)[i]<<" "<<(*reco_beam_calibrated_dEdX_SCE)[i]<<endl;
        //if (event == 78467) cout<<(*reco_beam_resRange)[i]<<" "<<(*reco_beam_calo_Z)[i]<<" "<<(*reco_beam_calibrated_dEdX_SCE)[i]<<endl;
        if ((*evt.reco_beam_resRange_SCE)[i]<5){
          vdEdx.push_back((*evt.reco_beam_calibrated_dEdX_SCE)[i]);
          //dEdx_5cm += (*reco_beam_calibrated_dEdX_SCE)[i];
          //++nhits;
        }
      }
      //if (nhits) dEdx_5cm/=nhits;
      //else dEdx_5cm = -1;
      if (!vdEdx.empty()){
        dEdx_5cm = TMath::Median(vdEdx.size(), &vdEdx[0]);
      }
    }
  //  if (!MC && reco_beam_PFP_trackScore_collection>=0 && reco_beam_PFP_trackScore_collection<0.01){
  //    cout<<run<<" "<<event<<endl;
  //  }

    // calculate true track length
    int start_idx = -1;
    true_trklen_accum.reserve(evt.true_beam_traj_Z->size()); // initialize true_trklen_accum
    for (int i=0; i<evt.true_beam_traj_Z->size(); i++){
      if ((*evt.true_beam_traj_Z)[i] >= 0){
        start_idx = i-1; // the trajectory point before entering the TPC
        if (start_idx < 0) start_idx = -1;
        break;
      }
      true_trklen_accum[i] = 0.; // initialize true_trklen_accum
    }
    true_trklen = -1999; // initialize
    if (start_idx >= 0){
      for (int i=start_idx+1; i<evt.true_beam_traj_Z->size(); i++){
        if (i == start_idx+1) {
          true_trklen = sqrt( pow( (*evt.true_beam_traj_X)[i]-(*evt.true_beam_traj_X)[i-1], 2)
                              + pow( (*evt.true_beam_traj_Y)[i]-(*evt.true_beam_traj_Y)[i-1], 2)
                              + pow( (*evt.true_beam_traj_Z)[i]-(*evt.true_beam_traj_Z)[i-1], 2)
                              ) * (*evt.true_beam_traj_Z)[i]/((*evt.true_beam_traj_Z)[i]-(*evt.true_beam_traj_Z)[i-1]);
        }
        else{
          true_trklen += sqrt( pow( (*evt.true_beam_traj_X)[i]-(*evt.true_beam_traj_X)[i-1], 2)
                              + pow( (*evt.true_beam_traj_Y)[i]-(*evt.true_beam_traj_Y)[i-1], 2)
                              + pow( (*evt.true_beam_traj_Z)[i]-(*evt.true_beam_traj_Z)[i-1], 2)
                              );
        }
        true_trklen_accum[i] = true_trklen;
      }
    }

    // calculate reco track length
    reco_trklen_accum.reserve(evt.reco_beam_calo_Z->size());
    reco_trklen = -1999;
    for (int i=1; i<evt.reco_beam_calo_Z->size(); i++){
      if (i == 1) reco_trklen = 0;
      reco_trklen += sqrt( pow( (*evt.reco_beam_calo_X)[i]-(*evt.reco_beam_calo_X)[i-1], 2)
                          + pow( (*evt.reco_beam_calo_Y)[i]-(*evt.reco_beam_calo_Y)[i-1], 2)
                          + pow( (*evt.reco_beam_calo_Z)[i]-(*evt.reco_beam_calo_Z)[i-1], 2)
                          );
      reco_trklen_accum[i] = reco_trklen;
    }
    // front-face energy
    true_ffKE = 999999.;
    if (start_idx >= 0) true_ffKE = (*evt.true_beam_traj_KE)[start_idx+1] + 2.18*(true_trklen_accum)[start_idx+1];
  }
  
  energy_calorimetry_SCE = 0; //MeV
  for (int i=0; i<evt.reco_beam_calibrated_dEdX_SCE->size(); i++){
    energy_calorimetry_SCE += (*evt.reco_beam_calibrated_dEdX_SCE)[i]*(*evt.reco_beam_TrkPitch_SCE)[i];
  }
  //cout<<evt.beam_particle_scores->size()<<"\t"<<(*evt.beam_particle_scores)[0]<<endl;
  if (evt.beam_particle_scores->size())
    beam_score = (*evt.beam_particle_scores)[0];
  else
    beam_score = -999.;
  //cout<<"$$$"<<evt.reco_beam_alt_len<<"\t"<<reco_trklen<<endl;//the two are the same
  // reco_trklen = evt.reco_beam_alt_len; // they should be the same
  if (fProtonCSDACheck)
    trklen_csda_proton = reco_trklen / csda_range_vs_mom_sm->Eval(evt.beam_inst_P);
  else
    trklen_csda_proton = -999;
}
