#include "Beam_Study.h"
#include "HadAna.h"
#include "anavar.h"
#include "TH1D.h"
#include "util.h"
#include <iostream>

using namespace std;

Beam_Study::Beam_Study(){
  hadana.Init();
}

double Beam_Study::Gaussian_Reweight(double mu1, double sigma1, double mu2, double sigma2, double x_min, double x_max, double x){
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
  //cout << "[Beam_Study::Gaussian_Reweight] true beam P : " << x << ", y1 : " << y1 << ", y2 : " << y2 << ", y1/y2 : " << y1/y2 << endl;
  if(out > 10.) return 10.;
  return y1/y2;
}

double Beam_Study::Beam_TrueP_Reweight(const anavar & evt){
  double out = 1.;
  /*
  if(evt.MC){

  }
  */
  return out;
}

void Beam_Study::BookHistograms(){
  //cout << "[Beam_Study::BookHistograms] Start" << endl;
  Hist.outfile = TFile::Open(fOutputFileName.c_str(), "recreate");
}

bool Beam_Study::Pass_Beam_PID(const anavar & evt, int PID){

  // == Study only proton and pion/muon beams
  if(PID != 2212 && PID != 211) return false;

  vector<int> PDGID_candiate_vec;
  PDGID_candiate_vec.clear();
  if(PID == 2212){
    PDGID_candiate_vec.push_back(2212);
  }

  if(PID == 211){
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
        if ((*evt.beam_inst_PDG_candidates)[j] == PDGID_candiate_vec[i]) return true; // Matching data PID result and PDGID candiates
      }
    }
    return false;
  }
}

bool Beam_Study::Pass_KE_diff_3sigma_Cut(double KE_RecoBeam, double KE_diff, TString suffix, TString key){

  bool is_scraper = false;
  if(suffix.Contains("pion") && key == "True"){
    if(KE_RecoBeam > 700.  && KE_RecoBeam < 800.  && KE_diff > 54.5 ) is_scraper = true;
    if(KE_RecoBeam > 800.  && KE_RecoBeam < 900.  && KE_diff > 71.5 ) is_scraper = true;
    if(KE_RecoBeam > 900.  && KE_RecoBeam < 1000. && KE_diff > 88.3 ) is_scraper = true;
    if(KE_RecoBeam > 1000. && KE_RecoBeam < 1100. && KE_diff > 102.3) is_scraper = true;
  }
  if(suffix.Contains("proton") && key == "True"){
    if(KE_RecoBeam > 300. && KE_RecoBeam < 400 && KE_diff > 58.8) is_scraper = true;
    if(KE_RecoBeam > 400. && KE_RecoBeam < 500 && KE_diff > 74.8) is_scraper = true;
    if(KE_RecoBeam > 500. && KE_RecoBeam < 600 && KE_diff > 92.1) is_scraper = true;
  }
  if(suffix.Contains("proton") && key == "FittedElas"){
    if(KE_RecoBeam > 300. && KE_RecoBeam < 400 && KE_diff > 61.9) is_scraper = true;
    if(KE_RecoBeam > 400. && KE_RecoBeam < 500 && KE_diff > 78.7) is_scraper = true;
    if(KE_RecoBeam > 500. && KE_RecoBeam < 600 && KE_diff > 93.4) is_scraper = true;
  }
  if(suffix.Contains("proton") && key == "FittedFakeData"){
    if(KE_RecoBeam > 300. && KE_RecoBeam < 400 && KE_diff > 67.7) is_scraper = true;
    if(KE_RecoBeam > 400. && KE_RecoBeam < 500 && KE_diff > 83.9) is_scraper = true;
    if(KE_RecoBeam > 500. && KE_RecoBeam < 600 && KE_diff > 100.9) is_scraper = true;
  }
  if(suffix.Contains("proton") && key == "FittedData"){
    // == 3 sigma

    if(KE_RecoBeam > 300. && KE_RecoBeam < 400 && KE_diff > 71.8) is_scraper = true;
    if(KE_RecoBeam > 400. && KE_RecoBeam < 500 && KE_diff > 84.9) is_scraper = true;
    if(KE_RecoBeam > 500. && KE_RecoBeam < 600 && KE_diff > 99.5) is_scraper = true;


    // == Based on distributions of elastic and misID:p
    /*
    if(KE_RecoBeam > 300. && KE_RecoBeam < 400 && KE_diff > 120.0) is_scraper = true;
    if(KE_RecoBeam > 400. && KE_RecoBeam < 500 && KE_diff > 140.0) is_scraper = true;
    if(KE_RecoBeam > 500. && KE_RecoBeam < 600 && KE_diff > 160.0) is_scraper = true;
    */
  }

  return is_scraper;
}

double Beam_Study::Convert_KE_Spectrometer_to_KE_ff(const anavar & evt, double KE_RecoBeam, TString suffix, TString key, int syst){

  double out = KE_RecoBeam;
  double delta_E = 0.;
  double p0 = 0., p1 = 0., p2 = 0.;
  if(suffix.Contains("proton")){
    if(key == "ElasTrue"){
      if(syst == 0){
	p0 = 39.78;
	p1 = -0.2396;
	p2 = 0.0004498;
      }
      else if(syst == -1){
	p0 = 17.92;
	p1 = -0.1334;
	p2 = 0.0003075;
      }
      else if(syst == 1){
	p0 = 60.52;
	p1 = -0.3404;
	p2 = 0.0005857;
      }
      else{
	return out;
      }
    }
    else if(key == "AllTrue"){
      if(syst == 0){
	p0 = 51.33;
	p1 = -0.2954;
	p2 = 0.0005164;
      }
      else if(syst == -1){
	p0 = 28.91;
	p1 = -0.1827;
	p2 = 0.0003602;
      }
      else if(syst == 1){
	p0 = 72.64;
	p1 = -0.4028;
	p2 = 0.0006661;
      }
      else{
	return out;
      }
    }
    else if(key == "ElasFitted"){
      if(syst == 0){
	p0 = 27.68;
	p1 = -0.1727;
	p2 = 0.0003779;
      }
      else if(syst == -1){
	p0 = -6.135;
	p1 = -0.01069;
	p2 = 0.0001675;
      }
      else if(syst == 1){
	p0 = 60.19;
	p1 = -0.3284;
	p2 = 0.0005808;
      }
      else{
	return out;
      }
    }
    else if(key == "AllFitted"){
      if(syst == 0){
	p0 = 37.57;
	p1 = -0.2144;
	p2 = 0.0004282;
      }
      else if(syst == -1){
	p0 = 12.66;
	p1 = -0.09451;
	p2 = 0.0002699;
      }
      else if(syst == 1){
	p0 = 61.69;
	p1 = -0.3305;
	p2 = 0.0005820;
      }
      else{
	return out;
      }
    }
    else if(key == "Data"){
      if(syst == 0){
	p0 = 47.35;
	p1 = -0.1748;
	p2 = 0.0003067;
      }
      else if(syst == -1){
	p0 = 28.27;
	p1 = -0.07865;
	p2 = 0.0001744;
      }
      else if(syst == 1){
	p0 = 65.99;
	p1 = -0.2688;
	p2 = 0.0004365;
      }
      else{
	return out;
      }
    }
    else {
      return out;
    }
  }
  else if(suffix.Contains("pion")){
    if(key == "AllTrue"){
      if(syst == 0){
	p0 = 222.3;
	p1 = -0.6348;
	p2 = 0.0004500;
      }
      else if(syst == -1){
	p0 = 178.3;
	p1 = -0.5222;
	p2 = 0.0003757;
      }
      else if(syst == 1){
	p0 = 289.8;
	p1 = -0.8049;
	p2 = 0.0005592;
      }
      else{
	return out;
      }
    }
  }
  delta_E = p0 + p1 * out + p2 * out * out;

  return out - delta_E;
}

bool Beam_Study::Pass_Beam_Scraper_Cut(double x, double y, TString suffix, TString key){

  bool out = false;

  double sigma = 4.8;
  if(key == "mc"){
    double center_x = -29.6;
    double center_y = 422;
    double this_distance = sqrt( pow(center_x - x, 2.) + pow(center_y - y, 2.) );
    if(this_distance < 1.4 * sigma) out = true;
  }
  if(key == "data"){
    double center_x = -32.16;
    double center_y = 422.7;
    double this_distance = sqrt( pow(center_x - x, 2.) + pow(center_y - y, 2.) );
    if(this_distance < 1.2 * sigma) out = true;
  }

  return out;
}

void Beam_Study::ProcessEvent(const anavar & evt){
  hadana.ProcessEvent(evt);
}

void Beam_Study::FillHistBeam(const anavar & evt, double weight, TString suffix, double KE_fit_gaussian){

  double P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  double mass_beam = 139.57;
  if(suffix.Contains("proton")) mass_beam = 938.272;
  double KE_beam_inst = sqrt(pow(P_beam_inst, 2) + pow(mass_beam, 2)) - mass_beam;

  // == Various KE_ff reco values
  double KE_ff_reco_ElasTrue = Convert_KE_Spectrometer_to_KE_ff(evt, KE_beam_inst, suffix, "ElasTrue");
  double KE_ff_reco_AllTrue = Convert_KE_Spectrometer_to_KE_ff(evt, KE_beam_inst, suffix, "AllTrue");
  double KE_ff_reco_ElasFitted = Convert_KE_Spectrometer_to_KE_ff(evt, KE_beam_inst, suffix, "ElasFitted");
  double KE_ff_reco_AllFitted = Convert_KE_Spectrometer_to_KE_ff(evt, KE_beam_inst, suffix, "AllFitted");
  double KE_ff_reco_Data = Convert_KE_Spectrometer_to_KE_ff(evt, KE_beam_inst, suffix, "Data");
  double KE_ff_reco = 0.;
  if(evt.MC) KE_ff_reco = KE_ff_reco_AllTrue;
  else KE_ff_reco = KE_ff_reco_Data;

  TString particle_type_str = Form("%d", hadana.pitype);
  if(suffix.Contains("proton")) particle_type_str = Form("%d", hadana.ptype);

  // ==== Beam
  // == Reco
  Hist.JSFillHist(suffix, "htrack_P_beam_inst_" + suffix + "_" + particle_type_str, P_beam_inst, weight, 5000., 0., 5000.);
  Hist.JSFillHist(suffix, "htrack_KE_ff_reco_" + suffix + "_" + particle_type_str, KE_ff_reco, weight, 5000., 0., 5000.);
  Hist.JSFillHist(suffix, "htrack_PandoraSlice_" + suffix + "_" + particle_type_str, hadana.PassPandoraSliceCut(evt), weight, 2., 0., 2.);
  Hist.JSFillHist(suffix, "htrack_CaloSize_" + suffix + "_" + particle_type_str, hadana.PassCaloSizeCut(evt),weight, 2., 0.,2.);
  Hist.JSFillHist(suffix, "htrack_beam_dx_" + suffix + "_" + particle_type_str, hadana.beam_dx, weight, 200., -10., 10.);
  Hist.JSFillHist(suffix, "htrack_beam_dy_" + suffix + "_" + particle_type_str, hadana.beam_dy, weight, 200., -10.,10.);
  Hist.JSFillHist(suffix, "htrack_beam_dz_" + suffix + "_" + particle_type_str, hadana.beam_dz, weight, 200., -10.,10.);
  Hist.JSFillHist(suffix, "htrack_beam_dxy_" + suffix + "_" + particle_type_str, hadana.beam_dxy, weight, 200., -10.,10.);
  Hist.JSFillHist(suffix, "htrack_beam_costh_" + suffix + "_" + particle_type_str, hadana.beam_costh, weight, 200., -1.,1.);
  Hist.JSFillHist(suffix, "htrack_daughter_michel_score_" + suffix + "_" + particle_type_str, hadana.daughter_michel_score, weight, 100., 0.,1.);
  Hist.JSFillHist(suffix, "htrack_chi2_proton_" + suffix + "_" + particle_type_str, hadana.chi2_proton, weight, 1000., 0., 100.);
  Hist.JSFillHist(suffix, "htrack_beam_alt_len_" + suffix + "_" + particle_type_str, evt.reco_beam_alt_len, weight, 600., 0., 600.);
  Hist.JSFillHist(suffix, "htrack_beam_inst_XY_" + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.); 
  Hist.JSFillHist(suffix, "htrack_beam_start_XY_" + suffix + "_" + particle_type_str, evt.reco_beam_startX, evt.reco_beam_startY, weight, 80., -50., -10., 80., 400., 440.);
  Hist.JSFillHist(suffix, "htrack_beam_start_Z_" + suffix + "_" + particle_type_str, evt.reco_beam_startZ, weight, 1000., -100., 900.);
  Hist.JSFillHist(suffix, "htrack_beam_end_Z_" + suffix + "_" + particle_type_str, evt.reco_beam_endZ, weight, 1000., -100., 900.);
  
  int hundered_low = KE_beam_inst / 100.;
  int KE_beam_inst_50bin = KE_beam_inst / 50.;
  TString KE_range_str_50MeV = Form("KE_beam_inst%dto%dMeV_", 50 * KE_beam_inst_50bin, 50 * (KE_beam_inst_50bin  + 1));
  TString KE_range_str = Form("KE_beam_inst%dto%dMeV_", 100 * hundered_low, 100 * (hundered_low + 1));

  Hist.JSFillHist(suffix, "htrack_beam_inst_XY_" + KE_range_str + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.);

  // == Truth
  if(evt.MC){
    double P_beam_true = evt.true_beam_startP * 1000.;
    double KE_beam_true = sqrt(pow(P_beam_true, 2) + pow(mass_beam, 2)) - mass_beam;
    double true_KE_ff = hadana.true_ffKE;
    double true_P_ff = hadana.map_BB[abs(evt.true_beam_PDG)] -> KEtoMomentum(true_KE_ff);

    double KE_diff_RecoBeam_TrueFF = KE_beam_inst - true_KE_ff;
    double KE_diff_TrueBeam_TrueFF = KE_beam_true - true_KE_ff;
    double KE_diff_RecoBeam_TrueBeam = KE_beam_inst - KE_beam_true;

    double KE_diff = KE_beam_inst - true_KE_ff;
    double KE_diff_true = KE_beam_true - true_KE_ff;
    double KE_diff_beam_true = KE_beam_inst - KE_beam_true;
    TString true_beam_PDG_str = Form("%d", abs(evt.true_beam_PDG));
    double P_beam_res = (P_beam_true - P_beam_inst) / P_beam_true;

    double diff_KE_ff_reco_ElasTrue_KE_ff_true = KE_ff_reco_ElasTrue - true_KE_ff;
    double diff_KE_ff_reco_AllTrue_KE_ff_true = KE_ff_reco_AllTrue - true_KE_ff;
    double diff_KE_ff_reco_ElasFitted_KE_ff_true = KE_ff_reco_ElasFitted - true_KE_ff;
    double diff_KE_ff_reco_AllFitted_KE_ff_true = KE_ff_reco_AllFitted - true_KE_ff;

    Hist.JSFillHist(suffix, "htrack_P_beam_Res_" + suffix + "_" + particle_type_str, P_beam_res, weight, 4000, -2, 2);
    Hist.JSFillHist(suffix, "htrack_P_beam_true_" + suffix + "_" + particle_type_str, P_beam_true, weight, 5000, 0, 5000);
    Hist.JSFillHist(suffix, "htrack_P_beam_true_" + suffix + "_PID" + true_beam_PDG_str, P_beam_true, weight, 5000, 0, 5000);
    Hist.JSFillHist(suffix, "htrack_KEffTruth_" + suffix + "_" + particle_type_str, true_KE_ff, weight, 5000, 0, 5000);
    Hist.JSFillHist(suffix, "htrack_PffTruth_" + suffix + "_" + particle_type_str, true_P_ff, weight, 5000, 0, 5000);
    Hist.JSFillHist(suffix, "htrack_KE_beam_inst_vs_KEffTruth_" + suffix + "_" + particle_type_str, KE_beam_inst, true_KE_ff, weight, 300., 0., 1500., 300., 0., 1500.);
    Hist.JSFillHist(suffix, "htrack_KE_beam_inst_vs_KEffTruth_" + suffix + "_PID" + true_beam_PDG_str, KE_beam_inst, true_KE_ff, weight, 300., 0., 1500., 300., 0., 1500.);
    Hist.JSFillHist(suffix, "htrack_KE_beam_inst_vs_KE_diff_" + suffix + "_" + particle_type_str, KE_beam_inst, KE_diff, weight, 300., 0., 1500., 500., -500., 500.);
    Hist.JSFillHist(suffix, "htrack_KE_beam_inst_vs_KE_diff_" + suffix + "_PID" + true_beam_PDG_str, KE_beam_inst, KE_diff, weight, 300., 0., 1500., 500., -500., 500.);

    int hundered_low = KE_beam_inst / 100.;
    int KE_beam_inst_50bin = KE_beam_inst / 50.;

    TString KE_range_str_50MeV = Form("KE_beam_inst%dto%dMeV_", 50 * KE_beam_inst_50bin, 50 * (KE_beam_inst_50bin  + 1));
    TString KE_range_str = Form("KE_beam_inst%dto%dMeV_", 100 * hundered_low, 100 * (hundered_low + 1));
    Hist.JSFillHist(suffix, "htrack_KE_diff_RecoBeam_TrueFF_" + KE_range_str + suffix + "_"    + particle_type_str, KE_diff_RecoBeam_TrueFF, weight, 500., -500., 500.);
    Hist.JSFillHist(suffix, "htrack_KE_diff_RecoBeam_TrueFF_" + KE_range_str + suffix + "_PID" + true_beam_PDG_str, KE_diff_RecoBeam_TrueFF, weight, 500., -500., 500.);
    Hist.JSFillHist(suffix, "htrack_KE_diff_TrueBeam_TrueFF_" + KE_range_str + suffix + "_"    + particle_type_str, KE_diff_TrueBeam_TrueFF, weight, 2000., -100., 100.);
    Hist.JSFillHist(suffix, "htrack_KE_diff_TrueBeam_TrueFF_" + KE_range_str + suffix + "_PID" + true_beam_PDG_str, KE_diff_TrueBeam_TrueFF, weight, 2000., -100., 100.);

    // == These cuts should be modified for protons
    bool is_scraper = Pass_KE_diff_3sigma_Cut(KE_beam_inst, KE_diff_RecoBeam_TrueFF, suffix, "True");
    
    Hist.JSFillHist(suffix, "htrack_beam_inst_XY_" + KE_range_str_50MeV + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.);
    Hist.JSFillHist(suffix, "htrack_beam_inst_XY_" + KE_range_str       + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.);
    if(is_scraper){
      Hist.JSFillHist(suffix, "htrack_beam_inst_XY_True_scraper_" + KE_range_str       + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.);
    }
    else{
      Hist.JSFillHist(suffix, "htrack_beam_inst_XY_True_nonscraper_" + KE_range_str       + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.);;
    }

    bool pass_beam_scraper_cut = Pass_Beam_Scraper_Cut(evt.beam_inst_X, evt.beam_inst_Y, suffix, "mc");
    if(pass_beam_scraper_cut){
      Hist.JSFillHist(suffix, "htrack_KE_diff_RecoBeam_TrueFF_nonscraper_"   + KE_range_str       + suffix + "_" + particle_type_str, KE_diff_RecoBeam_TrueFF,   weight, 500., -500., 500.);
      Hist.JSFillHist(suffix, "htrack_KE_diff_RecoBeam_TrueFF_nonscraper_"   + KE_range_str_50MeV + suffix + "_" + particle_type_str, KE_diff_RecoBeam_TrueFF,   weight, 500., -500., 500.);
      Hist.JSFillHist(suffix, "htrack_KE_diff_TrueBeam_TrueFF_nonscraper_"   + KE_range_str       + suffix + "_" + particle_type_str, KE_diff_TrueBeam_TrueFF,   weight, 2000., -100., 100.);
      Hist.JSFillHist(suffix, "htrack_KE_diff_TrueBeam_TrueFF_nonscraper_"   + KE_range_str_50MeV + suffix + "_" + particle_type_str, KE_diff_TrueBeam_TrueFF,   weight, 2000., -100., 100.);
      Hist.JSFillHist(suffix, "htrack_KE_diff_RecoBeam_TrueBeam_nonscraper_" + KE_range_str       + suffix + "_" + particle_type_str, KE_diff_RecoBeam_TrueBeam, weight, 500., -500., 500.);
      Hist.JSFillHist(suffix, "htrack_KE_diff_RecoBeam_TrueBeam_nonscraper_" + KE_range_str_50MeV + suffix + "_" + particle_type_str, KE_diff_RecoBeam_TrueBeam, weight, 500., -500., 500.);
    }
  }

  // == Utilize Hypothetical Track Length Fitting Method
  if(suffix.Contains("proton") && KE_fit_gaussian > 0.){
    double KE_diff_RecoBeam_FittedFF = KE_beam_inst - KE_fit_gaussian;
    Hist.JSFillHist(suffix, "htrack_KE_diff_RecoBeam_FittedFF_" + KE_range_str       + suffix + "_" + particle_type_str, KE_diff_RecoBeam_FittedFF, weight, 500., -500., 500.);
    Hist.JSFillHist(suffix, "htrack_KE_diff_RecoBeam_FittedFF_" + KE_range_str_50MeV + suffix + "_" + particle_type_str, KE_diff_RecoBeam_FittedFF, weight, 500., -500., 500.);
    bool pass_beam_scraper_cut = Pass_Beam_Scraper_Cut(evt.beam_inst_X, evt.beam_inst_Y, suffix, "data");
    if(pass_beam_scraper_cut){
      Hist.JSFillHist(suffix, "htrack_KE_diff_RecoBeam_FittedFF_nonscraper_" + KE_range_str       + suffix + "_" + particle_type_str, KE_diff_RecoBeam_FittedFF, weight, 500., -500., 500.);
      Hist.JSFillHist(suffix, "htrack_KE_diff_RecoBeam_FittedFF_nonscraper_" + KE_range_str_50MeV + suffix + "_" + particle_type_str, KE_diff_RecoBeam_FittedFF, weight, 500., -500., 500.);      
    }

    bool is_scraper_FittedData = Pass_KE_diff_3sigma_Cut(KE_beam_inst, KE_diff_RecoBeam_FittedFF, suffix, "FittedData");
    bool is_scraper_FittedFakeData = Pass_KE_diff_3sigma_Cut(KE_beam_inst, KE_diff_RecoBeam_FittedFF, suffix, "FittedFakeData");
    if(is_scraper_FittedData){
      Hist.JSFillHist(suffix, "htrack_beam_inst_XY_FittedData_scraper_" + KE_range_str + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.);
    }
    else{
      Hist.JSFillHist(suffix, "htrack_beam_inst_XY_FittedData_nonscraper_" + KE_range_str + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.);;
    }

    if(is_scraper_FittedFakeData){
      Hist.JSFillHist(suffix, "htrack_beam_inst_XY_FittedFakeData_scraper_" + KE_range_str + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.);
    }
    else{
      Hist.JSFillHist(suffix, "htrack_beam_inst_XY_FittedFakeData_nonscraper_" + KE_range_str + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.);;
    }

    if(evt.MC){
      TString true_beam_PDG_str = Form("%d", abs(evt.true_beam_PDG));

      double P_beam_true = evt.true_beam_startP * 1000.;
      double KE_beam_true = sqrt(pow(P_beam_true, 2) + pow(mass_beam, 2)) - mass_beam;
      double true_KE_ff = hadana.true_ffKE;

      double KE_diff_TrueBeam_FittedFF = KE_beam_true - KE_fit_gaussian;
      double KE_diff_TrueFF_FittedFF = true_KE_ff - KE_fit_gaussian;

      bool is_scraper_FittedElas = Pass_KE_diff_3sigma_Cut(KE_beam_inst, KE_diff_RecoBeam_FittedFF, suffix, "FittedElas");

      Hist.JSFillHist(suffix, "htrack_KE_diff_TrueBeam_FittedFF_" + KE_range_str       + suffix + "_"    + particle_type_str, KE_diff_TrueBeam_FittedFF, weight, 500., -500., 500.);
      Hist.JSFillHist(suffix, "htrack_KE_diff_TrueBeam_FittedFF_" + KE_range_str_50MeV + suffix + "_"    + particle_type_str, KE_diff_TrueBeam_FittedFF, weight, 500., -500., 500.);
      Hist.JSFillHist(suffix, "htrack_KE_diff_TrueBeam_FittedFF_" + KE_range_str       + suffix + "_PID" + true_beam_PDG_str, KE_diff_TrueBeam_FittedFF, weight, 500., -500., 500.);
      Hist.JSFillHist(suffix, "htrack_KE_diff_TrueFF_FittedFF_"   + KE_range_str       + suffix + "_"    + particle_type_str, KE_diff_TrueFF_FittedFF  , weight, 500., -500., 500.);
      Hist.JSFillHist(suffix, "htrack_KE_diff_TrueFF_FittedFF_"   + KE_range_str_50MeV + suffix + "_"    + particle_type_str, KE_diff_TrueFF_FittedFF  , weight, 500., -500., 500.);
      Hist.JSFillHist(suffix, "htrack_KE_diff_TrueFF_FittedFF_"   + KE_range_str       + suffix + "_PID" + true_beam_PDG_str, KE_diff_TrueFF_FittedFF  , weight, 500., -500., 500.);

      if(is_scraper_FittedElas){
	Hist.JSFillHist(suffix, "htrack_beam_inst_XY_FittedElas_scraper_" + KE_range_str + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.);
      }
      else{
	Hist.JSFillHist(suffix, "htrack_beam_inst_XY_FittedElas_nonscraper_" + KE_range_str + suffix + "_" + particle_type_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 80., -50., -10., 80., 400., 440.);;
      }
    }
  }
}


void Beam_Study::Run_Beam(const anavar & evt, double weight, TString suffix, int beam_PID){
  
  if(!Pass_Beam_PID(evt, beam_PID)) return;

  TString beam_particle_str = "";
  if(beam_PID == 2212) beam_particle_str = "proton";
  else if(beam_PID = 211) beam_particle_str = "pion";
  else return;

  if(!hadana.PassPandoraSliceCut(evt)) return;
  if(!hadana.PassCaloSizeCut(evt)) return;

  // == Measure energy of proton using Hyp. Fit.
  double length_gaus = -1.;
  if(beam_PID == 2212){
    int N_hits_X = (*evt.reco_beam_calibrated_dEdX_SCE).size();
    int N_hits_X_range = (*evt.reco_beam_resRange_SCE).size();
    vector<double> this_dEdx_vec;
    vector<double> this_range_vec;
    for(unsigned int i_hit = N_hits_X - 1; i_hit > 0; i_hit--){
      this_dEdx_vec.push_back((*evt.reco_beam_calibrated_dEdX_SCE).at(i_hit));
      this_range_vec.push_back((*evt.reco_beam_resRange_SCE).at(i_hit));
    }
    length_gaus = hadana.Fit_dEdx_Residual_Length(evt, this_dEdx_vec, this_range_vec, beam_PID, false, false);
  }  
  double KE_fit_gaussian = -9999.;
  if(length_gaus > 0.) KE_fit_gaussian = hadana.map_BB[beam_PID] -> KEFromRangeSpline(length_gaus);

  // == Not empty
  FillHistBeam(evt, weight, beam_particle_str + "_NotEmpty_" + suffix, KE_fit_gaussian);

  // == Beam Quality cut
  if(!hadana.PassBeamQualityCut()) return;
  FillHistBeam(evt, weight, beam_particle_str + "_BeamQuality_" + suffix, KE_fit_gaussian);

  if(beam_PID == 2212){
    // == chi2 cut : select elastic-scattering protons
    if(hadana.chi2_proton > 10.) return;
    FillHistBeam(evt, weight, beam_particle_str + "_Elastic_" + suffix, KE_fit_gaussian);
  }
  if(beam_PID == 211){
    // == APA3 cut
    if(!hadana.PassAPA3Cut(evt)) return;
    FillHistBeam(evt, weight, beam_particle_str + "_APA3_" + suffix, KE_fit_gaussian);

    // == Michel score cut
    if(!hadana.PassMichelScoreCut()) return;
    FillHistBeam(evt, weight, beam_particle_str + "_MichelScore_" + suffix, KE_fit_gaussian);

    // == Proton veto cut
    if(!hadana.PassProtonCut()) return;
    FillHistBeam(evt, weight, beam_particle_str + "_ProtonVeto_" + suffix, KE_fit_gaussian);
  }

  // == Beam window cut
  double P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  if(P_beam_inst < P_beam_inst_cut_lower || P_beam_inst > P_beam_inst_cut_upper) return;
  FillHistBeam(evt, weight, beam_particle_str + "_BeamWindow_" + suffix, KE_fit_gaussian);

  // == End of selection
  return;
}

void Beam_Study::SaveHistograms(){
  Hist.WriteHist();
}

void Beam_Study::Run(anavar & evt, Long64_t nentries=-1){

  BookHistograms();

  // == Beam Momentum Setting
  P_beam_str = "1.0";
  P_beam_inst_cut_upper = 9999.;
  P_beam_inst_cut_lower = 0.;
  if(P_beam_str == "2.0"){
    P_beam_inst_cut_upper = 2400.;
    P_beam_inst_cut_lower = 1600.;
  }
  else if(P_beam_str == "1.0"){
    P_beam_inst_cut_upper = 1200.;
    P_beam_inst_cut_lower = 800.;
  }
  else if(P_beam_str == "0.5"){
    P_beam_inst_cut_upper = 600.;
    P_beam_inst_cut_lower = 400.;
  }
  else{
    cout << "[Beam_Study::Run] Not Valid Beam Momentum" << endl;
    return;
  }
  cout << "[Beam_Study::Run] Beam Momentum is " << P_beam_str << " GeV" << endl;


  // == Run
  if (nentries == -1) nentries = evt.fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%1000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;

    // if (Cut(ientry) < 0) continue;
    //std::cout<<evt.run<<" "<<evt.event<<" "<<evt.MC<<" "<<evt.reco_beam_true_byE_matched<<" "<<evt.true_beam_PDG<<" "<<(*evt.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    Hist.FillHist("Cutflow", 0.5, 1., 10., 0., 10.);
    ProcessEvent(evt);
    if(evt.MC){
      // == Beam P scale
      if(P_beam_str == "2.0") P_beam_inst_scale = 2.0;
      if(P_beam_str == "1.0") P_beam_inst_scale = 1.0;
      if(P_beam_str == "0.5") P_beam_inst_scale = 0.5;
    }

    Run_Beam(evt, 1., "noweight", 2212);
    Run_Beam(evt, 1., "noweight", 211);
  }

  SaveHistograms();
}

