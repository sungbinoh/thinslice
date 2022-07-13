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
    htrack_length_reco[i] = new TH1D(Form("htrack_length_reco_%d",i),Form("%s;track_length", p::intTypeName[i]), 100, 0, 1000);
    htrack_length_BeamKEtoRange[i] = new TH1D(Form("htrack_length_BeamKEtoRange_%d",i),Form("%s;track_length", p::intTypeName[i]), 100, 0, 1000);
    htrack_length_fitted[i] = new TH1D(Form("htrack_length_fitted_%d",i),Form("%s;track_length_fitted", p::intTypeName[i]), 100, 0, 1000);
    htrack_BeamKE[i] = new TH1D(Form("htrack_BeamKE_%d",i),Form("%s;track_BeamKE", p::intTypeName[i]), 100, 0, 1000);
    htrack_BeamP[i] = new TH1D(Form("htrack_BeamP_%d",i),Form("%s;track_BeamP", p::intTypeName[i]), 500, 0, 2000);
    htrack_fittedP[i] = new TH1D(Form("htrack_fittedP_%d",i),Form("%s;track_fittedP", p::intTypeName[i]), 500, 0, 2000);
    htrack_fitted_dP[i] = new TH1D(Form("htrack_fitted_dP_%d",i),Form("%s;track_fitted_dP", p::intTypeName[i]), 2000, -1000, 1000);
    htrack_fittedKE[i] = new TH1D(Form("htrack_fittedKE_%d",i),Form("%s;track_fittedKE", p::intTypeName[i]), 500, 0, 2000);
    htrack_fitted_dKE[i] = new TH1D(Form("htrack_fitted_dKE_%d",i),Form("%s;track_fitted_dKE", p::intTypeName[i]), 2000, -1000, 1000);
    htrack_KECalo[i] = new TH1D(Form("htrack_KECalo_%d",i),Form("%s;track_KECalo", p::intTypeName[i]), 500, 0, 2000);
    htrack_dKE_fitted_vs_KECalo[i] = new TH1D(Form("htrack_dKE_fitted_vs_KECalo_%d",i),Form("%s;track_dKE_fitted_vs_KECalo", p::intTypeName[i]), 2000, -1000, 1000);
    htrack_KETruth[i] = new TH1D(Form("htrack_KETruth_%d",i),Form("%s;track_KETruth", p::intTypeName[i]), 500, 0, 2000);
    htrack_dKE_fitted_vs_Truth[i] = new TH1D(Form("htrack_dKE_fitted_vs_Truth_%d",i),Form("%s;track_dKE_fitted_vs_Truth", p::intTypeName[i]), 2000, -1000, 1000);
    htrack_length_ratio_fitted[i] = new TH1D(Form("htrack_length_ratio_fitted_%d",i),Form("%s;track_length_ratio_fitted", p::intTypeName[i]), 1000, 0, 2);
    htrack_length_ratio_eloss46[i] = new TH1D(Form("htrack_length_ratio_eloss46_%d",i),Form("E_{loss} = 46 MeV, %s;track_length_ratio", p::intTypeName[i]), 100, 0, 2);
    hend_energy[i] = new TH1D(Form("hend_energy_%d",i),Form("%s;End point energy", p::intTypeName[i]), 100, -500, 500);
  }

  // == Scanning Eloss
  for (int i = 0; i < 50; i++){
    htrack_length_ratio_elss_scan[i] = new TH1D(Form("htrack_length_ratio_eloss%dMeV", i), Form("E_{loss} = %d MeV;track_length_ratio", i), 100, 0, 2);
  }

}

void ProtonEnergy::ProcessEvent(const anavar & evt){
  hadana.ProcessEvent(evt);
}

void ProtonEnergy::FillHistograms(const anavar & evt){

  double beamKE = sqrt(pow(evt.beam_inst_P*1000, 2) + pow(938.272, 2)) - 938.272;
  double hypoth_length = hadana.map_BB[2212]->RangeFromKESpline(beamKE);
  double hypoth_length_eloss46 = hadana.map_BB[2212]->RangeFromKESpline(beamKE-46);
  double end_energy = hadana.map_BB[2212]->KEAtLength(beamKE, evt.reco_beam_alt_len);
  double KE_calo = hadana.Integrate_dEdx(*evt.reco_beam_dEdX_SCE, *evt.reco_beam_TrkPitch_SCE);
  //double KE_truth = ;
  //cout<<evt.beam_inst_P<<" "<<beamKE<<" "<<hypoth_length<<" "<<evt.reco_beam_alt_len<<endl;
  FillHistVec1D(htrack_length_ratio, evt.reco_beam_alt_len/hypoth_length, hadana.ptype);
  FillHistVec1D(htrack_length_ratio_eloss46, evt.reco_beam_alt_len/hypoth_length_eloss46, hadana.ptype);
  FillHistVec1D(htrack_length_reco, evt.reco_beam_alt_len, hadana.ptype);
  FillHistVec1D(htrack_length_BeamKEtoRange, hypoth_length, hadana.ptype);
  FillHistVec1D(htrack_BeamKE, beamKE, hadana.ptype);
  FillHistVec1D(htrack_BeamP, evt.beam_inst_P*1000., hadana.ptype);
  FillHistVec1D(htrack_KECalo, KE_calo, hadana.ptype);
  FillHistVec1D(hend_energy, end_energy, hadana.ptype);

  // == Fit hypothetical residual range
  if((*evt.reco_beam_dEdX_SCE).size() > 50){
    double fitted_length = hadana.Fit_dEdx_Residual_Length(evt, *evt.reco_beam_dEdX_SCE, *evt.reco_beam_resRange_SCE, 2212, false);
    double fitted_KE = hadana.map_BB[2212]->KEFromRangeSpline(fitted_length);
    double fitted_dKE = beamKE - fitted_KE;
    double dKE_calo_vs_fitted = KE_calo - fitted_KE;
    double fitted_P = hadana.map_BB[2212]->KEtoMomentum(fitted_KE);
    double fitted_dP = evt.beam_inst_P*1000. - fitted_P;
    double true_KE_ff = hadana.true_ffKE;
    double dKE_true_KE_ff_vs_fitted = true_KE_ff - fitted_KE;
    //cout << "[ProtonEnergy::FillHistograms] fitted_KE : " << fitted_KE << ", KE_calo : " << KE_calo << endl;

    FillHistVec1D(htrack_length_fitted,fitted_length, hadana.ptype);
    FillHistVec1D(htrack_length_ratio_fitted, evt.reco_beam_alt_len/fitted_length, hadana.ptype);
    FillHistVec1D(htrack_fittedP, fitted_P, hadana.ptype);
    FillHistVec1D(htrack_fitted_dP, fitted_dP, hadana.ptype);
    FillHistVec1D(htrack_fittedKE, fitted_KE, hadana.ptype);
    FillHistVec1D(htrack_fitted_dKE, fitted_dKE, hadana.ptype);
    FillHistVec1D(htrack_dKE_fitted_vs_KECalo, dKE_calo_vs_fitted, hadana.ptype);
    FillHistVec1D(htrack_KETruth, true_KE_ff, hadana.ptype);    
    FillHistVec1D(htrack_dKE_fitted_vs_Truth, dKE_true_KE_ff_vs_fitted, hadana.ptype);; 
   /*
    cout << "[ProtonEnergy::FillHistograms] fitted_length : " << fitted_length << endl;
    cout << "[ProtonEnergy::FillHistograms] reco length : " << evt.reco_beam_alt_len << endl;
    cout << "[ProtonEnergy::FillHistograms] hypoth_length : " << hypoth_length << endl;
    cout << "[ProtonEnergy::FillHistograms] hypoth_length_eloss46 : " << hypoth_length_eloss46 << endl;
    */
  }

  // == Scanning Eloss
  //if(hadana.ptype == p::kPElas){
  if(hadana.ptype == p::kData){
    for(int i = 0; i < 50; i++){
      double this_eloss = i + 0.;
      double hypoth_length_eloss = hadana.map_BB[2212]->RangeFromKESpline(beamKE - this_eloss);
      htrack_length_ratio_elss_scan[i] -> Fill(evt.reco_beam_alt_len/hypoth_length_eloss);
    }
  }

}

void ProtonEnergy::SaveHistograms(){
  outputFile->cd();
  outputFile->Write();
}

void ProtonEnergy::Make_dEdx_Range_Profile(){
  vector<double> dEdx_vec;
  vector<double> range_vec;
  double this_KE = 1.; // [MeV]
  double KE_step = 1.;
  double KE_max = 10000.;
  int i_max = (KE_max - this_KE) / KE_step;
  for(int i_profile = 0; i_profile < i_max; i_profile++){
    double this_dEdx = hadana.map_BB[2212]->meandEdx(this_KE);
    double this_range = hadana.map_BB[2212]->RangeFromKESpline(this_KE);
    cout << Form("[ProtonEnergy::Make_dEdx_Range_Profile] (this_range,this_dEdx) = (%f, %f)", this_range, this_dEdx) << endl;
    dEdx_vec.push_back(this_dEdx);
    range_vec.push_back(this_range);
    this_KE += KE_step;
  }

  TGraph *profile_gr = new TGraph(i_max - 1,&range_vec[0], &dEdx_vec[0]);
  profile_gr -> SetName("Proton_dEdx_range_Profile");
  profile_gr -> Write();

  dEdx_vec.clear();
  range_vec.clear();
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
  //Make_dEdx_Range_Profile();
}

