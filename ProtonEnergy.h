#ifndef PROTONENERGY_H
#define PROTONENERGY_H

#include "TFile.h"
#include "HadAna.h"

class anavar;
class HadAna;
class TH1D;
class TH2D;

class ProtonEnergy {

 public:

  ProtonEnergy();

  HadAna hadana;
  
  std::string fOutputFileName;
  void SetOutputFileName(std::string name){fOutputFileName = name;};
  void BookHistograms();
  void FillHistograms(const anavar & evt);
  void SaveHistograms();

  void ProcessEvent(const anavar & evt);

  void Make_dEdx_Range_Profile();

  void Run(anavar & evt, Long64_t nentries);

 private:
  
  TFile *outputFile;
  // ratio of the reconstructed track length and the predicted track length based on the beam momentum
  TH1D *htrack_length_ratio[p::nIntTypes+1];
  TH1D *htrack_length_ratio_fitted[p::nIntTypes+1];
  // same as above but consider a 46 MeV upstream energy loss
  TH1D *htrack_length_ratio_eloss46[p::nIntTypes+1];
  // Save beam KE and track length
  TH1D *htrack_fittedP[p::nIntTypes+1];
  TH1D *htrack_fitted_dP[p::nIntTypes+1];
  TH1D *htrack_fittedKE[p::nIntTypes+1];
  TH1D *htrack_fitted_dKE[p::nIntTypes+1];
  TH1D *htrack_BeamKE[p::nIntTypes+1];
  TH1D *htrack_BeamP[p::nIntTypes+1];
  TH1D *htrack_KECalo[p::nIntTypes+1];
  TH1D *htrack_dKE_fitted_vs_KECalo[p::nIntTypes+1];
  TH1D *htrack_KETruth[p::nIntTypes+1];
  TH1D *htrack_dKE_fitted_vs_Truth[p::nIntTypes+1];
  TH1D *htrack_dKE_KECalo_vs_Truth[p::nIntTypes+1];
  TH1D *htrack_length_fitted[p::nIntTypes+1];
  TH1D *htrack_length_reco[p::nIntTypes+1];
  TH1D *htrack_length_BeamKEtoRange[p::nIntTypes+1];
  TH2D *htrack_dKE_fitted_vs_Truth_2D[p::nIntTypes+1];
  // Scan Eloss for Elastic scattering
  TH1D *htrack_length_ratio_elss_scan[50];
  // Proton energy at the track end point based on the beam momentum and track length
  TH1D *hend_energy[p::nIntTypes+1];
};

#endif
