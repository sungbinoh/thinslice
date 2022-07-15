#ifndef PIONXSEC_H
#define PIONXSEC_H

#include "TFile.h"
#include "HadAna.h"

class anavar;
class HadAna;
class TH1D;
class TH2D;

class PionXsec {

 public:

  PionXsec();

  HadAna hadana;
  
  std::string fOutputFileName;
  void SetOutputFileName(std::string name){fOutputFileName = name;};
  void BookHistograms();
  void FillHistograms(const anavar & evt);
  void SaveHistograms();
  void ProcessEvent(const anavar & evt);
  void Run(anavar & evt, Long64_t nentries);

 private:
  
  TFile *outputFile;
  // Save beam KE and track length
  TH1D *htrack_BeamKE[pi::nIntTypes+1];
  TH1D *htrack_BeamP[pi::nIntTypes+1];
  TH1D *htrack_fittedP[pi::nIntTypes+1];
  TH1D *htrack_fitted_dP[pi::nIntTypes+1];
  TH1D *htrack_fittedKE[pi::nIntTypes+1];
  TH1D *htrack_fitted_dKE[pi::nIntTypes+1];
  TH1D *htrack_KECalo[pi::nIntTypes+1];
  TH1D *htrack_dKE_fitted_vs_KECalo[pi::nIntTypes+1];
  TH1D *htrack_KETruth[pi::nIntTypes+1];
  TH1D *htrack_dKE_fitted_vs_Truth[pi::nIntTypes+1];
  TH2D *htrack_dKE_fitted_vs_Truth_2D[pi::nIntTypes+1];
  // Proton energy at the track end point based on the beam momentum and track length
  TH1D *hend_energy[pi::nIntTypes+1];
};

#endif
