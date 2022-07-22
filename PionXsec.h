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
  void FillHistograms_precut(const anavar & evt);
  void FillHistograms(const anavar & evt);
  void SaveHistograms();
  void ProcessEvent(const anavar & evt);
  void Run(anavar & evt, Long64_t nentries);

 private:
  
  TFile *outputFile;

  // ==== Histograms for normalization
  TH1D *h_cutflow;

  // ==== Histograms for beam
  // == Before PiCuts
  TH1D *htrack_BeamKE_precut[pi::nIntTypes+1];
  TH1D *htrack_BeamP_precut[pi::nIntTypes+1];
  TH1D *htrack_PandoraSlice_precut[pi::nIntTypes+1];
  TH1D *htrack_CaloSize_precut[pi::nIntTypes+1];
  TH1D *htrack_beam_dx_precut[pi::nIntTypes+1]; // calo StartX
  TH1D *htrack_beam_dy_precut[pi::nIntTypes+1];
  TH1D *htrack_beam_dz_precut[pi::nIntTypes+1];
  TH1D *htrack_beam_dxy_precut[pi::nIntTypes+1];
  TH1D *htrack_beam_costh_precut[pi::nIntTypes+1];
  TH2D *htrack_beam_inst_XY_precut[pi::nIntTypes+1]; // beam Inst XY
  TH1D *htrack_daughter_michel_score_precut[pi::nIntTypes+1];
  TH1D *htrack_chi2_proton_precut[pi::nIntTypes+1];
  TH1D *htrack_KEffTruth_precut[pi::nIntTypes+1];
  TH1D *htrack_PffTruth_precut[pi::nIntTypes+1];


  // == After PiCuts
  TH1D *htrack_BeamKE[pi::nIntTypes+1];
  TH1D *htrack_BeamP[pi::nIntTypes+1];
  TH1D *htrack_PandoraSlice[pi::nIntTypes+1];
  TH1D *htrack_CaloSize[pi::nIntTypes+1];
  TH1D *htrack_beam_dx[pi::nIntTypes+1]; // calo StartX
  TH1D *htrack_beam_dy[pi::nIntTypes+1];
  TH1D *htrack_beam_dz[pi::nIntTypes+1];
  TH1D *htrack_beam_dxy[pi::nIntTypes+1];
  TH1D *htrack_beam_costh[pi::nIntTypes+1];
  TH2D *htrack_beam_inst_XY[pi::nIntTypes+1]; // beam Inst XY
  TH1D *htrack_daughter_michel_score[pi::nIntTypes+1];
  TH1D *htrack_chi2_proton[pi::nIntTypes+1];
  TH1D *htrack_KEffTruth[pi::nIntTypes+1];
  TH1D *htrack_PffTruth[pi::nIntTypes+1];

 
  // == More plots
  TH1D *htrack_BeamP_true[pi::nIntTypes+1];
  TH1D *htrack_BeamKE_loss[pi::nIntTypes+1];
  TH1D *htrack_BeamKE_loss_true_mass[pi::nIntTypes+1];
  TH1D *htrack_Beam_alt_length[pi::nIntTypes+1];
  TH1D *htrack_fittedP[pi::nIntTypes+1];
  TH1D *htrack_fitted_dP[pi::nIntTypes+1];
  TH1D *htrack_fittedKE[pi::nIntTypes+1];
  TH1D *htrack_fitted_dKE[pi::nIntTypes+1];
  TH1D *htrack_KECalo[pi::nIntTypes+1];
  TH1D *htrack_dKE_fitted_vs_KECalo[pi::nIntTypes+1];
  TH1D *htrack_dKE_fitted_vs_Truth[pi::nIntTypes+1];
  TH1D *htrack_BeamX_true[pi::nIntTypes+1];
  TH1D *htrack_BeamY_true[pi::nIntTypes+1];
  TH1D *htrack_BeamZ_true[pi::nIntTypes+1];

  TH2D *htrack_dKE_fitted_vs_Truth_2D[pi::nIntTypes+1];
  // Proton energy at the track end point based on the beam momentum and track length
  TH1D *hend_energy[pi::nIntTypes+1];

  // ==== Histograms for daughters
  // == Truth matched
  TH1D *hdaughter_proton_trackScore;
  TH1D *hdaughter_proton_emScore;
  TH1D *hdaughter_proton_chi2_proton;
  TH1D *hdaughter_pion_trackScore;
  TH1D *hdaughter_pion_emScore;
  TH1D *hdaughter_pion_chi2_proton;
  TH1D *hdaughter_other_trackScore;
  TH1D *hdaughter_other_emScore;
  TH1D *hdaughter_other_chi2_proton;

};

#endif
