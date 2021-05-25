#ifndef THINSLICE_H
#define THINSLICE_H

#include "TFile.h"
#include "SliceParams.h"
#include "HadAna.h"
#include "Unfold.h"

class ThinSlice {

 public:

  int reco_sliceID;
  int true_sliceID;

  bool isTestSample;

  TH1D *reco_incE[nthinslices];
  TH1D *true_incE[nthinslices];
  TH1D *reco_AngCorr;
  TH1D *true_AngCorr;

  TH1D *h_truesliceid_pion_all;
  TH1D *h_truesliceid_pion_uf;
  TH1D *h_truesliceid_pion_cuts;
  TH1D *h_truesliceid_pioninelastic_all;
  TH1D *h_truesliceid_pioninelastic_uf;
  TH1D *h_truesliceid_pioninelastic_cuts;
  TH1D *h_recosliceid_allevts_cuts;
  TH1D *h_recosliceid_pion_cuts;
  TH1D *h_recosliceid_pioninelastic_cuts;

  double true_interactions[nthinslices];
  double true_incidents[nthinslices];

  TH1D *htrue_beam_endZ[nCuts][nParTypes+1];
  TH1D *hreco_beam_endZ[nCuts][nParTypes+1];
  TH1D *hreco_true_beam_endZ[nCuts][nParTypes+1];
  TH2D *hreco_vs_true_beam_endZ[nCuts][nParTypes+1];
  TH2D *hreco_true_vs_true_beam_endZ[nCuts][nParTypes+1];

  TH1D *htrue_beam_endZ_SCE[nCuts][nParTypes+1];
  TH1D *hreco_beam_endZ_SCE[nCuts][nParTypes+1];
  TH1D *hreco_true_beam_endZ_SCE[nCuts][nParTypes+1];
  TH2D *hreco_vs_true_beam_endZ_SCE[nCuts][nParTypes+1];
  TH2D *hreco_true_vs_true_beam_endZ_SCE[nCuts][nParTypes+1];

  TH1D *htrue_sliceID[nCuts][nParTypes+1];
  TH1D *hreco_sliceID[nCuts][nParTypes+1];
  TH1D *hreco_true_sliceID[nCuts][nParTypes+1];
  TH2D *hreco_vs_true_sliceID[nCuts][nParTypes+1];
  TH2D *hreco_true_vs_true_sliceID[nCuts][nParTypes+1];

  TH1D *hmediandEdx[nCuts][nParTypes+1];
  TH1D *hdaughter_michel_score[nCuts][nParTypes+1];

  TH1D *hdeltax[nCuts][nParTypes+1];
  TH1D *hdeltay[nCuts][nParTypes+1];
  TH1D *hdeltaz[nCuts][nParTypes+1];
  TH1D *hcostheta[nCuts][nParTypes+1];

  std::string fOutputFileName;
  TFile *outputFile;
  void SetOutputFileName(std::string name){fOutputFileName = name;};
  void BookHistograms();
  void FillHistograms(int cut, const HadAna & evt);
  void SaveHistograms();

  void ProcessEvent(const HadAna & evt, Unfold & uf);
  void CalcXS(const Unfold & uf);

  void Run(HadAna & evt, Unfold & uf);

};

#endif
