#ifndef THINSLICE_H
#define THINSLICE_H

#include "SliceParams.h"
#include "HadAna.h"

class ThinSlice {

 public:

  int reco_sliceID;
  int true_sliceID;

  TH1D *reco_incE[nslices];
  TH1D *true_incE[nslices];

  TH1D *h_truesliceid_pion_all;
  TH1D *h_truesliceid_pion_cuts;
  TH1D *h_truesliceid_pioninelastic_all;
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

  void BookHistograms();
  void FillHistograms(int cut, const HadAna & evt);
  void ProcessEvent(const HadAna & evt);
  void CalcXS();

};

#endif
