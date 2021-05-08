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

  void BookHistograms();
  void ProcessEvent(const HadAna & evt);
  void CalcXS();

};

#endif
