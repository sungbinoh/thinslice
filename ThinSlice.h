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

  double true_interactions[nthinslices];
  double true_incidents[nthinslices];

  void BookHistograms();
  void ProcessEvent(const HadAna & evt);

};

#endif
