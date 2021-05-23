#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

#include "anavar.h"

const unsigned int nCuts = 7;

const char cutName[nCuts][100] = {"Nocut",
                                  "PandoraSlice",
                                  "CaloSize",
                                  "BeamQuality",
                                  "APA3",
                                  "MichelScore",
                                  "MediandEdx"};

enum cut{
  kNocut = 0,
  kPandoraSlice,
  kCaloSize,
  kBeamQuality,
  kAPA3,
  kMichelScore,
  kMediandEdx
};

#endif
