#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

#include "anavar.h"

const unsigned int nCuts = 6;

const char cutName[nCuts][100] = {"Nocut",
                                  "PandoraSlice",
                                  "BeamQuality",
                                  "APA3",
                                  "CaloSize",
                                  "TrunMeandEdx"};

enum cut{
  kNocut = 0,
  kPandoraSlice,
  kBeamQuality,
  kAPA3,
  kCaloSize,
  kTrunMeandEdx
};

#endif
