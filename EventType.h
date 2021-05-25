#ifndef EVENTTYPE_H
#define EVENTTYPE_H

#include "anavar.h"

const unsigned int nParTypes = 9;

const char parTypeName[nParTypes+1][100] = {"Data",
                                            "PiInel",
                                            "PiElas",
                                            "Muon",
                                            "misID:cosmic",
                                            "misID:p",
                                            "misID:pi",
                                            "misID:mu",
                                            "misID:e/#gamma",
                                            "misID:other"};



enum parType{
  kData,
  kPiInel,
  kPiElas,
  kMuon,
  kMIDcosmic,
  kMIDp,
  kMIDpi,
  kMIDmu,
  kMIDeg,
  kMIDother
};

#endif
