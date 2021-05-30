#ifndef EVENTTYPE_H
#define EVENTTYPE_H

#include "anavar.h"

const unsigned int nIntTypes = 9;

const char intTypeName[nIntTypes+1][100] = {"Data",
                                            "PiInel",
                                            "PiElas",
                                            "Muon",
                                            "misID:cosmic",
                                            "misID:p",
                                            "misID:pi",
                                            "misID:mu",
                                            "misID:e/#gamma",
                                            "misID:other"};



enum intType{
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
