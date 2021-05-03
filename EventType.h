#ifndef EVENTTYPE_H
#define EVENTTYPE_H

#include "anavar.h"

const unsigned int nParTypes = 13;

const char parTypeName[nParTypes+1][100] = {"all",
                                            "MisID",
                                            "PrimMuP",
                                            "PrimPiPInEl",
                                            "PrimPiPEl",
                                            "PrimKaPInEl",
                                            "PrimKaPEl",
                                            "PrimProInEl",
                                            "PrimProEl",
                                            "PrimMuM",
                                            "PrimPiMInEl",
                                            "PrimPiMEl",
                                            "PrimKaMInEl",
                                            "Unknown"};

enum parType{
  kMisID = 1,
  kPrimMuP,
  kPrimPiPInEl,
  kPrimPiPEl,
  kPrimKaPInEl,
  kPrimKaPEl,
  kPrimProInEl,
  kPrimProEl,
  kPrimMuM,
  kPrimPiMInEl,
  kPrimPiMEl,
  kPrimKaMInEl,
  kUnknown
};

#endif
