#ifndef UTIL_H
#define UTIL_H

#include "TH1D.h"
#include "TH2D.h"
#include "EventType.h"

void FillHistVec1D(TH1D *hist[nParTypes+1], const double &value, const int &partype);

void FillHistVec2D(TH2D *hist[nParTypes+1], const double &value1, const double &value2, const int &partype);

#endif
