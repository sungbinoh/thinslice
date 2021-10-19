#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

namespace pi{
  const unsigned int nCuts = 7;

  const char cutName[nCuts][100] = {"Nocut",
                                    "PandoraSlice",
                                    "CaloSize",
                                    "BeamQuality",
                                    "MediandEdx",
                                    "MichelScore",
                                    "APA3"};
  
  enum cut{
    kNocut = 0,
    kPandoraSlice,
    kCaloSize,
    kBeamQuality,
    kMediandEdx,
    kMichelScore,
    kAPA3,
  };
}

namespace p{
  const unsigned int nCuts = 4;

  const char cutName[nCuts][100] = {"Nocut",
                                    "PandoraSlice",
                                    "CaloSize",
                                    "BeamQuality"};

  enum cut{
    kNocut = 0,
    kPandoraSlice,
    kCaloSize,
    kBeamQuality
  };
}

const double beam_startX_data = -27.911;
const double beam_startY_data = 424.364;
const double beam_startZ_data = 3.77836;
const double beam_startX_rms_data = 4.71128;
const double beam_startY_rms_data = 5.16472;
const double beam_startZ_rms_data = 1.10265;

const double beam_startX_mc = -30.8075;
const double beam_startY_mc = 422.41;
const double beam_startZ_mc = 0.11171;
const double beam_startX_rms_mc = 5.01719;
const double beam_startY_rms_mc = 4.50862;
const double beam_startZ_rms_mc = 0.217733;

const double beam_angleX_data = 100.454;
const double beam_angleY_data = 103.523;
const double beam_angleZ_data = 17.8288;

const double beam_angleX_mc = 101.578;
const double beam_angleY_mc = 101.189;
const double beam_angleZ_mc = 16.5942;

#endif
