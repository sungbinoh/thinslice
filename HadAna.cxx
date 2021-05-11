#include "HadAna.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "util.h"
#include <iostream>

bool HadAna::isTrueSelectedPart() const{
  for (size_t i = 0; i<truepdglist.size(); ++i){
    if (true_beam_PDG == truepdglist[i]) return true;
  }
  return false;
}

int HadAna::GetParType(){

  if (!reco_beam_true_byE_matched){
    return kMisID;
  }
  else if (true_beam_PDG == -13){
    return kPrimMuP;
  }
  else if (true_beam_PDG == 13){
    return kPrimMuM;
  }
  else if (true_beam_PDG == 211){
    if ((*true_beam_endProcess) == "pi+Inelastic"){
      return kPrimPiPInEl;
    }
    else return kPrimPiPEl;
  }
  else if (true_beam_PDG == 2212){
    if ((*true_beam_endProcess) == "protonInelastic"){
      return kPrimProInEl;
    }
    else return kPrimProEl;
  }
  
  return kOther;
}


bool HadAna::PassPandoraSliceCut() const{
  return reco_beam_type == pandora_slice_pdg;
}

bool HadAna::PassBeamQualityCut() const{

  double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
  double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

  double projectX = (true_beam_startX + -1*true_beam_startZ*(true_beam_startDirX/true_beam_startDirZ) );
  double projectY = (true_beam_startY + -1*true_beam_startZ*(true_beam_startDirY/true_beam_startDirZ) );
  double cos = true_beam_startDirX*reco_beam_trackDirX + true_beam_startDirY*reco_beam_trackDirY + true_beam_startDirZ*reco_beam_trackDirZ;
  
  if ( (reco_beam_startX - projectX) < xlow )
    return false;
  
  if ( (reco_beam_startX - projectX) > xhigh )
    return false;
  
  if ( (reco_beam_startY - projectY) < ylow )
    return false;
  
  if ( (reco_beam_startY - projectY) > yhigh )
    return false;
  
  if (reco_beam_startZ < zlow || zhigh < reco_beam_startZ)
    return false;
  
  if ( cos < coslow)
    return false;
  
  return true;
  
};

bool HadAna::PassAPA3Cut() const{

  double cutAPA3_Z = 226.;

  return reco_beam_endZ < cutAPA3_Z;
}

bool HadAna::PassCaloSizeCut() const{
  
  return !(reco_beam_calo_wire->empty());
}

bool HadAna::PassMichelScoreCut() const{
  
  return daughter_michel_score < 0.55;
}

bool HadAna::PassMediandEdxCut() const{

  return median_dEdx < 2.4;
}

bool HadAna::PassAllCuts() const{
  return PassPandoraSliceCut()&&
    PassBeamQualityCut()&&
    PassAPA3Cut()&&
    PassCaloSizeCut()&&
    PassMichelScoreCut()&&
    PassMediandEdxCut();
}

void HadAna::ProcessEvent(){

  partype = -1;

  median_dEdx = -1;
  daughter_michel_score = -1;
  if (!reco_beam_calo_wire->empty()){

    median_dEdx = TMath::Median(reco_beam_calibrated_dEdX_SCE->size(), &(*reco_beam_calibrated_dEdX_SCE)[0]);
    if (!reco_daughter_PFP_michelScore_collection->empty()){
      daughter_michel_score = TMath::MaxElement(reco_daughter_PFP_michelScore_collection->size(), &(*reco_daughter_PFP_michelScore_collection)[0]);
      //std::cout<<daughter_michel_score<<std::endl;
    }
    else{
      daughter_michel_score = 0;
    }
  }

  if (MC){
    partype = GetParType();
  }
}

